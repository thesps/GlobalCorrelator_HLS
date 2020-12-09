#include "algo.h"
#include "TreeReduce.h"
#include <cmath>
#include <cassert>
#include <hls_math.h>
#include <iostream>

#ifndef __SYNTHESIS__
#include <cstdio>
#endif

template<class data_T, int N>
inline float real_val_from_idx(unsigned i){
    // Treat the index as the top N bits
    static constexpr int NB = ceillog2(N); // number of address bits for table
    data_T x(0);
    x(x.width-1, x.width-NB) = i;
    return (float) x;
}

template<class data_T, int N>
inline unsigned idx_from_real_val(data_T x){
    // Slice the top N bits to get an index into the table
    static constexpr int NB = ceillog2(N); // number of address bits for table
    ap_uint<NB> y = x(x.width-1, x.width-NB); // slice the top N bits of input
    return (unsigned) y(NB-1, 0);
}

template<class data_T, class table_T, int N>
void init_invert_table(table_T table_out[N]){
    // The template data_T is the data type used to address the table
    for(unsigned i = 0; i < N; i++){
        float x = real_val_from_idx<data_T, N>(i);
        table_T inv_x = 1 / x;
        table_out[i] = inv_x;
    }
}

template<class in_t, class table_t, int N>
table_t invert_with_shift(in_t in){
    table_t inv_table[N];
    init_invert_table<in_t, table_t, N>(inv_table);

    // find the first '1' in the denominator
    int msb = 0;
    for(int b = 0; b < in.width; b++){
        #pragma HLS unroll
        if(in[b]) msb = b;
    }
    // shift up the denominator such that the left-most bit (msb) is '1'
    in_t in_shifted = in << (in.width-msb-1);
    // lookup the inverse of the shifted input
    int idx = idx_from_real_val<in_t,N>(in_shifted);
    table_t inv_in = inv_table[idx];
    // shift the output back
    table_t out = inv_in << (in.width-msb-1);
#ifndef __SYNTHESIS__
    std::cout << "           x " << in << ", msb = " << msb << ", shift = " << (in.width-msb) << ", idx = " << idx << std::endl;
    std::cout << "     pre 1 / " << in_shifted << " = " << inv_in << "(" << 1/(float)in_shifted << ")" << std::endl;
    std::cout << "    post 1 / " << in << " = " << out << "(" << 1/(float)in << ")" << std::endl;
#endif 
    return out;
}

void updateAxis(etaphi_t seed_eta, etaphi_t seed_phi,
                     pt_t sum_pt, pt_etaphi_t sum_pt_eta, pt_etaphi_t sum_pt_phi, count_t count, Jet & jet) {
    inv_pt_t inv_pt = invert_with_shift<pt_t, inv_pt_t, N_table_inv_pt>(sum_pt);
    
    etaphi_t jet_eta = seed_eta + etaphi_t(sum_pt_eta * inv_pt);
    etaphi_t jet_phi = seed_phi + etaphi_t(sum_pt_phi * inv_pt);
#ifndef __SYNTHESIS__
    std::cout << " uncorr eta: " << seed_eta << ", phi: " << seed_phi << std::endl;
    std::cout << "   corr eta: " << jet_eta << ", phi: " << jet_phi << std::endl;
#endif
    bool fiducial = (sum_pt >= JET_PT_CUT) &&
        (-FIDUCIAL_ETA_PHI <= jet_eta && jet_eta <= FIDUCIAL_ETA_PHI) &&
        (-FIDUCIAL_ETA_PHI <= jet_phi && jet_phi <= FIDUCIAL_ETA_PHI);
    if (fiducial) {
#ifndef __SYNTHESIS__
        //printf("HW algo         found jet pt %.2f, eta %+5.3f, phi %+5.3f\n", sum_pt * 0.25, jet_eta*0.01, jet_phi*0.01);
#endif
        jet.hwPt  = sum_pt;
        jet.hwEta = jet_eta;
        jet.hwPhi = jet_phi;
        jet.nCand = count;
        jet.iSeed = 0;
    } else {
        clear(jet);
    }

}

void addJetPtSorted(const Jet currentJet, Jet myjet[NJETS]){
    // Add one jet to the list, keeping them pt sorted
    #pragma HLS latency min=5 max=6
    #pragma HLS array_partition variable=myjet complete
    // tmp should get shifted along the list of jets
    Jet forward = currentJet;
    JetStreamSort:
    for(int i = 0; i < NJETS; i++){
        Jet tmp;
        bool GT = forward.hwPt >= myjet[i].hwPt;
        tmp = GT ? myjet[i] : forward;
        myjet[i] = GT ? forward : myjet[i];
        forward = tmp;
    }
}

void formJet(Particle work[NPARTICLES], bool incone[NPARTICLES], etaphi_t seed_eta, etaphi_t seed_phi, pt_t & sum_pt, pt_etaphi_t & sum_pt_eta, pt_etaphi_t & sum_pt_phi, count_t & count){

	#pragma HLS pipeline
#ifndef __SYNTHESIS__
    std::cout << " Seed: pt " << sum_pt << ", eta: " << seed_eta << ", phi: " << seed_phi << std::endl;
#endif
    pt_t sum_pts[NPARTICLES];
    pt_etaphi_t sum_pt_etas[NPARTICLES];
    pt_etaphi_t sum_pt_phis[NPARTICLES];
	#pragma HLS array_partition variable=sum_pts complete
	#pragma HLS array_partition variable=sum_pt_etas complete
	#pragma HLS array_partition variable=sum_pt_phis complete

    // this block builds the jet out of the seed, and marks the used candidates
    // Reset seed pT to 0, as it will be clustered from work
	sum_pt = 0; sum_pt_eta = 0; sum_pt_phi = 0;
	count = (sum_pt > 0) ? 1 : 0;
    JetsLoopParticleLoop:
    for (unsigned int i = 0; i < NPARTICLES; ++i) {
        #pragma HLS unroll
        detaphi_t deta = work[i].hwEta - seed_eta;
        detaphi_t dphi = work[i].hwPhi - seed_phi;
        // phi wrap
        detaphi_t dphi0 = dphi > PI ? (detaphi_t) (TWOPI - dphi) : (detaphi_t) dphi;
        detaphi_t dphi1 = dphi < -PI ? (detaphi_t) (TWOPI + dphi) : (detaphi_t) dphi;
        detaphi_t dphiw = dphi > 0 ? dphi0 : dphi1;
        incone[i] = deta*deta + dphiw*dphiw < R2CONE;
        pt_t maybePt =  incone[i] ? work[i].hwPt : pt_t(0);
        sum_pts[i]    = maybePt;
        sum_pt_etas[i] = maybePt * (detaphi_t)deta; // range is bounded since |deta| < 0.04 = 40 units
        sum_pt_phis[i] = maybePt * (detaphi_t)dphiw;
        count += (maybePt > 0);
#ifndef __SYNTHESIS__
        if(incone[i]){
            std::cout << " part: pt " << maybePt << ", eta: " << work[i].hwEta << ", phi: " << work[i].hwPhi << std::endl;
            std::cout << "    d: deta: " << deta << ", dphi: " << dphi << std::endl;
        }
#endif
    }
    // sum the in-cone pt and pt-weighted delta-eta, delta-phis
    Op_add<pt_t> op_add_pt;
    sum_pt = reduce<pt_t, NPARTICLES, Op_add<pt_t>>(sum_pts, op_add_pt);
    Op_add<pt_etaphi_t> op_add_etaphi;
    sum_pt_eta = reduce<pt_etaphi_t, NPARTICLES, Op_add<pt_etaphi_t>>(sum_pt_etas, op_add_etaphi);
    sum_pt_phi = reduce<pt_etaphi_t, NPARTICLES, Op_add<pt_etaphi_t>>(sum_pt_phis, op_add_etaphi);
}

void findSeed(Particle work[NPARTICLES], Particle & seedp, etaphi_t & seed_eta, etaphi_t & seed_phi){
	#pragma HLS pipeline
    // Pick the highest pT particle as the seed
    Op_max<Particle> op_max;
    seedp = reduce<Particle, NPARTICLES, Op_max<Particle>>(work, op_max);
    seed_eta = seedp.hwPt > 0 ? seedp.hwEta : etaphi_t (0);
    seed_phi = seedp.hwPt > 0 ? seedp.hwPhi : etaphi_t (0);
}

void updateWork(Particle work[NPARTICLES], bool incone[NPARTICLES]){
    for (unsigned int i = 0; i < NPARTICLES; ++i) {
        #pragma HLS unroll
    	if(incone[i]){
    		work[i].hwPt = pt_t(0);
    	}
    }
}

void copyInput(const Particle particles[NPARTICLES], Particle work[NPARTICLES]){
    CopyInLoop:
    for (unsigned int i = 0; i < NPARTICLES; ++i) {
        #pragma HLS unroll
        work[i] = particles[i];
    }
}

void clearJets(Jet myjet[NJETS]){
    ClearJets:
    for(int i = 0; i < NJETS; i++){
        #pragma HLS unroll
        clear(myjet[i]);
    }
}

void copyOutput(Jet myjet[NJETS], Jet jet[NJETS]){
    CopyOutLoop:
    for(int i = 0; i < NJETS; i++){
        #pragma HLS unroll
        jet[i] = myjet[i];
    }
}

void algo_main(const Particle particles[NPARTICLES], Jet jet[NJETS]) {
    #pragma HLS array_partition variable=particles complete
    #pragma HLS array_partition variable=jet complete
    #pragma HLS interface ap_none port=jet 
    #pragma HLS data_pack variable=particles
    #pragma HLS data_pack variable=jet

    Jet myjet[NJETS];
    #pragma HLS array_partition variable=myjet complete
    clearJets(myjet);
    
    Particle work[NPARTICLES];
    #pragma HLS array_partition variable=work complete
    copyInput(particles, work);

    bool incone [NPARTICLES];
    #pragma HLS array_partition variable=incone complete

    JetsLoop:
    for (unsigned int j = 0; j < NJETS; ++j) {
        #pragma HLS pipeline
        Jet currentJet;
        Particle seedp;
        etaphi_t seed_eta, seed_phi;
        pt_t sum_pt;
        pt_etaphi_t sum_pt_eta, sum_pt_phi;
        count_t count;
        findSeed(work, seedp, seed_eta, seed_phi);
        formJet(work, incone, seed_eta, seed_phi, sum_pt, sum_pt_eta, sum_pt_phi, count);
        updateAxis(seed_eta, seed_phi, sum_pt, sum_pt_eta, sum_pt_phi, count, currentJet);
        updateWork(work, incone);
        addJetPtSorted(currentJet, myjet);
    }

    copyOutput(myjet, jet);
}

