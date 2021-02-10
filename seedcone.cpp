#include "algo.h"
#include "TreeReduce.h"
#include <cmath>
#include <cassert>
#include <hls_math.h>
#include <iostream>

#ifndef __SYNTHESIS__
#include <cstdio>
#include <vector>
#endif

template<class data_T, int N>
inline float real_val_from_idx(unsigned i){
    // Treat the index as the top N bits
    static constexpr int NB = ceillog2(N); // number of address bits for table
    data_T x(0);
    // The MSB of 1 is implicit in the table
    x[x.width-1] = 1;
    // So we can use the next NB bits for real data
    x(x.width-2, x.width-NB-1) = i;
    return (float) x;
}

template<class data_T, int N>
inline unsigned idx_from_real_val(data_T x){
    // Slice the top N bits to get an index into the table
    static constexpr int NB = ceillog2(N); // number of address bits for table
    // Slice the top-1 NB bits of the value
    // the MSB of '1' is implicit, so only slice below that
    ap_uint<NB> y = x(x.width-2, x.width-NB-1);
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
    std::cout << " sum_pt_eta: " << sum_pt_eta << ", sum_pt_eta * 1 / pt: " << etaphi_t(sum_pt_eta * inv_pt) << std::endl;
    std::cout << " sum_pt_phi: " << sum_pt_phi << ", sum_pt_phi * 1 / pt: " << etaphi_t(sum_pt_phi * inv_pt) << std::endl;
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

// from https://github.com/Licheng-Guo/vivado-hls-broadcast-optimization/
template<class T>
T dreg(T in){
    #pragma HLS pipeline
    #pragma HLS inline off
    #pragma HLS interface port=return register
    return in;
}

/*template<class T>
void broadcast128(T in, T out[128]){
    #pragma HLS pipeline
    T tmp[16];
    for(int i = 0; i < 16; i++){
        #pragma HLS unroll
        tmp[i] = dreg<T>(in);
    }
    for(int i = 0; i < 128; i++){
        #pragma HLS unroll
        out[i] = dreg<T>(tmp[i/8]);
    }
}*/

/*template<class T, int NOUT, int MAXFANOUT>
void broadcast(T in, T out[NOUT]){
    #pragma HLS pipeline
    static constexpr int NSTAGES = ceillog2(MAXFANOUT);
    // Assigning NOUT for each stages is possibly wasteful
    // but we won't write to all, so some should be trimmed
    T tmp[NSTAGES+1][NOUT]; 
    tmp[0][0] = in;
    int lim = 1;
    for(int i = 1; i <= NSTAGES; i++){
        #pragma HLS unroll
        lim *= MAXFANOUT;
        for(int j = 0; j < lim; j++){
            #pragma HLS unroll
            tmp[i][j] = dreg<T>(tmp[i-1][j/MAXFANOUT]);
        }
    }
    for(int i = 0; i < NOUT; i++){
        #pragma HLS unroll
        out[i] = tmp[NSTAGES+1][i];
    }
}*/

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

void findJet(const Particle work[NPARTICLES], bool incone[NPARTICLES], etaphi_t seed_eta, etaphi_t seed_phi,
             PartialParticle partialParts[NPARTICLES]){

	#pragma HLS pipeline
    #pragma HLS array_partition variable=work complete
    #pragma HLS array_partition variable=incone complete
    #pragma HLS array_partition variable=partialParts complete

    etaphi_t seed_eta_broadcast[NPARTICLES];
    etaphi_t seed_phi_broadcast[NPARTICLES];
    #pragma HLS array_partition variable=seed_eta_broadcast complete
    #pragma HLS array_partition variable=seed_phi_broadcast complete
    //broadcast128<etaphi_t>(seed_eta, seed_eta_broadcast);
    //broadcast128<etaphi_t>(seed_phi, seed_phi_broadcast);

    // this block finds the jet constituents, and marks them
    JetsLoopParticleLoop:
    for (unsigned int i = 0; i < NPARTICLES; i++) {
        #pragma HLS unroll
        detaphi_t deta = dreg<detaphi_t>(work[i].hwEta - seed_eta);
        detaphi_t dphi = dreg<detaphi_t>(work[i].hwPhi - seed_phi);
        //detaphi_t deta = work[i].hwEta - dreg(seed_eta);
        //detaphi_t dphi = work[i].hwPhi - dreg(seed_phi);

        // phi wrap
        detaphi_t dphi0 = dphi > PI ? (detaphi_t) (TWOPI - dphi) : (detaphi_t) dphi;
        detaphi_t dphi1 = dphi < -PI ? (detaphi_t) (TWOPI + dphi) : (detaphi_t) dphi;
        detaphi_t dphiw = dphi > 0 ? dphi0 : dphi1;
        bool ic = deta*deta + dphiw*dphiw < R2CONE;
        incone[i] = ic;
        pt_t maybePt =  ic ? work[i].hwPt : pt_t(0);
        partialParts[i].hwPt = maybePt;
        partialParts[i].hwEta = (detaphi_t) deta;
        partialParts[i].hwPhi = (detaphi_t) dphiw;
/*#ifndef __SYNTHESIS__
        if(incone[i]){
            std::cout << " part: pt " << maybePt << ", eta: " << work[i].hwEta << ", phi: " << work[i].hwPhi << std::endl;
            std::cout << "    d: deta: " << deta << ", dphi: " << dphi << std::endl;
        }
#endif*/
    }
}

void formJet(const PartialParticle partialParts[NPARTICLES], pt_t & sum_pt, pt_etaphi_t & sum_pt_eta, pt_etaphi_t & sum_pt_phi, count_t & count){
    #pragma HLS pipeline
    pt_t sum_pts[NPARTICLES];
    pt_etaphi_t sum_pt_etas[NPARTICLES];
    pt_etaphi_t sum_pt_phis[NPARTICLES];
    #pragma HLS array_partition variable=sum_pts complete
    #pragma HLS array_partition variable=sum_pt_etas complete
    #pragma HLS array_partition variable=sum_pt_phis complete

    for(int i = 0; i < NPARTICLES; i++){
        #pragma HLS unroll
        sum_pts[i]     = partialParts[i].hwPt;
        sum_pt_etas[i] = partialParts[i].hwPt * partialParts[i].hwEta; // range is bounded since |deta| < 0.04 = 40 units
        sum_pt_phis[i] = partialParts[i].hwPt * partialParts[i].hwPhi;
    }
#ifndef __SYNTHESIS__
    for(int i = 0; i < NPARTICLES; i++){
        if(sum_pts[i] > 0){
            std::cout << "pt_deta: " << sum_pt_etas[i] << std::endl;
        }
    }
    for(int i = 0; i < NPARTICLES; i++){
        if(sum_pts[i] > 0){
            std::cout << "pt_dhi: " << sum_pt_phis[i] << std::endl;
        }
    }
#endif
    // sum the in-cone pt and pt-weighted delta-eta, delta-phis
    Op_add<pt_t> op_add_pt;
    sum_pt = reduce<pt_t, NPARTICLES, Op_add<pt_t>>(sum_pts, op_add_pt);
    Op_add<pt_etaphi_t> op_add_etaphi;
    sum_pt_eta = reduce<pt_etaphi_t, NPARTICLES, Op_add<pt_etaphi_t>>(sum_pt_etas, op_add_etaphi);
    sum_pt_phi = reduce<pt_etaphi_t, NPARTICLES, Op_add<pt_etaphi_t>>(sum_pt_phis, op_add_etaphi);

    // Compute the multiplicity
    count_t counts[NPARTICLES];
    #pragma HLS array_partition variable=counts complete
    for(int i = 0; i < NPARTICLES; i++){
        #pragma HLS unroll
        counts[i] = (count_t) (sum_pts[i] > 0);
    }
    Op_add<count_t> op_add_count;
    count = reduce<count_t, NPARTICLES, Op_add<count_t>>(counts, op_add_count);
}

void findSeed(const Particle work[NPARTICLES], Particle & seedp, etaphi_t & seed_eta, etaphi_t & seed_phi){
	#pragma HLS pipeline
    // Pick the highest pT particle as the seed
    Op_max<Particle> op_max;
    seedp = reduce<Particle, NPARTICLES, Op_max<Particle>>(work, op_max);
    seed_eta = seedp.hwPt > 0 ? seedp.hwEta : etaphi_t (0);
    seed_phi = seedp.hwPt > 0 ? seedp.hwPhi : etaphi_t (0);
#ifndef __SYNTHESIS__
    std::cout << " Seed: pt " << seedp.hwPt << ", eta: " << seed_eta << ", phi: " << seed_phi << std::endl;
#endif
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


void jet_loop(const Particle particles_in[NPARTICLES], Particle particles_out[NPARTICLES],
              PartialParticle partialParts[NPARTICLES], etaphi_t & seed_eta, etaphi_t & seed_phi) {
    #pragma HLS array_partition variable=particles_in complete
    #pragma HLS array_partition variable=particles_out complete
    #pragma HLS array_partition variable=partialParts complete
    #pragma HLS data_pack variable=particles_in
    #pragma HLS data_pack variable=particles_out
    #pragma HLS data_pack variable=partialParts
    #pragma HLS interface ap_none port=particles_out
    #pragma HLS interface ap_none port=partialParts
    #pragma HLS pipeline

    Particle seedp;
    bool incone[NPARTICLES];

    #pragma HLS array_partition variable=incone complete
    PartialParticle partialParts_int[NPARTICLES];
    etaphi_t seed_eta_int, seed_phi_int;
    for(int i = 0; i < NPARTICLES; i++){
        #pragma HLS unroll
        incone[i] = false;
        clear(partialParts_int[i]);
    }

    findSeed(particles_in, seedp, seed_eta_int, seed_phi_int);
    findJet(particles_in, incone, seed_eta_int, seed_phi_int, partialParts_int);
    copyInput(particles_in, particles_out);
    updateWork(particles_out, incone);

    for(int i = 0; i < NPARTICLES; i++){
        #pragma HLS unroll
        partialParts[i] = partialParts_int[i];
    }
    seed_eta = seed_eta_int;
    seed_phi = seed_phi_int;
}

void jet_compute(const PartialParticle particles[NPARTICLES], 
                 const etaphi_t seed_eta, const etaphi_t seed_phi, Jet & jet){
    #pragma HLS array_partition variable=particles complete
    #pragma HLS data_pack variable=particles
    #pragma HLS data_pack variable=jet
    #pragma HLS pipeline

    pt_t sum_pts[NPARTICLES];
    pt_etaphi_t sum_pt_etas[NPARTICLES];
    pt_etaphi_t sum_pt_phis[NPARTICLES];
    #pragma HLS array_partition variable=sum_pts complete
    #pragma HLS array_partition variable=sum_pt_etas complete
    #pragma HLS array_partition variable=sum_pt_phis complete

    pt_t sum_pt;
    pt_etaphi_t sum_pt_eta, sum_pt_phi;
    count_t count;
    formJet(particles, sum_pt, sum_pt_eta, sum_pt_phi, count);
    updateAxis(seed_eta, seed_phi, sum_pt, sum_pt_eta, sum_pt_phi, count, jet);
}

#ifndef __SYNTHESIS__
void algo_main(const Particle particles[NPARTICLES], Jet jets[NJETS]){
    // Function calling both Loop and Compute functions, intended for simulation
    Particle parts_in[NPARTICLES];
    Particle parts_out[NPARTICLES];
    PartialParticle partialParts[NPARTICLES];

    std::vector<Jet> jets_int;
    jets_int.resize(NJETS);
    for(int j = 0; j < NPARTICLES; j++){ parts_in[j] = particles[j]; }

    for(int i = 0; i < NJETS; i++){
        etaphi_t seed_eta, seed_phi;
        jet_loop(parts_in, parts_out, partialParts, seed_eta, seed_phi);
        Jet jet;
        jet_compute(partialParts, seed_eta, seed_phi, jet);
        jets_int[i] = jet;
        for(int j = 0; j < NPARTICLES; j++){ parts_in[j] = parts_out[j]; }
    }
    std::sort(jets_int.begin(), jets_int.end(), [](Jet i, Jet j) {return (i.hwPt > j.hwPt);});
    for(int i = 0; i < NJETS; i++){ jets[i] = jets_int.at(i); }
}
#endif
