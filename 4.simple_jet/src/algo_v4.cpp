#include "algo.h"
#include <cmath>
#include <cassert>

#ifndef __SYNTHESIS__
#include <cstdio>
#endif


template<unsigned int N>
void find_seed_rec(const Particle particles[N], const bool used[N], ap_uint<5> i0, ap_uint<5> & iseed, etaphi_t & seed_eta, etaphi_t & seed_phi, pt_t & seed_pt)
{
    ap_uint<5> iseed_up, iseed_down; 
    etaphi_t seed_eta_up, seed_phi_up, seed_eta_down, seed_phi_down;
    pt_t seed_pt_up, seed_pt_down;
    find_seed_rec< N/2 >(&particles[ 0 ], &used[ 0 ], i0+0,   iseed_up,   seed_eta_up,   seed_phi_up,   seed_pt_up);
    find_seed_rec<N-N/2>(&particles[N/2], &used[N/2], i0+N/2, iseed_down, seed_eta_down, seed_phi_down, seed_pt_down);
    if (seed_pt_up >= seed_pt_down) {
        iseed = iseed_up;
        seed_pt = seed_pt_up;
        seed_eta = seed_eta_up;
        seed_phi = seed_phi_up;
    } else {
        iseed = iseed_down;
        seed_pt = seed_pt_down;
        seed_eta = seed_eta_down;
        seed_phi = seed_phi_down;
    }

}
template<>
void find_seed_rec<1>(const Particle particles[1], const bool used[1], ap_uint<5> i0, ap_uint<5> & iseed, etaphi_t & seed_eta, etaphi_t & seed_phi, pt_t & seed_pt)
{
    seed_pt  = !used[0] ? particles[0].hwPt  : pt_t(0);
    seed_eta = !used[0] ? particles[0].hwEta : etaphi_t(0);
    seed_phi = !used[0] ? particles[0].hwPhi : etaphi_t(0);
    iseed    = !used[0] ? i0                 : ap_uint<5>(0);
}

bool find_seed(const Particle particles[NPARTICLES], const bool used[NPARTICLES], ap_uint<5> & iseed, etaphi_t & seed_eta, etaphi_t & seed_phi) {
    pt_t seed_pt;
    find_seed_rec<NPARTICLES>(particles, used, 0, iseed, seed_eta, seed_phi, seed_pt);
    return seed_pt > 0;
}

void _invert_lut_init(ap_uint<18> table[1024]) {
    for (int i = 0; i < 1024; ++i) {
        int val = i > 16 ? (1 << 22)/i : 0; 
        assert(val >= 0 && val < (1 << 18));
        table[i] = val;
    }
}
void reconstruct_jet(bool valid_seed, ap_uint<5> iseed, etaphi_t seed_eta, etaphi_t seed_phi, 
                     ap_uint<15> sum_pt, ap_int<22> sum_pt_eta, ap_int<22> sum_pt_phi, ap_uint<5> count, Jet & jet) {
    ap_uint<18> _inv_table[1024]; _invert_lut_init(_inv_table);

    ap_uint<10> den; ap_int<17> num_eta, num_phi;
    if (sum_pt[14] || sum_pt[13] || sum_pt[12]) { // sum_pt > 4*1024
        den     = sum_pt >> 5;
        num_eta = sum_pt_eta >> 5;
        num_phi = sum_pt_phi >> 5;
    } else if (sum_pt[11] || sum_pt[10]) {
        den     = sum_pt >> 2;
        num_eta = sum_pt_eta >> 2;
        num_phi = sum_pt_phi >> 2;
    } else {
        den     = sum_pt;
        num_eta = sum_pt_eta;
        num_phi = sum_pt_phi;
    }
    assert((sum_pt < JET_PT_CUT) || den > 16);
    ap_uint<18> inv_den = _inv_table[den];
    etaphi_t jet_eta = seed_eta + etaphi_t((num_eta * inv_den) >> 22);
    etaphi_t jet_phi = seed_phi + etaphi_t((num_phi * inv_den) >> 22);
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
        jet.iSeed = iseed;
    } else {
        clear(jet);
    }

}

template<typename T, int NIn, int NOut>
void sort_and_crop(const T in[NIn], T out[NOut]) {
    T tmp[NOut];

    for (int iout = 0; iout < NOut; ++iout) {
        clear(tmp[iout]);
    }

    for (int it = 0; it < NIn; ++it) {
        for (int iout = NOut-1; iout >= 0; --iout) {
            if (tmp[iout].hwPt <= in[it].hwPt) {
                if (iout == 0 || tmp[iout-1].hwPt > in[it].hwPt) {
                    tmp[iout] = in[it];
                } else {
                    tmp[iout] = tmp[iout-1];
                }
            }
        }
    }

    for (int iout = 0; iout < NOut; ++iout) {
        out[iout] = tmp[iout];
    }

}


void algo_main(const Particle particles[NPARTICLES], Jet jet[NJETS]) {
    #pragma HLS array_partition variable=particles complete
    #pragma HLS array_partition variable=jet complete
    #pragma HLS interface ap_none port=jet 
    #pragma HLS pipeline II=1
    bool used[NPARTICLES];
    #pragma HLS array_partition variable=used complete
    for (unsigned int i = 0; i < NPARTICLES; ++i) {
        used[i] = false;
    }

    Jet myjet[NJETS+MOREJETS];
    #pragma HLS array_partition variable=myjet complete

    for (unsigned int j = 0; j < NJETS+MOREJETS; ++j) {
        bool found = false;
        ap_uint<5> iseed; etaphi_t seed_eta, seed_phi;
        found = find_seed(particles, used, iseed, seed_eta, seed_phi);
#ifndef __SYNTHESIS__
        //printf("HW algo iter %d: seed %2d of pt %.2f, valid? %1d\n", j, int(iseed), particles[iseed].hwPt*0.25, int(found));
#endif
        ap_uint<15> sum_pt = 0; 
        ap_int<22> sum_pt_eta = 0, sum_pt_phi = 0;
        ap_uint<5> count = 0;
        for (unsigned int i = 0; i < NPARTICLES; ++i) {
            int deta = particles[i].hwEta - seed_eta;
            int dphi = particles[i].hwPhi - seed_phi;
            bool incone = deta*deta + dphi*dphi < R2CONE;
            bool addme = found && incone && !used[i];
            ap_uint<15> maybePt =  addme ? ap_uint<15>(particles[i].hwPt(14,0)) : ap_uint<15>(0);
            sum_pt     += maybePt;
            sum_pt_eta += maybePt * ap_int<7>(deta); // range is bounded since |deta| < 0.04 = 40 units
            sum_pt_phi += maybePt * ap_int<7>(dphi);
            used[i] = used[i] || addme;
            count += (maybePt > 0);
#ifndef __SYNTHESIS__
            //if (addme) printf("                add cand %d, cluster pt now %.2f, ncand %d\n", i, sum_pt*0.25, count);
#endif
        }
#ifndef __SYNTHESIS__
        //printf("HW algo iter %d: cluster of pt %.2f, ncand %d, valid? %1d\n", j, sum_pt*0.25, count, int(found));
#endif
        reconstruct_jet(found, iseed, seed_eta, seed_phi, sum_pt, sum_pt_eta, sum_pt_phi, count, myjet[j]);
    }
    sort_and_crop<Jet,NJETS+MOREJETS,NJETS>(myjet, jet);

}


