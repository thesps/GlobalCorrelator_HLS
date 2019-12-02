#include "algo.h"
#include <cmath>

#ifndef __SYNTHESIS__
#include <cstdio>
#endif

inline bool inCone(etaphi_t eta1, etaphi_t phi1, etaphi_t eta2, etaphi_t phi2) {
    int deta = (eta1-eta2);
    int dphi = (phi1-phi2);
    return deta*deta + dphi*dphi < R2CONE;
}

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
#ifndef __SYNTHESIS__
    //printf("Found seed with pt %+7.2f i = %d\n", pt_seed*0.25, iseed);
#endif
    return seed_pt > 0;
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
    #pragma HLS pipeline II=2
    bool used[NPARTICLES];
    #pragma HLS array_partition variable=used complete
    for (unsigned int i = 0; i < NPARTICLES; ++i) {
        used[i] = false;
    }

    Jet myjet[NJETS+MOREJETS];
    #pragma HLS array_partition variable=myjet complete

    for (unsigned int j = 0; j < NJETS+MOREJETS; ++j) {
#ifndef __SYNTHESIS__
        //printf("Jet finding, iteration %d\n", j);
#endif
        bool found = false;
        ap_uint<5> iseed; etaphi_t seed_eta, seed_phi;
        found = find_seed(particles, used, iseed, seed_eta, seed_phi);
        int sum_pt = 0, sum_pt_eta = 0, sum_pt_phi = 0, count = 0;
        for (unsigned int i = 0; i < NPARTICLES; ++i) {
            bool addme = found && inCone(particles[i].hwEta, particles[i].hwPhi, seed_eta, seed_phi);
            pt_t maybePt =  addme ? particles[i].hwPt : pt_t(0);
            sum_pt     += maybePt;
            sum_pt_eta += maybePt * (particles[i].hwEta - seed_eta);
            sum_pt_phi += maybePt * (particles[i].hwPhi - seed_phi);
            used[i] = used[i] || addme;
            count += addme && (particles[i].hwPt > 0);
        }
        etaphi_t jet_eta = found ? int(seed_eta) + (sum_pt_eta / sum_pt) : 0;
        etaphi_t jet_phi = found ? int(seed_phi) + (sum_pt_phi / sum_pt) : 0;
        bool fiducial = (sum_pt >= JET_PT_CUT) &&
                        (-FIDUCIAL_ETA_PHI <= jet_eta && jet_eta <= FIDUCIAL_ETA_PHI) &&
                        (-FIDUCIAL_ETA_PHI <= jet_phi && jet_phi <= FIDUCIAL_ETA_PHI);
        if (fiducial) {
            myjet[j].hwPt  = sum_pt;
            myjet[j].hwEta = jet_eta;
            myjet[j].hwPhi = jet_phi;
            myjet[j].nCand = count;
            myjet[j].iSeed = iseed;
        } else {
            clear(myjet[j]);
        }
    }
    sort_and_crop<Jet,NJETS+MOREJETS,NJETS>(myjet, jet);

}


