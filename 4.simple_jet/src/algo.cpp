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

bool find_seed(const Particle particles[NPARTICLES], const bool used[NPARTICLES], int & iseed, etaphi_t & seed_eta, etaphi_t & seed_phi) {
    bool ret = false;
    iseed = 0;
    seed_eta = 0;
    seed_phi = 0;
    pt_t pt_seed = 0; 
    for (unsigned int i = 0; i < NPARTICLES; ++i) {
        if (!used[i] && particles[i].hwPt > pt_seed) {
            iseed = i;
            pt_seed = particles[i].hwPt;
            seed_eta = particles[i].hwEta;
            seed_phi = particles[i].hwPhi;
            ret = true;
        }
    }
#ifndef __SYNTHESIS__
    /*if (ret) { 
        printf("Found seed with pt %+7.2f i = %d\n", pt_seed*0.25, iseed);
    } else {
        printf("No seed found\n");
    }*/
#endif
    return ret;
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
    #pragma HLS pipeline II=9
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
        int iseed; etaphi_t seed_eta, seed_phi;
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


