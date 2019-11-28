#include "src/data.h"
#include "src/algo.h"
#include <cmath>

inline int dR2(const Particle & p1, const Particle &p2) {
    int deta = (p1.hwEta-p2.hwEta);
    int dphi = (p1.hwPhi-p2.hwPhi);
    return deta*deta + dphi*dphi;
}

void algo_main_ref(const Particle particles[NPARTICLES], Jet jet[NJETS]) {
    bool used[NPARTICLES];
    for (unsigned int i = 0; i < NPARTICLES; ++i) {
        used[i] = false;
    }

    Jet myjet[NJETS+MOREJETS];
    for (unsigned int j = 0; j < NJETS+MOREJETS; ++j) {
        clear(myjet[j]);
    }

    for (unsigned int j = 0; j < NJETS+MOREJETS; ++j) {
        // find a seed
        int iseed = -1;
        for (unsigned int i = 0; i < NPARTICLES; ++i) {
            if (used[i] || particles[i].hwPt == 0) continue;
            if (iseed == -1 || particles[i].hwPt > particles[iseed].hwPt) {
                iseed = i;
            }
        }
        // if not found, stop
        if (iseed == -1) break;
        int sum_pt = 0, sum_pt_eta = 0, sum_pt_phi = 0, count = 0;
        for (unsigned int i = 0; i < NPARTICLES; ++i) {
            if (dR2(particles[i], particles[iseed]) < R2CONE) {
                sum_pt     += particles[i].hwPt;
                sum_pt_eta += particles[i].hwPt * (particles[i].hwEta - particles[iseed].hwEta);
                sum_pt_phi += particles[i].hwPt * (particles[i].hwPhi - particles[iseed].hwPhi);
                used[i] = true;
                count += (particles[i].hwPt > 0);
            }

        }
        if (sum_pt < JET_PT_CUT) continue;
        int jet_eta = particles[iseed].hwEta + (sum_pt_eta / sum_pt);
        int jet_phi = particles[iseed].hwPhi + (sum_pt_phi / sum_pt);
        if (std::abs(jet_eta) > FIDUCIAL_ETA_PHI || std::abs(jet_phi) > FIDUCIAL_ETA_PHI) continue;
        myjet[j].hwPt  = sum_pt;
        myjet[j].hwEta = jet_eta;
        myjet[j].hwPhi = jet_phi;
        myjet[j].nCand = count;
        myjet[j].iSeed = iseed;
    }

    // sort the jets
    for (unsigned int j = 0; j < NJETS+MOREJETS; ++j) {
        for (unsigned int j2 = j+1; j2 < NJETS+MOREJETS; ++j2) {
            if (myjet[j2].hwPt > myjet[j].hwPt) {
                std::swap(myjet[j], myjet[j2]);
            }
        }
    }
    for (unsigned int j = 0; j < NJETS; ++j) {
        jet[j] = myjet[j];
    }
}

