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
        printf("Seed:  pt %7.2f  eta %+5.2f  phi %+5.2f\n",
                            iseed, jet[iseed].hwPt * 0.25, jet[iseed].hwEta * 0.01, jet[iseed].hwPhi * 0.01, int(jet[iseed].iSeed), int(jet[iseed].nCand));
        // if not found, stop
        if (iseed == -1) break;
        //printf("Ref algo iter %d: found seed %2d of pt %.2f, eta %+.2f, phi %+.2f\n", j, iseed, particles[iseed].hwPt * 0.25, particles[iseed].hwEta*0.01, particles[iseed].hwPhi*0.01);
        int sum_pt = 0, sum_pt_eta = 0, sum_pt_phi = 0, count = 0;
        for (unsigned int i = 0; i < NPARTICLES; ++i) {
            if (dR2(particles[i], particles[iseed]) < R2CONE && !used[i]) {
                sum_pt     += particles[i].hwPt;
                sum_pt_eta += particles[i].hwPt * (particles[i].hwEta - particles[iseed].hwEta);
                sum_pt_phi += particles[i].hwPt * (particles[i].hwPhi - particles[iseed].hwPhi);
                used[i] = true;
                count += (particles[i].hwPt > 0);
                //printf("                add cand %d, pt %.2f, eta %+.2f, phi %+.2f, cluster pt now %.2f, ncand %d\n", i, particles[i].hwPt*0.25, particles[i].hwEta*0.01, particles[i].hwPhi*0.01, sum_pt*0.25, count);
            }

        }
        //printf("Ref algo iter %d: found cluster of total pt %.2f, ncand %d\n", j, sum_pt * 0.25, count);
        if (sum_pt < JET_PT_CUT) continue;
        // --- with a division
        //int jet_eta = particles[iseed].hwEta + (sum_pt_eta / sum_pt);
        //int jet_phi = particles[iseed].hwPhi + (sum_pt_phi / sum_pt);
        // --- with the approx division we use in the FPGA
        int sum_pt_den = sum_pt;
        if (sum_pt > 4096) {
            sum_pt_den = sum_pt >> 5; sum_pt_eta = sum_pt_eta >> 5; sum_pt_phi = sum_pt_phi >> 5; 
        } else if (sum_pt > 1024) {
            sum_pt_den = sum_pt >> 2; sum_pt_eta = sum_pt_eta >> 2; sum_pt_phi = sum_pt_phi >> 2; 
        }
        int jet_eta = particles[iseed].hwEta + etaphi_t((sum_pt_eta * ((1<<22)/sum_pt_den)) >> 22);
        int jet_phi = particles[iseed].hwPhi + etaphi_t((sum_pt_phi * ((1<<22)/sum_pt_den)) >> 22);
        if (std::abs(jet_eta) > FIDUCIAL_ETA_PHI || std::abs(jet_phi) > FIDUCIAL_ETA_PHI) continue;
        //printf("Ref algo iter %d: found jet pt %.2f, eta %+5.3f, phi %+5.3f\n", j, sum_pt * 0.25, jet_eta*0.01, jet_phi*0.01);
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

