#include <cstdio>
#include <cstdlib>
#include <iostream>
#include "algo.h"

#define NTEST 6

int main() {
    Particle particles[NPARTICLES];
    Jet jet[NJETS], ref[NJETS], ak4[NJETS];
    int test = 0;
    FILE *dump = fopen("pfcands_ttbar.txt","r"); 
    if (!dump) return 3;
    printf("Opened dump file\n");
    while (!feof(dump) && (++test <= NTEST)) {

        for (unsigned int i = 0; i < NPARTICLES; ++i) clear(particles[i]);
        for (unsigned int j = 0; j < NJETS; ++j) clear(ak4[j]);

        int ncands, nak4, id, ntrunk = 0; float pt, eta, phi;
        if (fscanf(dump, "Event with %d candidates, %d jets in the selected region\n", &ncands, &nak4) != 2) return 2; 
        for (int i = 0; i < ncands; ++i) {
            if (fscanf(dump, "   cand pt %f eta %f phi %f  id %d", &pt, &eta, &phi, &id) != 4) return 2;
            if (i < NPARTICLES) {
                particles[i].hwPt =  pt; //int(pt/0.25 + 0.49999);
                particles[i].hwEta = (eta * 100 / 128); //int(eta/0.01);
                particles[i].hwPhi = (phi * 100 / 128); //int(phi/0.01);
            } else { ntrunk++; }
        }
        for (int i = 0; i < nak4; ++i) {
            if (fscanf(dump, "   jet  pt %f eta %f phi %f  constituents %d", &pt, &eta, &phi, &id) != 4) return 2;
            if (i < NJETS) {
                ak4[i].hwPt =  pt; //int(pt/0.25 + 0.49999);
                ak4[i].hwEta = (eta * 100 / 128); //int(eta/0.01);
                ak4[i].hwPhi = (phi * 100 / 128); //int(phi/0.01);
                ak4[i].iSeed = 0;
                ak4[i].nCand = id;
            }
        }
        fscanf(dump, "\n");

        // run the algorithm
        algo_main(particles, jet);
        algo_main_ref(particles, ref);
        // check the output
        printf("Event with %d particles (%d truncated away), %d ak4 jets\n", ncands, ntrunk, nak4);
        bool ok = true;
        for (int i = 0; i < NJETS; ++i) {
            std::cout << "Jet    HW:" << i << ": pt " << jet[i].hwPt << " eta " << 
                jet[i].hwEta << " phi " << jet[i].hwPhi << " iseed " << jet[i].iSeed << " ncand " <<
                jet[i].nCand << std::endl;
            /*std::cout << "      ref:" << i << ": pt " << ref[i].hwPt << " eta " << 
                ref[i].hwEta << " phi " << ref[i].hwPhi << " iseed " << ref[i].iSeed << " ncand " <<
                ref[i].nCand << std::endl;*/
            std::cout << "      ak4:" << i << ": pt " << ak4[i].hwPt << " eta " << 
                ak4[i].hwEta << " phi " << ak4[i].hwPhi << " iseed " << ak4[i].iSeed << " ncand " <<
                ak4[i].nCand << std::endl;
            if (jet[i].hwPt != ref[i].hwPt || jet[i].hwEta != ref[i].hwEta || jet[i].hwPhi != ref[i].hwPhi || jet[i].nCand != ref[i].nCand) {
                ok = false;
            } 
        }
        if (!ok) { printf("MISMATCH\n"); return 1; }
        printf("\n");
    }
    printf("Passed all %d tests\n", NTEST);
    return 0;
}
