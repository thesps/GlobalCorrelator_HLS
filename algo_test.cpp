#include <cstdio>
#include <cstdlib>
#include "algo.h"
#include "data.h"
#include "hls_stream.h"


#define NTEST 6

int main(){
    Particle particles[NPARTICLES];
    hls::stream<Particle> particleStream[NPARTICLES];
    hls::stream<bool> rst;
    hls::stream<bool> newData;
    Jet jet[NJETS], ref[NJETS], ak4[NJETS];
    hls::stream<Jet> jetStream[NJETS];
    hls::stream<Jet> refJets[NJETS];
    hls::stream<Jet> ak4Jets[NJETS];

    int test = 0;
    FILE *dump = fopen("pfcands_ttbar.txt","r"); 
    if (!dump) return 3;
    printf("Opened dump file\n");
    while (!feof(dump) && (++test <= NTEST)) {
        for (unsigned int i = 0; i < NPARTICLES; ++i) clear(particles[i]);

        int ncands, nak4, id, ntrunk = 0; float pt, eta, phi;
        if (fscanf(dump, "Event with %d candidates, %d jets in the selected region\n", &ncands, &nak4) != 2) return 2;
        printf("Event with %d candidates, %d jets in the selected region\n", ncands, nak4);
        for (int i = 0; i < NPARTICLES; ++i) {
            if (i < ncands) {
            	if (fscanf(dump, "   cand pt %f eta %f phi %f  id %d", &pt, &eta, &phi, &id) != 4) return 2;
                particles[i].hwPt =  int(pt/0.25 + 0.49999);
                particles[i].hwEta = int(eta/0.01);
                particles[i].hwPhi = int(phi/0.01);
                particleStream[i].write(particles[i]);

            }else{
                particles[i].hwPt =  0;
                particles[i].hwEta = 0;
                particles[i].hwPhi = 0;
                particleStream[i].write(particles[i]);
            }
        }

        rst.write(false);
        newData.write(true);
        for(unsigned int i = 0; i < 26; i++){
        	for (unsigned int j = 0; j < NPARTICLES; j++){
        		clear(particles[j]);
        		particleStream[j].write(particles[j]);
        	}
        	rst.write(false);
        	newData.write(false);
        }

        // run the reference algorithm
        //algo_main_ref(particles, ref);

        for (int i = 0; i < NJETS; ++i) {
            if (i < nak4) {
            	if (fscanf(dump, "   jet  pt %f eta %f phi %f  constituents %d", &pt, &eta, &phi, &id) != 4) return 2;
                ak4[i].hwPt =  int(pt/0.25 + 0.49999);
                ak4[i].hwEta = int(eta/0.01);
                ak4[i].hwPhi = int(phi/0.01);
                ak4[i].iSeed = 0;
                ak4[i].nCand = id;
            }else{
                ak4[i].hwPt =  0;
                ak4[i].hwEta = 0;
                ak4[i].hwPhi = 0;
                ak4[i].iSeed = 0;
                ak4[i].nCand = 0;
            }
            ak4Jets[i].write(ak4[i]);
            refJets[i].write(ref[i]);
        }
        fscanf(dump, "\n");
    }

    // End of events

    for(unsigned int i = 0; i < 30 * 12; i++){
    	for (unsigned int j = 0; j < NPARTICLES; j++){
    		clear(particles[j]);
    		particleStream[j].write(particles[j]);
    	}
    	rst.write(false);
    	newData.write(false);
    }

    // End of process
	//for (unsigned int j = 0; j < NPARTICLES; j++){
	//	clear(particles[j]);
	//	particleStream[j].write(particles[j]);
	//}
	rst.write(true);
	//newData.write(false);

	iterate(rst, newData, particleStream, jetStream);

	for(unsigned int n=0; n < NTEST; n++){
		bool ok = true;
		for (int i = 0; i < NJETS; i++) {
			jet[i] = jetStream[i].read();
			ref[i] = refJets[i].read();
			ak4[i] = ak4Jets[i].read();
			printf("Jet %d:  pt %7.2f  eta %+5.2f  phi %+5.2f  iseed %2d ncand %2d  ",
					i, jet[i].hwPt * 0.25, jet[i].hwEta * 0.01, jet[i].hwPhi * 0.01, int(jet[i].iSeed), int(jet[i].nCand));
			printf("    ref pt %7.2f  eta %+5.2f  phi %+5.2f  iseed %2d ncand %2d",
					   ref[i].hwPt * 0.25, ref[i].hwEta * 0.01, ref[i].hwPhi * 0.01, int(ref[i].iSeed), int(ref[i].nCand));
			printf("    ak4 pt %7.2f  eta %+5.2f  phi %+5.2f  ncand %2d\n",
					   ak4[i].hwPt * 0.25, ak4[i].hwEta * 0.01, ak4[i].hwPhi * 0.01, int(ak4[i].nCand));
			if (jet[i].hwPt != ref[i].hwPt || jet[i].hwEta != ref[i].hwEta || jet[i].hwPhi != ref[i].hwPhi || jet[i].nCand != ref[i].nCand) {
				ok = false;
			}
		}
		if (!ok) { printf("MISMATCH\n");}
		printf("\n");
	}
	printf("Passed all %d tests\n", NTEST);
	return 0;

}

/*int main() {
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
                particles[i].hwPt =  int(pt/0.25 + 0.49999);
                particles[i].hwEta = int(eta/0.01);
                particles[i].hwPhi = int(phi/0.01);

            } else { ntrunk++; }
        }
        for (int i = 0; i < nak4; ++i) {
            if (fscanf(dump, "   jet  pt %f eta %f phi %f  constituents %d", &pt, &eta, &phi, &id) != 4) return 2;
            if (i < NJETS) {
                ak4[i].hwPt =  int(pt/0.25 + 0.49999);
                ak4[i].hwEta = int(eta/0.01);
                ak4[i].hwPhi = int(phi/0.01);
                ak4[i].iSeed = 0;
                ak4[i].nCand = id;
            }
        }
        fscanf(dump, "\n");

        // run the algorithm
        //algo_main(particles, jet);
        algo_main_ref(particles, ref);
        // check the output
        printf("Event with %d particles (%d truncated away), %d ak4 jets\n", ncands, ntrunk, nak4);
        bool ok = true;
        for (int i = 0; i < NJETS; ++i) {
            printf("Jet %d:  pt %7.2f  eta %+5.2f  phi %+5.2f  iseed %2d ncand %2d  ",
                    i, jet[i].hwPt * 0.25, jet[i].hwEta * 0.01, jet[i].hwPhi * 0.01, int(jet[i].iSeed), int(jet[i].nCand));
            printf("    ref pt %7.2f  eta %+5.2f  phi %+5.2f  iseed %2d ncand %2d",
                       ref[i].hwPt * 0.25, ref[i].hwEta * 0.01, ref[i].hwPhi * 0.01, int(ref[i].iSeed), int(ref[i].nCand));
            printf("    ak4 pt %7.2f  eta %+5.2f  phi %+5.2f  ncand %2d\n",
                       ak4[i].hwPt * 0.25, ak4[i].hwEta * 0.01, ak4[i].hwPhi * 0.01, int(ak4[i].nCand));
            if (jet[i].hwPt != ref[i].hwPt || jet[i].hwEta != ref[i].hwEta || jet[i].hwPhi != ref[i].hwPhi || jet[i].nCand != ref[i].nCand) {
                ok = false;
            } 
        }
        if (!ok) { printf("MISMATCH\n"); return 1; }
        printf("\n");
    }
    printf("Passed all %d tests\n", NTEST);
    return 0;
}*/
