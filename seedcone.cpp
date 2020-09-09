#include "algo.h"
#include "TreeReduce.h"
#include <cmath>
#include <cassert>
#include "hls_stream.h"
#include "data.h"

#ifndef __SYNTHESIS__
#include <cstdio>
#endif


void _invert_lut_init(ap_uint<18> table[1024]) {
    for (int i = 0; i < 1024; ++i) {
        int val = i > 16 ? (1 << 22)/i : 0; 
        assert(val >= 0 && val < (1 << 18));
        table[i] = val;
    }
}

void updateAxis(etaphi_t seed_eta, etaphi_t seed_phi,
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

void formJet(Particle work[NPARTICLES], bool incone[NPARTICLES], etaphi_t seed_eta, etaphi_t seed_phi, ap_uint<15> & sum_pt, ap_int<22> & sum_pt_eta, ap_int<22> & sum_pt_phi, ap_uint<5> & count){
	#pragma HLS pipeline
    // this block builds the jet out of the seed, and marks the used candidates
    // Reset seed pT to 0, as it will be clustered from work
	sum_pt = 0; sum_pt_eta = 0; sum_pt_phi = 0;
	count = (sum_pt > 0) ? 1 : 0;
    JetsLoopParticleLoop:
    for (unsigned int i = 0; i < NPARTICLES; ++i) {
        #pragma HLS unroll
        int deta = work[i].hwEta - seed_eta;
        int dphi = work[i].hwPhi - seed_phi;
        incone[i] = deta*deta + dphi*dphi < R2CONE;
        ap_uint<15> maybePt =  incone[i] ? ap_uint<15>(work[i].hwPt(14,0)) : ap_uint<15>(0);
        sum_pt     += maybePt;
        sum_pt_eta += maybePt * ap_int<7>(deta); // range is bounded since |deta| < 0.04 = 40 units
        sum_pt_phi += maybePt * ap_int<7>(dphi);
        count += (maybePt > 0);
    }
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
		ap_uint<15> sum_pt;
		ap_int<22> sum_pt_eta, sum_pt_phi;
		ap_uint<5> count;

		findSeed(work, seedp, seed_eta, seed_phi);
		formJet(work, incone, seed_eta, seed_phi, sum_pt, sum_pt_eta, sum_pt_phi, count);
		updateAxis(seed_eta, seed_phi, sum_pt, sum_pt_eta, sum_pt_phi, count, currentJet);
        updateWork(work, incone);
        addJetPtSorted(currentJet, myjet);
    }

    copyOutput(myjet, jet);
}

void jetLoopBody(Particle work[NPARTICLES], Jet myjet[NJETS]){
	#pragma HLS array_partition variable=myjet complete

	#pragma HLS pipeline II=LOOPII
	Jet currentJet;
	Particle seedp;
	etaphi_t seed_eta, seed_phi;
	ap_uint<15> sum_pt;
	ap_int<22> sum_pt_eta, sum_pt_phi;
	ap_uint<5> count;

	bool incone [NPARTICLES];
	#pragma HLS array_partition variable=incone complete

	findSeed(work, seedp, seed_eta, seed_phi);
	formJet(work, incone, seed_eta, seed_phi, sum_pt, sum_pt_eta, sum_pt_phi, count);
	updateAxis(seed_eta, seed_phi, sum_pt, sum_pt_eta, sum_pt_phi, count, currentJet);
	updateWork(work, incone);
	addJetPtSorted(currentJet, myjet);
}

void iterate0(Particle inParticles[NPARTICLES], Jet outJets[NJETS]){
	Particle workParticles[LOOPII][NPARTICLES];
	Jet workJets[NJETS];

	#pragma HLS array_partition variable=inParticles complete
	#pragma HLS array_partition variable=outJets complete
	#pragma HLS array_partition variable=workParticles complete
	#pragma HLS array_partition variable=workJets complete

	JetsLoop:
	for(int i = 0; i < NJETS; i++){
		#pragma HLS pipeline
		IILOOP:
		for(int j = 0; j < LOOPII; j++){
			Particle loopWorkParticles[NPARTICLES];
			#pragma HLS array_partition variable=loopWorkParticles complete
			for(int k = 0; k < NPARTICLES; k++){
				#pragma HLS unroll
				loopWorkParticles[k] = (i == 0) ? inParticles[k] : workParticles[j][k];
			}

			jetLoopBody(loopWorkParticles, workJets);

			if(i < NJETS-1){
				for(int k = 0; k < NPARTICLES; k++){
					#pragma HLS unroll
					workParticles[j][k] = loopWorkParticles[k];
				}
			}else{
				for(int k = 0; k < NJETS; k++){
					#pragma HLS unroll
					outJets[k] = workJets[k];
				}
			}
		}
	}

}

void iterate(hls::stream<bool> &reset, hls::stream<bool> &newInput, hls::stream<Particle> inParticles[NPARTICLES], hls::stream<Jet> outJets[NJETS]){
	#pragma HLS pipeline II=1
	#pragma HLS array_partition variable=inParticles complete
	#pragma HLS array_partition variable=outJets complete
	//static const int LOOPII=16;

	static Particle workParticles[LOOPII][NPARTICLES];
	static Particle tmp[NPARTICLES];
	static Particle store[NPARTICLES];
	static Jet workJets[LOOPII][NJETS];
	static Jet tmpJets[NJETS];
	#pragma HLS array_partition variable=workParticles complete
	#pragma HLS array_partition variable=tmp complete
	#pragma HLS array_partition variable=store complete
	#pragma HLS array_partition variable=workJets complete

	// Control signals
	static int nJetsDone[LOOPII];
	static bool working[LOOPII];
	static bool newIn[LOOPII];
	static int phase = 0;
	static bool newInputData = false;
	static bool resetData = false;

	MainLoop:
	//while(!resetData){
	for(int nEv = 0; nEv < 2; nEv++){
		#pragma HLS pipeline II=1 rewind

		// Input section
		// Read the input streams
		newInput.read_nb(newInputData);
		for(int i = 0; i < NPARTICLES; i++){
			inParticles[i].read_nb(tmp[i]);
		}

		// If the input is valid, copy it to the internal registers to wait for a slot in the algo loop pipeline
		if(newInputData){
			for(int i = 0; i < NPARTICLES; i++){
				store[i] = tmp[i];
			}
		}

		// Only copy from the store array to work when that slot is not in use
		// Rate of data in should be low enough that tmp is never clobbered
		if(!working[phase]){
			for(int i = 0; i < NPARTICLES; i++){
				workParticles[phase][i].hwPt = tmp[i].hwPt;
				workParticles[phase][i].hwEta = tmp[i].hwEta;
				workParticles[phase][i].hwPhi = tmp[i].hwPhi;

			}
		}

		// Signal that there's new input
		if(newInputData && !working[phase]){
			newIn[phase] = true;
		}else{
			newIn[phase] = false;
		}

		// Algorithm section
		jetLoopBody(workParticles[phase], workJets[phase]);

		// Control section
		reset.read_nb(resetData);
		bool workingN = working[phase];
		working[phase] = newIn[phase] || (workingN && nJetsDone[phase] < NJETS);
		int nJetsDoneN = nJetsDone[phase];
		nJetsDone[phase] = workingN ? nJetsDoneN + 1 : 0;

		// Output section
		if(nJetsDoneN == NJETS){
			for(int i = 0; i < NJETS; i++){
				#pragma HLS unroll
				tmpJets[i] = workJets[phase][i];
				outJets[i].write_nb(tmpJets[i]);
			}
		}

		phase = phase >= LOOPII ? 0 : phase + 1;
	}


}

void multipleEvents(hls::stream<bool> newInput, hls::stream<Particle> inParticles[NPARTICLES], hls::stream<Jet> outJets[NJETS]){
	#pragma HLS pipeline II=1
	static const int NEVENTS_INFLIGHT = 6;
	Particle workParticles[NEVENTS_INFLIGHT][NPARTICLES];
	Particle tmpParticles[NPARTICLES];
	Jet tmpJets[NEVENTS_INFLIGHT][NJETS];
	#pragma HLS array_partition variable=inParticles complete
	#pragma HLS array_partition variable=workParticles dim=1 complete
	#pragma HLS array_partition variable=tmpParticles dim=0 complete
	#pragma HLS array_partition variable=tmpJets dim=1 complete
	#pragma HLS array_partition variable=outJets complete

	bool newIn;
	int iEvt = 0;

	while(true){
		newInput.read_nb(newIn);
		if(newIn){
			for(int np = 0; np < NPARTICLES; np++){
				#pragma HLS unroll
				inParticles[np].read_nb(tmpParticles[np]);
			}
			iEvt++;
		}

		ProcessLoop:
		for(int NEVT = 0; NEVT < NEVENTS_INFLIGHT; NEVT++){
			#pragma HLS pipeline II=1
			algo_main(workParticles[NEVT], tmpJets[NEVT]);
		}
	}
}


