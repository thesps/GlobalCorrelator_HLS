template<int N>
void streamSortParticles(const Particle particlesIn[N], Particle particlesOut[N]){
    // Add one jet to the list, keeping them pt sorted
    #pragma HLS latency min=N-1 max=N
    #pragma HLS array_partition variable=particlesOut complete
    // tmp should get shifted along the list of jets
	for(int i = 0; i < N; i++){
		Particle p = particlesIn[i];
	    Particle forward = p;
	    JetStreamSort:
	    for(int i = 0; i < NJETS; i++){
			#pragma HLS pipeline
	        Particle tmp;
	        bool GT = forward.hwPt >= particlesOut[i].hwPt;
	        tmp = GT ? particlesOut[i] : forward;
	        particlesOut[i] = GT ? forward : particlesOut[i];
	        forward = tmp;
	    }
	}

}
