#include "algo.h"
#include <cmath>


pt_t algo_main(Particle particles[NPARTICLES]) {
    #pragma HLS ARRAY_PARTITION variable=particles complete
    #pragma HLS pipeline II=1

    int sum = 0;
    for (unsigned int i = 0; i < NPARTICLES; ++i) {
        #pragma HLS unroll
        //====  Naive version, with an if guarding the sum
        if (-240 <= particles[i].hwEta && particles[i].hwEta <= 240) {
            sum += particles[i].hwPt;
        }
        //====  Version where the sum is always computed, but we add zero in some entries
        //bool central = (-240 <= particles[i].hwEta && particles[i].hwEta <= 240);
        //sum += (central ? particles[i].hwPt : pt_t(0));
    }
    return pt_t(sum);
}


