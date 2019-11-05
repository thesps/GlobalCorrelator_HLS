#include "algo.h"
#include <cmath>
#include <cassert>
#ifndef __SYNTHESIS__
#include <cstdio>
#endif

pt_t algo_main(Particle particles[NPARTICLES]) {
    #pragma HLS ARRAY_PARTITION variable=particles complete
    #pragma HLS pipeline II=1
    int sum1 = 0;
    for (unsigned int i = 0; i < NPARTICLES/2; ++i) {
        #pragma HLS unroll
        if (-240 <= particles[i].hwEta && particles[i].hwEta <= 240) {
            sum1 += particles[i].hwPt;
        }
    }
    int sum2 = 0;
    for (unsigned int i = NPARTICLES/2; i < NPARTICLES; ++i) {
        #pragma HLS unroll
        if (-240 <= particles[i].hwEta && particles[i].hwEta <= 240) {
            sum2 += particles[i].hwPt;
        }
    }
    return pt_t(sum1+sum2);
}

