#include <cstdio>
#include "src/algo.h"

#define NTEST 500

int main() {
    srand(42);
    Particle particles[NPARTICLES];
    for (int test = 1; test <= NTEST; ++test) {
        for (int i = 0; i < NPARTICLES; ++i) {
            particles[i].hwPt  = int((rand()/float(RAND_MAX) * 20 + 4)/0.25);
            particles[i].hwEta = int((rand()/float(RAND_MAX/2) - 1)*5 / 0.01);
        }
        // run the algorithm
        pt_t hw  = algo_main(particles);
        pt_t ref = algo_main_ref(particles);
        // check the output
        if (hw != ref) {
            printf("Error in test %d\n", test);
            for (int i = 0; i < NPARTICLES; ++i) {
                printf("   Particle %2d hwPt %5d  hwEta %+4d\n", i, 
                            int(particles[i].hwPt), int(particles[i].hwEta));
            }
            printf("  Algo result     :  %6d\n", int(hw));
            printf("  Reference result:  %6d\n", int(ref));
            printf("\n");
            return 1;
        }
    }
    printf("Passed all %d tests\n", NTEST);
    return 0;
}
