#include "src/data.h"
#include "src/algo.h"
#include <cmath>

pt_t algo_main_ref(Particle particles[NPARTICLES]) {
    int sum = 0;
    for (unsigned int i = 0; i < NPARTICLES; ++i) {
        if (std::abs(particles[i].hwEta) <= 240) {
            sum += particles[i].hwPt;
        }
    }
    return pt_t(sum);
}

