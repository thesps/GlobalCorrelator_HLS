#include "src/data.h"
#include "src/algo.h"
#include <cmath>

int approx_sqrt(int i) {
    assert(i >= 0);
    int multiplier = 1;
    while (i >= 4*1024) {
        i /= 4;
        multiplier *= 2;
    }
    return int(std::sqrt(float(i)))*multiplier;
}

void algo_main_ref(const Particle particles[NPARTICLES], pxy_t & met_px, pxy_t & met_py, pt_t & met_pt) {
    int sumpx = 0, sumpy = 0;
    for (unsigned int i = 0; i < NPARTICLES; ++i) {
        float phi = int(particles[i].hwPhi)*0.01f;
        sincos_t sinphi = std::max(std::min<int>(std::sin(phi) * 256, +255), -255);
        sincos_t cosphi = std::max(std::min<int>(std::cos(phi) * 256, +255), -255);
        pxy_t px = (particles[i].hwPt * cosphi) >> 6; // throw away some bits already here
        pxy_t py = (particles[i].hwPt * sinphi) >> 6; // throw away some bits already here
        sumpx += px;
        sumpy += py;
    }
    met_px = sumpx;
    met_py = sumpy;
    // throw away some more bits before squaring
    sumpx = sumpx >> 4; // now 1 unit = 1 GeV
    sumpy = sumpy >> 4; // max (8.191 TeV) is 2^12-1
    
    int sumpt2 = (sumpx * sumpx + sumpy * sumpy); 
    // 1 unit = 1 GeV^2, max is 2^24-1

    int met = sumpt2 > 512*512-1 ? 512 : approx_sqrt(sumpt2);
    met_pt = met << 2;
}

void algo_main_ref_float(const Particle particles[NPARTICLES], pxy_t & met_px, pxy_t & met_py, pt_t & met_pt) {
    float sumpx = 0, sumpy = 0;
    for (unsigned int i = 0; i < NPARTICLES; ++i) {
        float phi = int(particles[i].hwPhi)*0.01f, pt = 0.25f*int(particles[i].hwPt);
        sumpx += std::cos(phi) * pt;
        sumpy += std::sin(phi) * pt;
    }
    met_px = pxy_t(int(sumpx/0.25f*4));
    met_py = pxy_t(int(sumpy/0.25f*4));

    float pt = std::hypot(sumpx, sumpy);
    met_pt = pt_t(int(pt/0.25f));
}

