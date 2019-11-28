#include <cstdio>
#include "src/algo.h"

#define NTEST 500

int main() {
    srand(42);
    Particle particles[NPARTICLES];
    double sumdiff = 0, sumdiff2 = 0, sumdiff_sqrtht = 0, sumdiff2_sqrtht = 0, sum_sqrtht = 0;
    unsigned int nfail = 0;
    for (int test = 1; test <= NTEST; ++test) {
        float ht = 0;
        for (int i = 0; i < NPARTICLES; ++i) {
            // generate an exponentially falling pt spectrum staring at 10 GeV and falling
            float pt = 10.f - 40.f*std::log((rand()+1.f)/RAND_MAX);
            ht += pt;
            particles[i].hwPt  = int(std::min<float>(pt*4,32767));
            particles[i].hwPhi = int((rand()/float(RAND_MAX/2) - 1)*M_PI / 0.01);
        }
        // run the algorithm
        pt_t hw, ref, flt_ref; pxy_t hw_px, hw_py, ref_px, ref_py, flt_refpx, flt_refpy;
        algo_main(particles, hw_px, hw_py, hw);
        algo_main_ref(particles, ref_px, ref_py, ref);
        algo_main_ref_float(particles, flt_refpx, flt_refpy, flt_ref);

        float dmet = 0.25*int(flt_ref) - 0.25*int(ref);
        sumdiff  += dmet;
        sumdiff2 += dmet*dmet;
        sumdiff_sqrtht  += dmet/std::sqrt(ht);
        sumdiff2_sqrtht += dmet*dmet/ht;
        sum_sqrtht += 0.25*int(flt_ref)/std::sqrt(ht);

        // check the output
        if (hw != ref) {
            nfail++;
            if (nfail > 3) continue;
            printf("Error in test %d\n", test);
            for (int i = 0; i < NPARTICLES; ++i) {
                printf("   Particle %2d hwPt %8d  hwPhi %+4d\n", i, 
                            int(particles[i].hwPt), int(particles[i].hwPhi));
            }
            printf("  Algo result     :  %6d   (x %+6d  y %+6d)\n", int(hw), int(hw_px), int(hw_py));
            printf("  Reference result:  %6d   (x %+6d  y %+6d)\n", int(ref), int(ref_px), int(ref_py));
            printf("\n");
            //return 1;
        }
        //printf("MET X %.3f  %.3f Y %.3f  %.3f   MET %.3f  %.3f  HT %.2f\n", 0.0625*int(flt_refpx), 0.0625*int(ref_px), 0.0625*int(flt_refpy), 0.0625*int(ref_py), 0.25*int(flt_ref), 0.25*int(ref), ht);
    }
    printf("Average MET(int) - MET(flt) = %9.1f\n", sumdiff/NTEST);
    printf("RMS     MET(int) - MET(flt) = %9.1f\n", std::sqrt(std::abs((sumdiff2/NTEST-std::pow(sumdiff/NTEST,2)))));
    printf("Average       MET(flt)          / sqrt(HT) = %9.3f\n", sum_sqrtht/NTEST);
    printf("Average [ MET(int) - MET(flt) ] / sqrt(HT) = %9.3f\n", sumdiff_sqrtht/NTEST);
    printf("RMS     [ MET(int) - MET(flt) ] / sqrt(HT) = %9.3f\n", std::sqrt(std::abs((sumdiff2_sqrtht/NTEST-std::pow(sumdiff_sqrtht/NTEST,2)))));
    if (nfail == 0) {
        printf("Passed all %d tests\n", NTEST);
        return 0;
    } else {
        printf("Failed %d / %d tests\n", nfail, NTEST);
        return 1;
    }
}
