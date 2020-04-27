#ifndef ALGO_DATA_H
#define ALGO_DATA_H

#include "ap_int.h"

typedef ap_uint<16>  pt_t;      // 1 unit = 0.25 GeV; max = 16 TeV
typedef ap_int<9>    etaphi_t;  // 1 unit = 0.01;     max = +/- 2.55
// we're assuming anyway that all the particles are in a region of smallish size,
// so that there's no need to wrap-aroud the phi coordinate


struct Particle {
    pt_t hwPt;
    etaphi_t hwEta; 
    etaphi_t hwPhi; 
};
inline void clear(Particle & p) {
    p.hwPt = 0;
    p.hwEta = 0;
    p.hwPhi = 0;
}

struct Jet : public Particle {
    ap_uint<5> iSeed;
    ap_uint<5> nCand;
};

inline void clear(Jet & jet) {
    jet.hwPt = 0;
    jet.hwEta = 0;
    jet.hwPhi = 0;
    jet.iSeed = 0;
    jet.nCand = 0;
}

#define NPARTICLES 128
#define NJETS 8
#define MOREJETS 2

#define RCONE 40
#define R2CONE (RCONE*RCONE)

#define JET_PT_CUT 100
#define FIDUCIAL_ETA_PHI 125

#endif
