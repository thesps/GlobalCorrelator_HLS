#ifndef ALGO_DATA_H
#define ALGO_DATA_H

#include "ap_int.h"
#include "ap_fixed.h"

typedef ap_ufixed<16,14>  pt_t;      // 1 unit = 0.25 GeV; max = 16 TeV
typedef ap_fixed<10,4>    etaphi_t;   // 1 unit = 0.01;     max = +/- 5.12
typedef ap_fixed<11,5>    detaphi_t;  // 1 unit = 0.01;     max = +/- 10.24 
typedef ap_fixed<22,16> pt_etaphi_t; // type for product of pt with eta or phi
typedef ap_uint<5> count_t; // type for multiplicity

// constants for the axis update
typedef ap_ufixed<18,-2> inv_pt_t;
static constexpr int N_table_inv_pt = 1024;

static const detaphi_t TWOPI = 3.14159 * 0.78125 * 2; // 0.78125 is 100 / 128
static const detaphi_t PI = 3.14159 * 0.78125; // 0.78125 is 100 / 128
static const detaphi_t HALFPI = 3.14159 * 0.78125 / 2; // 0.78125 is 100 / 128
static const detaphi_t RCONE = 0.4 * 100 / 128;
static const detaphi_t R2CONE = RCONE * RCONE;

static const etaphi_t FIDUCIAL_ETA_PHI = 5.11 * 100 / 128;
static const pt_t JET_PT_CUT = 5;


template<class pt_T, class etaphi_T>
class TemplateParticle {
public:
    pt_T hwPt;
    etaphi_T hwEta; 
    etaphi_T hwPhi; 

    bool operator >= (const TemplateParticle<pt_T,etaphi_T> &b){
        return hwPt >= b.hwPt;
    }
};

template<class pt_T, class etaphi_T>
inline void clear(TemplateParticle<pt_T, etaphi_T> & p) {
    p.hwPt = 0;
    p.hwEta = 0;
    p.hwPhi = 0;
}

typedef TemplateParticle<pt_t, etaphi_t> Particle;
typedef TemplateParticle<pt_t, pt_etaphi_t> PartialParticle;

struct Jet : public Particle {
    ap_uint<5> iSeed;
    count_t nCand;
};

inline void clear(Jet & jet) {
    jet.hwPt = 0;
    jet.hwEta = 0;
    jet.hwPhi = 0;
    jet.iSeed = 0;
    jet.nCand = 0;
}

#define NPARTICLES 128
#define NJETS 12

//#define RCONE 40
//#define R2CONE (RCONE*RCONE)

//#define JET_PT_CUT 20
//#define FIDUCIAL_ETA_PHI 511

#endif
