#include "algo.h"
#include <cmath>
#include <hls_math.h>
#include <cassert>
#ifndef __SYNTHESIS__
#include <cstdio>
#endif

//=== fill in a table with sin[phi] and cos[phi] * 2^8 = 256 
//--- version with two lookup tables
void _lut_sincos_init(sincos_t table_sin[512], sincos_t table_cos[512]) {
//--- version with one lookup table includign both sets of numbers
// void _lut_sincos_init(ap_uint<18> table_sincos[512]) {
    for (int i = 0; i < 512; ++i) {
        float alpha = 0.01f*i;
        int is = sin(alpha) * 256, ic = cos(alpha) * 256;
        if (is > 255) is = 255; if (is < -255) is = -255;
        if (ic > 255) ic = 255; if (ic < -255) ic = -255;
        //-- two BRAMs
        table_sin[i] = is; table_cos[i] = ic;
        //-- one BRAM
        //table_sincos[i] = ap_int<9>(is).concat(ap_int<9>(ic));
    }
}

void pt_phi_to_px_py(pt_t pt, etaphi_t phi, pxy_t & px, pxy_t & py) {
    // --- version with 2 BRAMs
    sincos_t _table_sin[512], _table_cos[512]; _lut_sincos_init( _table_sin, _table_cos );
    // --- version with 1 BRAM
    // ap_uint<18> _table_sincos[512]; _lut_sincos_init( _table_sincos );
    int iphi = phi; 
    if (phi < 0) iphi = -iphi;
    // --- version with 2 BRAMs
    sincos_t sinphi = _table_sin[iphi], cosphi = _table_cos[iphi];
    // --- version with 1 BRAM
    // ap_uint<18> record = _table_sincos[iphi];
    // sincos_t sinphi = sincos_t(record(17,9)), cosphi = sincos_t(record(8, 0));
    if (phi < 0) sinphi = -sinphi;
    px = (pt * cosphi) >> 6; // throw away some bits already here
    py = (pt * sinphi) >> 6; // throw away some bits already here
}

#define SQRT_TABLE_SIZE 4*1024
void _lut_sqrt_init(ap_uint<6> table[SQRT_TABLE_SIZE]) {
    for (int i = 0; i < SQRT_TABLE_SIZE; ++i) {
        table[i] = int(std::sqrt(float(i)));
    }
}
int _short_sqrt(int x) {
    ap_uint<6> sqrt_table[SQRT_TABLE_SIZE];
    _lut_sqrt_init(sqrt_table);
    assert(x >= 0 && x < SQRT_TABLE_SIZE);
    int ret = sqrt_table[x];
    return ret;
}
pt_t px_py_to_pt(int px, int py) {

    ap_int<13> pxGeV = px >> 4;
    ap_int<13> pyGeV = py >> 4;
    ap_uint<25> pt2 = (pxGeV * pxGeV) + (pyGeV * pyGeV); // 25 bits, always positive

    pt_t ret; 
    if (pt2.range(24,18).or_reduce()) {//pt2 >= 512*512) {
        ret = pt_t(512 << 2);
    } else if (pt2[17] || pt2[16]) {// pt2 >=  64*1024) {
        ret = pt_t(_short_sqrt(pt2 >> 6) << (2+3));
    } else if (pt2[15] || pt2[14]) {// pt2 >=  16*1024) {
        ret = pt_t(_short_sqrt(pt2 >> 4) << (2+2));
    } else if (pt2[13] || pt2[12]) {// pt2 >=   4*1024) {
        ret = pt_t(_short_sqrt(pt2 >> 2) << (2+1));
    } else {
        ret = pt_t(_short_sqrt(   pt2  ) <<   2  );
    }

    return ret;
}
void algo_main(const Particle particles[NPARTICLES], pxy_t & met_px, pxy_t & met_py, pt_t & met_pt) {
    #pragma HLS INTERFACE ap_none port=met_px
    #pragma HLS INTERFACE ap_none port=met_py
    #pragma HLS INTERFACE ap_none port=met_pt
    #pragma HLS ARRAY_PARTITION variable=particles complete
    #pragma HLS pipeline II=1
    int sumx = 0, sumy = 0;
    for (unsigned int i = 0; i < NPARTICLES; ++i) {
        pxy_t px, py;
        pt_phi_to_px_py(particles[i].hwPt, particles[i].hwPhi, px, py);
        sumx += px;
        sumy += py;
    }
    met_px = sumx;
    met_py = sumy;
    met_pt = px_py_to_pt(sumx,sumy);
}

