#ifndef ALGO_TAU_LAYER2_H
#define ALGO_TAU_LAYER2_H

#include <ap_int.h>
#include <hls_stream.h>
#include "data.h"

#define DATA_SIZE 64
#define NTAU  6
#define NREGIONS 36
#define NPART 25
#define DEPTH 2
#define NTAUPARTS  10
#define DRCONE 8410
#define DR2MAX 10000
#define MP7_NCHANN 72
//typedef ap_uint<32> MP7DataWord;

static float PT_SCALE = 4.0;     // quantize in units of 0.25 GeV (can be changed)
static float ETAPHI_FACTOR = 4;  // size of an ecal crystal in phi in integer units (our choice)
static float ETAPHI_SCALE = ETAPHI_FACTOR*(180./M_PI);  // M_PI/180 is the size of an ECal crystal; we make a grid that is 4 times that size
static int16_t PHI_WRAP = 360*ETAPHI_FACTOR;            // what is 3.14 in integer

void algo_tau_layer2_v3(MP7DataWord input[MP7_NCHANN],MP7DataWord output[MP7_NCHANN]);
//void algo_tau_layer2(hls::stream<PFChargedObj > allparts_in [DATA_SIZE],hls::stream<axi_t> &link_out);

#endif