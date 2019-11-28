#ifndef ALGO_H
#define ALGO_H

#include "data.h"

// implementation to be synthethised
void algo_main(const Particle particles[NPARTICLES], pxy_t & met_px, pxy_t & met_py, pt_t & met_pt) ;

// reference implementation for validation
void algo_main_ref(const Particle particles[NPARTICLES], pxy_t & met_px, pxy_t & met_py, pt_t & met_pt) ;

// reference implementation for validation
void algo_main_ref_float(const Particle particles[NPARTICLES], pxy_t & met_px, pxy_t & met_py, pt_t & met_pt) ;

#endif
