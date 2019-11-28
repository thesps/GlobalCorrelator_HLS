#ifndef ALGO_H
#define ALGO_H

#include "data.h"

// implementation to be synthethised
void algo_main(const Particle particles[NPARTICLES], Jet jet[NJETS]) ;

// reference implementation for validation
void algo_main_ref(const Particle particles[NPARTICLES], Jet jet[NJETS]) ;

#endif
