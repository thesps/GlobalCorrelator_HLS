#ifndef ALGO_H
#define ALGO_H

#include "data.h"
#include "hls_stream.h"

// implementation to be synthethised
void seedcone(const Particle particles[NPARTICLES], Jet jet[NJETS]) ;

// reference implementation for validation
void algo_main_ref(const Particle particles[NPARTICLES], Jet jet[NJETS]) ;

void iterate(hls::stream<bool> &reset, hls::stream<bool> &newInput, hls::stream<Particle> inParticles[NPARTICLES], hls::stream<Jet> outJets[NJETS]);

#endif
