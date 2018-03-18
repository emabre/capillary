#ifndef FREEZE_FLUID_H
#define FREEZE_FLUID_H


void ZeroHypFlux (const State_1D *state, int beg, int end, 
                 double *cmax, Grid *grid);

#endif