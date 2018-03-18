#include"pluto.h"
#include "freeze_fluid.h"

/* ********************************************************************* */
void ZeroHypFlux (const State_1D *state, int beg, int end, 
                 double *cmax, Grid *grid)
{
/*!
 * This function is a replacement for a Riemann solver in case one wants
 * to keep the fluid evolution frozen (no advection at all).
 * I created this function modifying the hll Riemann solver.
 * 
 * \param[in,out] state   pointer to State_1D structure
 * \param[in]     beg     initial grid index
 * \param[out]    end     final grid index
 * \param[out]    cmax    1D array of maximum characteristic speeds
 * \param[in]     grid    pointer to array of Grid structures.
 *
 *********************************************************************** */
{
  int    nv, i, xdface;
  double scrh;
  double *vL, *vR, *uL, *uR, *SR, *SL;
  static double **VL, **VR, **UL, **UR;
  static double *a2L, *a2R;
  double **bgf;
  /* **fL is useless, but it is used as check to know whehter
  other variables have to be allocate.
  Since I do not want to loose time debugging, I keep it.
  */
  static double **fL;
  // static double **fL, **fR;
  // static double *pL, *pR;
  // static double **Uhll;
    
  if (fL == NULL){
    fL = ARRAY_2D(NMAX_POINT, NFLX, double);
    // fR = ARRAY_2D(NMAX_POINT, NFLX, double);

    //pL  = ARRAY_1D(NMAX_POINT, double);
    //pR  = ARRAY_1D(NMAX_POINT, double);
    a2L = ARRAY_1D(NMAX_POINT, double);
    a2R = ARRAY_1D(NMAX_POINT, double);

    // Uhll = ARRAY_2D(NMAX_POINT, NFLX, double);
  }

/* ------------------------------------------------
     Solve 2x2 Riemann problem with GLM Cleaning
   ------------------------------------------------ */
   VL = state->vL; UL = state->uL;
   VR = state->vR; UR = state->uR;

/* ----------------------------------------------------
     compute sound speed & fluxes at zone interfaces
   ---------------------------------------------------- */

  SoundSpeed2 (VL, a2L, NULL, beg, end, FACE_CENTER, grid);
  SoundSpeed2 (VR, a2R, NULL, beg, end, FACE_CENTER, grid);

  // Next two lines to be deleted
  //Flux (UL, VL, a2L, bgf, fL, pL, beg, end);
  //Flux (UR, VR, a2R, bgf, fR, pR, beg, end);

/* ----------------------------------------
     get max and min signal velocities
   ---------------------------------------- */
             
  SL = state->SL; SR = state->SR;
  HLL_Speed (VL, VR, a2L, a2R, bgf, SL, SR, beg, end);

/* ----------------------------------------
           compute HLL flux
   ---------------------------------------- */	     

  for (i = beg; i <= end; i++) {

    scrh  = MAX(fabs(SL[i]), fabs(SR[i]));
    cmax[i] = scrh;
  
    for (nv = 0; nv < NFLX; nv++) {
        state->flux[i][nv] = 0.0;
      }
    state->press[i] = 0.0;
/* Following lines to be deleted...
    if (SL[i] > 0.0){
    
      for (nv = 0; nv < NFLX; nv++) {
        state->flux[i][nv] = fL[i][nv];
      }
      state->press[i] = pL[i];
      
    }else if (SR[i] < 0.0){
    
      for (nv = 0; nv < NFLX; nv++) {
        state->flux[i][nv] = fR[i][nv];
      }
      state->press[i] = pR[i];
      
    }else{
    
      uL = UL[i]; uR = UR[i];

      scrh = 1.0 / (SR[i] - SL[i]);
    
      for (nv = 0; nv < NFLX; nv++) {
        state->flux[i][nv] = SL[i]*SR[i]*(uR[nv] - uL[nv]) +
                             SR[i]*fL[i][nv] - SL[i]*fR[i][nv];
        state->flux[i][nv] *= scrh;
      }
      state->press[i] = (SR[i]*pL[i] - SL[i]*pR[i])*scrh;
    }
*/
  }

}