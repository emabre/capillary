#include "pluto.h"
#include "gamma_transp.h"
#include "current_table.h"

#define KAPPAMAX 1e24

void TC_kappa(double *v, double x1, double x2, double x3,
              double *kpar, double *knor, double *phi)
{
  double mu=0.0, z=0.0, T=0.0;
  double k=0.0;
  // double unit_Mfield;

  if (GetPV_Temperature(v, &(T) )!=0) {
    print1("\nTC_kappa:[Ema] Error computing temperature!");
  }
  // print1("\nI just assigned %g to T[%d][%d][%d] for output",T[k][j][i], k,j,i);
  GetMu(T, v[RHO], &mu);

  if (g_inputParam[KAPPA_GAUBOB] > 0.0) {
    
    // Fixed value from pluto.ini
    *kpar = g_inputParam[KAPPA_GAUBOB]*CONST_kB;
    *knor = g_inputParam[KAPPA_GAUBOB]*CONST_kB;

  } else {
    z = 1/mu - 1;

    // unit_Mfield = COMPUTE_UNIT_MFIELD(UNIT_VELOCITY, UNIT_DENSITY);
    // k = thermCond_norm(z, v[RHO]*UNIT_DENSITY, T*CONST_kB, 1, v[iBPHI]*unit_Mfield);
    k = thermCond_norm_DUED(z, v[RHO]*UNIT_DENSITY, T*CONST_kB);

    #ifdef KAPPAMAX
      if (k > KAPPAMAX) {
        k = KAPPAMAX;
      }
    #endif

    *knor = k;
    *kpar = k; //This should be useless

    // simplified formula present in the documentation
    //*kpar = 5.6e-7*T*T*sqT;
    //*knor = 5.6e-7*T*T*sqT;
  }

  /**************************************************
   [Ema] Adimensionalization. It is slightly modified from
   what is adviced in PLUTO's userguide:
   I do not normalize by mu, since it is not constant over the domain,
   I belive that it is not correct in general to normalize by mu.
   Normalizing by mu could be a good idea when the EOS is IDEAL and no
   change in mu over space/time is possible [at the time of writing it is
   indeed necessary when EOS==IDEAL, due to the way temperature is computed in
   the STS and EXPLICIT algorithms (T=p/rho). If one decides not to normalize by
   mu then he should change the computation of T for the ideal case with
   T=mu*p/rho]. For using ADI scheme I implemented, one must not normalize by mu.
   *************************************************/
  *kpar *= CONST_mp/(UNIT_DENSITY*UNIT_VELOCITY*UNIT_LENGTH*CONST_kB);
  *knor *= CONST_mp/(UNIT_DENSITY*UNIT_VELOCITY*UNIT_LENGTH*CONST_kB);
  /***************************************************/

  *phi = 1; //[Ema] this should be useless with the saturation of the thermal conduction off

  if (*kpar>1e15 || *knor>1e15) {
      print1("\nDid you take the wrong convention for kappa?\n");
      QUIT_PLUTO(1);
  }

}
