#include "pluto.h"
#include "gamma_transp.h"
#include "current_table.h"
#include "transport_tables.h"
#include "capillary_wall.h"

#define KAPPAMAX 1e7
#define KAPPA_LOW 1e3
#define LOW_RHO_REL_TOLL 1e-2

void TC_kappa(double *v, double x1, double x2, double x3,
              double *kpar, double *knor, double *phi)
{
  #if KAPPA_TABLE
  static int tc_tab_not_done = 1;
  #else
  double mu=0.0, z=0.0;
  #endif
  double T=0.0;
  double k=0.0;
  // double unit_Mfield;

  if (g_inputParam[KAPPA_GAU] > 0.0) {
    // Fixed value from pluto.ini
    *kpar = g_inputParam[KAPPA_GAU];
    *knor = g_inputParam[KAPPA_GAU];

  } else {
    if (GetPV_Temperature(v, &(T) )!=0) {
      print1("\nTC_kappa:[Ema]Err.comp.temp");
    }
    #if KAPPA_TABLE
      if (tc_tab_not_done) {
        MakeThermConductivityTable();
        tc_tab_not_done = 0;
      }
      k = GetThermConductivityFromTable(v[RHO]*UNIT_DENSITY, T);
    #else
      GetMu(T, v[RHO], &mu);
      z = fmax(1/mu - 1, IONIZMIN);

      // unit_Mfield = COMPUTE_UNIT_MFIELD(UNIT_VELOCITY, UNIT_DENSITY);
      // k = thermCond_norm(z, v[RHO]*UNIT_DENSITY, T*CONST_kB, 1, v[iBPHI]*unit_Mfield);
      k = thermCond_norm_DUED(z, v[RHO]*UNIT_DENSITY, T*CONST_kB);
      // I increase the value of k to take into account the conductivity of a hydrogen gas at 2000K
      // (according to Timrot value of H2 conductivity, cited in Mehl et al., "Ab initio transport...",2010 ),
      // anyway, I am not sure this correction is ok, especially because I don't know what is the "low-density
      // limit" where this value is applicable. Moreover, I think it might be that I should also add the
      // conductivity of hydrogen ions (negligible at high T).
      k = k + 8e4;
    #endif

    #ifdef KAPPAMAX
      if (k > KAPPAMAX) {
        k = KAPPAMAX;
      }
    #endif

    #ifdef CONE_LOW_TCKAPPA
      if (IsOutCone(CONE_LOW_TCKAPPA, x1, x2) && v[RHO] < LOW_RHO_REL_TOLL*g_inputParam[DENS0]/UNIT_DENSITY)
        k = KAPPA_LOW;
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
  *kpar /= UNIT_KAPPA;
  *knor /= UNIT_KAPPA;
  /***************************************************/

  *phi = 1; //[Ema] this should be useless with the saturation of the thermal conduction off

  if (*kpar>1e15 || *knor>1e15) {
      print1("\nDid you take the wrong convention for kappa?\n");
      QUIT_PLUTO(1);
  }

}
