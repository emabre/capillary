#include "pluto.h"
#include "gamma_transp.h"

void TC_kappa(double *v, double x1, double x2, double x3,
              double *kpar, double *knor, double *phi)
{
  double mu, T, sqT;
  double gamma1Bob=1; /*Just for now, as I do not want to debug too much*/
  double ne, ioniz, freq_coll_e;
  double const q_elem4 = Q_ELEM*Q_ELEM*Q_ELEM*Q_ELEM;

  //T = 20000; //[Ema] almost 2eV
  GetPV_Temperature(v, &T);

  #if EOS==IDEAL
      mu = MeanMolecularWeight(v);
  #elif EOS==PVTE_LAW
      GetMu(T, v[RHO], &mu);
  #endif

  sqT  = sqrt(T);

  // /***************************************************/
  // // Bobrova/esaulov's formulas (so that I can compare the results with the one published)
  // ioniz = 1/mu - 1;
  // ne = ioniz * v[RHO] / CONST_mp;
  // /* Collision frequency of electrons according to Esaulov/Bobrova*/
  // freq_coll_e = 4*sqrt(2*CONST_PI)/3*q_elem4/sqrt(CONST_me)*ne/(T*T*sqT)*(cl_ei(ne,T)+cl_en(ioniz,T));
  // *knor = ne*T/(CONST_me*freq_coll_e)*gamma1Bob;
  // // This kpar expression might be uncorrect,
  // // but it is unused in the axis-symmetric approx
  // // with B = (0,Bphi,0)
  // *kpar = ne*T/(CONST_me*freq_coll_e)*gamma1Bob;
  // /***************************************************/

  /***************************************************/
  // simplified formula present in the documentation
  //*kpar = 5.6e-7*T*T*sqT;
  //*knor = 5.6e-7*T*T*sqT;
  /***************************************************/
  *kpar = g_inputParam[KAPPA_GAU];
  *knor = g_inputParam[KAPPA_GAU];
  /***************************************************/
// [Ema] adimensionalization (it should be correct, I didn't change it)
  *kpar *= CONST_mp*mu/(UNIT_DENSITY*UNIT_VELOCITY*UNIT_LENGTH*CONST_kB);
  *knor *= CONST_mp*mu/(UNIT_DENSITY*UNIT_VELOCITY*UNIT_LENGTH*CONST_kB);
  /***************************************************/

  *phi = 1; //[Ema] this should be useless with the saturation of the thermal conduction off

    //this is just a check of good compilation of an additiona file
    //gammaBob = Gamma_i(3,0.01,0.3);

}
