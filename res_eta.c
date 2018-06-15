#include "pluto.h"
#include "gamma_transp.h"
#include "current_table.h"
#include "capillary_wall.h"

#define RESMAX 1.0e-9

void Resistive_eta(double *v, double x1, double x2, double x3, double *J, double *eta)
{
  const double eta0 = 4*CONST_PI/(CONST_c*CONST_c)*UNIT_VELOCITY*UNIT_LENGTH;/*unit of eta for adimensionalization*/
  double mu=0.0, z=0.0, T=0.0;
  double res=0.0;
  // double unit_Mfield;

  if (g_inputParam[ETAX_GAU] > 0.0) {
    
    // Fixed value from pluto.ini
    res = g_inputParam[ETAX_GAU];

    //[Err] delete next test if lines
    // if (x2 >= zcap_real-0.5*dzcap_real && x2 <= zcap_real)
    //   res = g_inputParam[ETAX_GAU] * (1 + 100*(x2-(zcap_real-0.5*dzcap_real))/(0.5*dzcap_real));
    // else if (x2 > zcap_real)
    //   res = g_inputParam[ETAX_GAU] * 100;
    // else if (x2 < zcap_real-0.5*dzcap_real)
    //   res = g_inputParam[ETAX_GAU];
    // else {
    //   print1("[Resistive eta]Something wrong!");
    //   QUIT_PLUTO(1);
    // }

    //[Err] delete next test if lines
    // if (x2 >= zcap_real && x2 <= zcap_real+0.5*dzcap_real)
    //   res = g_inputParam[ETAX_GAU] * (1 + 100*(x2-zcap_real)/(0.5*dzcap_real));
    // else if (x2 > zcap_real+0.5*dzcap_real)
    //   res = g_inputParam[ETAX_GAU] * 100;
    // else if (x2 < zcap_real)
    //   res = g_inputParam[ETAX_GAU];
    // else {
    //   print1("[Resistive eta]Something wrong!");
    //   QUIT_PLUTO(1);
    // }
    //[Err] end test

  } else {

    if (GetPV_Temperature(v, &(T) )!=0) {
      print1("\nResistive_eta:[Ema] Error computing temperature!");
    }
    // print1("\nI just assigned %g to T[%d][%d][%d] for output",T[k][j][i], k,j,i);
    GetMu(T, v[RHO], &mu);
    z = 1/mu - 1;

    // unit_Mfield = COMPUTE_UNIT_MFIELD(UNIT_VELOCITY, UNIT_DENSITY);
    // res = elRes_norm(z, v[RHO]*UNIT_DENSITY, T*CONST_kB, 1, v[iBPHI]*unit_Mfield);
    res = elRes_norm_DUED(z, v[RHO]*UNIT_DENSITY, T*CONST_kB);
    
    #ifdef RESMAX
      // This is a test limiter
      if (res > RESMAX) {
        res = RESMAX;
      }
    #endif
  }

  // print1("\nT:%g, z:%g, elRes:%g, eta0:%g", T, z, res, eta0);

  /***************************************************/
  /* [Ema] adimensionalization (it should be correct, I didn't change it)*/
  /**************************************************/
  eta[IDIR] =  res / eta0;
  eta[JDIR] =  res / eta0;
  eta[KDIR] =  res / eta0;
  /**************************************************/

}
