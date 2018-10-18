#include "pluto.h"
#include "gamma_transp.h"
#include "current_table.h"
#include "capillary_wall.h"
#include "transport_tables.h"

#define RESMAX_PLASMA 1.0e-9

void Resistive_eta(double *v, double x1, double x2, double x3, double *J, double *eta)
{
  static int res_tab_not_done = 1;
  double mu=0.0, z=0.0, T=0.0;
  double res=0.0;
  double const res_copper = 7.8e-18; // Roughly: resisitivity of warm copper
  double const res_wall = 1.0e-7; // Roughly: resistivity of glass at 1000-2000°C
  // double unit_Mfield;

  if (res_tab_not_done) {
    MakeElecResistivityTable();
    res_tab_not_done = 0;
  }

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
    z = fmax(1/mu - 1, IONIZMIN);

    // if (g_time>1.65e-1) {
    //   print1("\n--------\n");
    //   print1("\nz:%g\n",z);
    //   print1("\nv[RHO]*UNIT_DENSITY:%g\n",v[RHO]*UNIT_DENSITY);
    //   print1("\nmu:%g\n",mu);
    //   print1("\nres:%g\n",res);      
    // }

    // unit_Mfield = COMPUTE_UNIT_MFIELD(UNIT_VELOCITY, UNIT_DENSITY);
    // res = elRes_norm(z, v[RHO]*UNIT_DENSITY, T*CONST_kB, 1, v[iBPHI]*unit_Mfield);
    res = elRes_norm_DUED(z, v[RHO]*UNIT_DENSITY, T*CONST_kB);

    #ifdef RESMAX_PLASMA
      // This is a test limiter
      if (res > RESMAX_PLASMA) {
        res = RESMAX_PLASMA;
      }
    #endif

    // I check if I am on the wall or electrode (or in an intermediate region where I smooth the res)
    if (x2>zcap_real-dzcap_real+dzcap_real*0.1 && x2<=zcap_real && x1>=rcap_real)
      res = 0.5*(res_copper+res);
    else if (x2<zcap_real-dzcap_real-dzcap_real*0.1 && x1>=rcap_real)
      res = 0.5*(res_wall+res);
    else if (x2>=zcap_real-dzcap_real-dzcap_real*0.1 && x2<=zcap_real-dzcap_real+dzcap_real*0.1 && x1>=rcap_real)
      res = (res_copper + res_wall + res)/3;
  }

  // print1("\nT:%g, z:%g, elRes:%g, eta0:%g", T, z, res, eta0);

  /***************************************************/
  /* [Ema] adimensionalization (it should be correct, I didn't change it)*/
  /**************************************************/
  eta[IDIR] =  res / UNIT_ETA;
  eta[JDIR] =  res / UNIT_ETA;
  eta[KDIR] =  res / UNIT_ETA;
  /**************************************************/
}
