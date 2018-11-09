#include "pluto.h"
#include "gamma_transp.h"
#include "current_table.h"
#include "capillary_wall.h"
#include "transport_tables.h"

#define RESMAX_PLASMA 1.0e-9
#define REALISTIC_WALL_ETA NO

void Resistive_eta(double *v, double x1, double x2, double x3, double *J, double *eta)
{
  #if ETA_TABLE
  static int res_tab_not_done = 1;
  #else
  double mu=0.0, z=0.0;
  #endif
  double T=0.0;
  double res=0.0;
  #if REALISTIC_WALL_ETA
    double const res_copper = 7.8e-18; // Roughly: resisitivity of warm copper
    double const res_wall = 1.0e-7; // Roughly: resistivity of glass at 1000-2000°C
  #endif
  // double unit_Mfield;

  if (g_inputParam[ETAX_GAU] > 0.0) {
    // Fixed value from pluto.ini
    res = g_inputParam[ETAX_GAU];
    
  } else {
    if (GetPV_Temperature(v, &(T) )!=0) {
      #if WARN_ERR_COMP_TEMP
        print1("\nResistive_eta:[Ema]Err.comp.temp");
      #endif
    }
    #if ETA_TABLE
      if (res_tab_not_done) {
        MakeElecResistivityTable();
        res_tab_not_done = 0;
      }
      res = GetElecResisitivityFromTable(v[RHO]*UNIT_DENSITY, T);
    #else
      GetMu(T, v[RHO], &mu);
      z = fmax(1/mu - 1, IONIZMIN);
      res = elRes_norm_DUED(z, v[RHO]*UNIT_DENSITY, T*CONST_kB);
    #endif

    #ifdef RESMAX_PLASMA
      // This is a test limiter
      if (res > RESMAX_PLASMA) {
        res = RESMAX_PLASMA;
      }
    #endif

    #if REALISTIC_WALL_ETA
      // I check if I am on the wall or electrode (or in an intermediate region where I smooth the res)
      if (x2>zcap_real-dzcap_real+dzcap_real*0.1 && x2<=zcap_real && x1>=rcap_real)
        res = 0.5*(res_copper+res);
      else if (x2<zcap_real-dzcap_real-dzcap_real*0.1 && x1>=rcap_real)
        res = 0.5*(res_wall+res);
      else if (x2>=zcap_real-dzcap_real-dzcap_real*0.1 && x2<=zcap_real-dzcap_real+dzcap_real*0.1 && x1>=rcap_real)
        res = (res_copper + res_wall + res)/3;
    #endif
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
