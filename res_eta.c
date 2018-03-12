#include "pluto.h"
#include "gamma_transp.h"
#include "current_table.h"

void Resistive_eta(double *v, double x1, double x2, double x3, double *J, double *eta)
{
  const double eta0 = 4*CONST_PI/(CONST_c*CONST_c)*UNIT_VELOCITY*UNIT_LENGTH;/*unit of eta for adimensionalization*/
  double mu, z, T, unit_Mfield;
  double res;

  unit_Mfield = COMPUTE_UNIT_MFIELD(UNIT_VELOCITY, UNIT_DENSITY);

  if (GetPV_Temperature(v, &(T) )!=0) {
    print1("Resistive_eta:[Ema] Error computing temperature!");
  }
  // print1("\nI just assigned %g to T[%d][%d][%d] for output",T[k][j][i], k,j,i);
  GetMu(T, v[RHO], &mu);
  z = 1/mu - 1;

  // res = elRes_norm(z, v[RHO]*UNIT_DENSITY, T*CONST_kB, 1, v[iBPHI]*unit_Mfield);
  res = elRes_norm_DUED(z, v[RHO]*UNIT_DENSITY, T*CONST_kB);

  // print1("\nT:%g, z:%g, elRes:%g, eta0:%g", T, z, res, eta0);

  eta[IDIR] =  res / eta0;
  eta[JDIR] =  res / eta0;
  eta[KDIR] =  res / eta0;
  
  // To set a fixed eta from input pluto.ini
  //eta[IDIR] = g_inputParam[ETAX_GAU] / eta0;
  //eta[JDIR] = g_inputParam[ETAY_GAU] / eta0;
  //eta[KDIR] = g_inputParam[ETAZ_GAU] / eta0;

}
