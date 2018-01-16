#include "pluto.h"

void Resistive_eta(double *v, double x1, double x2, double x3, double *J, double *eta)
{
  double eta0 = 4*CONST_PI/(CONST_c*CONST_c)*UNIT_VELOCITY*UNIT_LENGTH;/*unit of eta for adimensionalization*/

  eta[IDIR] = g_inputParam[ETAX_GAU] / eta0;
  eta[JDIR] = g_inputParam[ETAY_GAU] / eta0;
  eta[KDIR] = g_inputParam[ETAZ_GAU] / eta0;
}
