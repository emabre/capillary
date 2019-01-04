/* /////////////////////////////////////////////////////////////////// */
/*! \file  
 *  \brief Specification of explicit first and second viscosity coefficients*/
/* /////////////////////////////////////////////////////////////////// */
#include "pluto.h"
/* ************************************************************************** */
void Visc_nu(double *v, double x1, double x2, double x3,
                        double *nu1, double *nu2)
/*! 
 *
 *  \param [in]      v  pointer to data array containing cell-centered quantities
 *  \param [in]      x1 real, coordinate value 
 *  \param [in]      x2 real, coordinate value 
 *  \param [in]      x3 real, coordinate value 
 *  \param [in, out] nu1  pointer to first viscous coefficient
 *  \param [in, out] nu2  pointer to second viscous coefficient
 *
 *  \return This function has no return value.
 * ************************************************************************** */
{
  
  *nu2 = 0.0;
  if (g_inputParam[DYN_VISC_GAU] > 0.0) {
    // Fixed value from pluto.ini
    *nu1  = g_inputParam[DYN_VISC_GAU]/(UNIT_DENSITY*UNIT_LENGTH*UNIT_VELOCITY);
    // MAGARI FARE CHE SE LA RHO SCENDE SOTTO UN CERTO VALORE, ALLORA LA nu1 VIENE RIDOTTA PROPORZIONALMENTE?
    // MOTIVO: COSÌ NON SI ALZA TROPPO IL PARAMETRO DI DIFFUSIONE, TANTO, CHI SE NE CURA DI MODELLARE BENE LE ZONE A DENSITÀ MOLTO BASSA!!
  } else {
    print1("Accurate viscosity not implemented yet");
    QUIT_PLUTO(1);
  }

}
