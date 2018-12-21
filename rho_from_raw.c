#include "pluto.h"

#define INITIAL_RHO_SCRIPT  "initial_mass_dens/initial_condition_4PLUTO.py"
#define NCELLS_RHO_X1  1000
#define NCELLS_RHO_X2  2000

/*******************************************************************************
  Function to write a binary file with rho and a ascii file with the corresponding grid,
  (from which later pluto will interpolate),
  starting from a ASCII file with pressure and a file with temperature.
*******************************************************************************/
void WriteRhoGridFromRaw(char rho_ic_fi[], char grid_rho_ic_fi[]){

  char command[300];

  sprintf(command, "python3 %s %e %e %d %d %e %e %e %e %e %d %s %s", INITIAL_RHO_SCRIPT,
          (double)(UNIT_LENGTH), (double)(UNIT_DENSITY),
          (int) (NCELLS_RHO_X1), (int)(NCELLS_RHO_X2),
          (double)(UNIT_LENGTH*(g_domBeg[IDIR] - 0.001*fabs(g_domBeg[IDIR]))),
          (double)(UNIT_LENGTH*(g_domEnd[IDIR] + 0.001*fabs(g_domEnd[IDIR]))),
          (double)(UNIT_LENGTH*(g_domBeg[JDIR] - 0.001*fabs(g_domBeg[JDIR]))),
          (double)(UNIT_LENGTH*(g_domEnd[JDIR] + 0.001*fabs(g_domEnd[JDIR]))),
          (double)((RHO_TAB_MIN)*1.002),
          (int) REPLOT_P_RHO,
          (char*)(rho_ic_fi), (char*)(grid_rho_ic_fi) );
  system(command);
}