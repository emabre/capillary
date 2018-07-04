#include "pluto.h"
#include "current_table.h"

/*Electrical current table, in Ampere, seconds*/
// Current as in commit 5e85a9c57...
double i_curr_tab[N_CURR_TAB] = { 0.0001e+00,
                                  200.00e+00,
                                  200.00e+00};
double t_curr_tab[N_CURR_TAB] = { -1.000000e-20,
                                  200.0e-9,
                                  2000.000000e-09};
//
/*------------------------------------------------------
  Function to compute current at the present time step
------------------------------------------------------*/
double current_from_time(double t){
  return interp_lin_table(t, t_curr_tab, i_curr_tab, N_CURR_TAB);
}


/*------------------------------------------------------
  Function computing something by linear-interpolating values from a table
  Note: the table is assumed to be ordered in ascending order (x[i]<x[i+1])
--------------------------------------------------------*/
double interp_lin_table(double xpoint, double *x, double *y, const int N){
  int i = 0;

  if (xpoint<x[0] || xpoint>x[N-1]) {
    print1("\n[Ema]Error! Data ouside table!");
    print1("\nxpoint=%d, x[0]=%d", xpoint, x[0]);
    QUIT_PLUTO(1);
  }

  while(x[i]<xpoint) i++;
  return y[i-1] + (y[i]-y[i-1])/(x[i]-x[i-1])*(xpoint-x[i-1]);
}
