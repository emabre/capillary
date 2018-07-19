#include "pluto.h"
#include "current_table.h"

/*Electrical current table, in Ampere, seconds*/
double i_curr_tab[N_CURR_TAB] = { 0.0001e+00,
                                  9.600000e+01,
                                  1.700000e+02,
                                  1.400000e+02,
                                  6.700000e+01,
                                  1.200000e+01};
double t_curr_tab[N_CURR_TAB] = { -1.000000e-20,
                                  7.500000e-08,
                                  2.500000e-07,
                                  4.650000e-07,
                                  7.260000e-07,
                                  1.260000e-06};
// Current as in commit 5e85a9c57...
// double i_curr_tab[N_CURR_TAB] = { 0.0001e+00,
//                                   200.00e+00,
//                                   200.00e+00};
// double t_curr_tab[N_CURR_TAB] = { -1.000000e-20,
//                                   200.0e-9,
//                                   2000.000000e-09};
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
