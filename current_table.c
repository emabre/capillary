#include "pluto.h"
#include "current_table.h"

/*Electrical current table, in Ampere, seconds*/
double i_curr_tab[N_CURR_TAB] = { 0.000000e+00,
                                  2.400000e+01,
                                  9.500000e+01,
                                  8.500000e+01,
                                  4.500000e+01,
                                  1.000000e+01,
                                  2.700000e+00,
                                  0.000000e+00};
double t_curr_tab[N_CURR_TAB] = { -1.000000e-20,
                                  1.200000e-07,
                                  2.300000e-07,
                                  3.500000e-07,
                                  5.500000e-07,
                                  9.000000e-07,
                                  1.400000e-06,
                                  1.900000e-06};
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

  if (xpoint<x[0] || xpoint>x[N-1]) print1("[Ema]Error! Data ouside table!");
  while(x[i]<xpoint) i++;
  return y[i-1] + (y[i]-y[i-1])/(x[i]-x[i-1])*(xpoint-x[i-1]);
}
