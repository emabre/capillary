#include "pluto.h"
#include "current_table.h"

/*Electrical current table, in Ampere, seconds*/
// This is the 95A peak current discharge (same as jets_new/base1/jet/ but converted current units to Ampere)
double i_curr_tab[N_CURR_TAB] = { 0.000000e00,
                                  4.0,
                                  24.0,
                                  95.0,
                                  85.0,
                                  45.0,
                                  10.0,
                                  2.7,
                                  0.0};

double t_curr_tab[N_CURR_TAB] = {-1.000000e-20,
                                  1.500000e-07,
                                  2.200000e-07,
                                  3.300000e-07,
                                  4.500000e-07,
                                  6.500000e-07,
                                  1.000000e-06,
                                  1.500000e-06,
                                  2.000000e-06};
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
