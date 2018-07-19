#include "pluto.h"
#include "current_table.h"

/*Electrical current table, in Ampere, seconds*/
double i_curr_tab[N_CURR_TAB];
double t_curr_tab[N_CURR_TAB];
int Ncurr;
// double i_curr_tab[N_CURR_TAB] = { 0.0001e+00,
//                                   9.600000e+01,
//                                   1.700000e+02,
//                                   1.400000e+02,
//                                   6.700000e+01,
//                                   1.200000e+01};
// double t_curr_tab[N_CURR_TAB] = { -1.000000e-20,
//                                   7.500000e-08,
//                                   2.500000e-07,
//                                   4.650000e-07,
//                                   7.260000e-07,
//                                   1.260000e-06};
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
  static int first_call=1;
  int i;
  char filename[] = "current_table.dat";

  // I initialize the table reading it from file, if this is the first call
  if (first_call) {
    // Initialize to zero
    for (i=0; i<N_CURR_TAB; i++) {
      t_curr_tab[i] = 0.0;
      i_curr_tab[i] = 0.0;
    }
    if(read_table_file(t_curr_tab, i_curr_tab, &Ncurr, N_CURR_TAB, filename)) {
      print1("Error reading electric current table file: %s\n", filename);
      QUIT_PLUTO(1);
    }
    first_call = 0;
  }

  return interp_lin_table(t, t_curr_tab, i_curr_tab, Ncurr);
}


/*------------------------------------------------------
  Function computing something by linear-interpolating values from a table
  Note: the table is assumed to be ordered in ascending order (x[i]<x[i+1])
--------------------------------------------------------*/
double interp_lin_table(double xpoint, double *x, double *y, const int N){
  int i = 0;

  if (xpoint<x[0] || xpoint>x[N-1]) {
    print1("\n[interp_lin_table]Error! Data ouside table!");
    print1("\nx[0]=%g, xpoint=%g, x[N-1]=%g", x[0], xpoint, x[N-1]);
    QUIT_PLUTO(1);
  }

  while(x[i]<xpoint) i++;
  return y[i-1] + (y[i]-y[i-1])/(x[i]-x[i-1])*(xpoint-x[i-1]);
}

/*------------------------------------------------------
  Function to read a 2xN table from file
------------------------------------------------------*/
int read_table_file(double *x, double *y, int *N, int Nmax, char *filename) {
  FILE *table;
  char header1[300], header2[300];
  int i;
  int n_read;

  table = fopen("current_table.dat", "r");

  if (table == NULL) {
    print1("\nError reading current table file\n");
  }

  // fscanf(table,"%s \n", header);
  fscanf(table,"%[^\n] \n", header1);
  fscanf(table,"%[^\n] \n", header2);

  #ifdef DEBUG_EMA
    printf("\ntable %s header1: %s", filename, header1);
    printf("\ntable %s header2: %s", filename, header2);
  #endif

  i = 0;
  while( 2 == (n_read=fscanf(table, "%lf %lf", &(x[i]), &(y[i]))) ) {
    if (i>=Nmax) {
      print1("\n[read_table_file] More entry in table then allowed!");
      return -1;
    }
    if (n_read==1 || n_read>2) {
      print1("\n[read_table_file] Wrong number of entries in table line!");
      print1("\nentries: %d", n_read);
      return n_read;
    }
    i++;
  }

  *N = i;

  #ifdef DEBUG_EMA
    printf("\nI read %d lines in table", *N);
    for (i=0; i<*N; i++) {
      printf("\nx:%e,y:%e", x[i], y[i]);
    }
  #endif

  return 0;
}