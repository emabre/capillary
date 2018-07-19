#ifndef CURRENT_TABLE_H
#define CURRENT_TABLE_H

#define COMPUTE_UNIT_MFIELD(unit_vel, unit_dens) ((unit_vel)*sqrt(4*CONST_PI*(unit_dens)))

// This is a particular form of biot-savart law, to compute B in gauss,
// from current in Ampere and radius (g_domEnd[0]*UNIT_LENGTH) in cm
#define BIOTSAV_GAUSS_A_CM(curr_A, rad_cm) (0.2*(curr_A)/(rad_cm))

#define N_CURR_TAB 6

double i_curr_tab[N_CURR_TAB];  /*should I also say extern?*/
double t_curr_tab[N_CURR_TAB]; /*should I also say extern?*/

double current_from_time(double t);
double interp_lin_table(double xpoint, double *x, double *y, const int N);

#endif
