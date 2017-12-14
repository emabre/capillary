#ifndef GAMMA_TRANSP_H
#define GAMMA_TRANSP_H

#define Q_ELEM 4.8032e-10 /*Elementary (electron) charge in statcoulomb*/
/*Epsilon 0 constant in Esaulov's paper Plasma Phys control Fusion 43 571 */
#define EPS_0_ESAU 5.81e-11 /* in erg*/

double Gamma_jil[2][10][4];
int extern idx_rcap, idx_zcap;
double extern rcap,zcap;

double Gamma_i(int i, double xe, double w);
double Gamma_ij(int i, int j, double w);
double cl_ei(double ne, double T);
double cl_en(double ioniz, double T);
// double cl_ee(double ne, double T);

int find_idx_closest(double *vec, int Nvec, double v);

#endif
