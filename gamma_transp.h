#ifndef GAMMA_TRANSP_H
#define GAMMA_TRANSP_H


#define Q_ELEM 4.8032e-10 /*Elementary (electron) charge in statcoulomb*/
/*Epsilon 0 constant in Esaulov's paper Plasma Phys control Fusion 43 571 */

#define EPS_0_ESAU 5.81e-11 /* in erg*/

#define ELEC_DENS(rho,z)  ((z) * (rho) / (CONST_mp+CONST_me))

#define IONIZMIN 1.0e-9

const double extern unit_Mfield;

const double extern Gamma_jil[2][10][4];

double Gamma_i(int i, double xe, double w);
double Gamma_ij(int i, int j, double w);
double cl_ei(double ne, double kT);
double cl_en(double ioniz, double kT);
double cl_ee(double ne, double kT);
double xe(double B, double f_coll_e);
double w(double cl_ee, double cl_ei, double cl_en);
double freq_coll_e(double z, double rho, double kT);
double elRes_norm(double z, double rho, double kT, int corr, double B);
double thermCond_norm(double z,double rho, double kT, int corr, double B);

double cl_ei_el_DUED(double ne, double kT, double z);
double cl_ei_th(double ne, double kT, double z);
double erg2keV(double kT);
double elRes_norm_DUED(double z, double rho, double kT);
double thermCond_norm_DUED(double z, double rho, double kT);


#endif
