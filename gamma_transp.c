#include <math.h>
#include "pluto.h"
#include "gamma_transp.h"

/*
This file contains:
- the definition of the gamma parameters
as in:
[1] N. A. Bobrova and P. V. Sasorov, Plasma Phys. Rep. 19, 409 1993
This file contains also some other quantities and parameter useful for computing quantities like thermal conductivity, electrical resistivity.
- other plasma parameters, computed as in Ref[1] and in:
[2] Plasma Phys. Control. Fusion 43 (2001) 571–588
[3] Physical Review E, Volume 65, 016407
- plasma parameters as defined inside dued
*/

#define ERG2KEV(kT) ((kT)/CONST_eV/1e3)

const double sigma_ea = 1.4e-15;  // See Golant book


/*GAMMA_ilj[i,l,j] is the parameter "capital gamma"
with subscripts i+1,j+1 and superscript l*/
const double Gamma_jil[2][10][4] = {
                     {{13/4,            2,          0,             0},
                      {5/2,             0,          0,             0},
                      {212.7/56,        103.2/56,   0,             0},
                      {3/2,             0,          0,             0},
                      {3630.33/(49*16), 123.705/49, 0,             0},
                      {477/280,         0,          0,             0},
                      {5866.01/(49*16), 412.69/49,  132.52/49,     0},
                      {1.2,             1.2,        0,             0},
                      {1,               0,          0,             0},
                      {265.32/49,       440.94/49,  793.21/(49*4), 0}},

                     {{939.61/(49*16),  375.74/49,  427.68/49,     129.6/49},
                      {32079.7/(49*64), 632.025/49, 227.7/49,      0},
                      {7.161/49,        59.394/49,  127.728/49,    51.84/49},
                      {42.9675/49,      100.665/49, 70.92/49,      0},
                      {3.3201/49,       34.1064/49, 95.7888/49,    41.472/49},
                      {4.608/49,        21.744/49,  36.432/49,     0},
                      {0.31,            12.08/7,    5.76/7,        0},
                      {58.752/49,       243.253/49, 294.051/49,    109.47/49},
                      {149.4/49,        254.46/49,  465.61/(49*4), 0},
                      {576/700,         1806/700,   1068/700,      0}}
                  };

/*---------------------------------------------*/
/*Returns the value of the parameter "capital gamma (xe,w)" with subscript i.
  (i=1,2,3,4,5,6,8,9)*/
/*---------------------------------------------*/
double Gamma_i(int i, double xe, double w){
    int I;
    if (i>=1 && i<=6){
        I=7;
    }else if (i==8 || i==9){
        I=10;
    }else{
        print1("[Ema] Gamma_i is defined only for i=1,2,3,4,5,6,8,9, no 7 or 10 or others");
    }
    return (Gamma_ij(i,1,w)*pow(xe,2) + Gamma_ij(i,2,w)) / (pow(xe,4) + Gamma_ij(I,1,w)*pow(xe,2) + pow((Gamma_ij(I,2,w)),2));
}

/*---------------------------------------------*/
/*Returns the value of the parameter "capital gamma (w)" with subscripts i,j
  (i=1,2,3,4,5,6,7,8,9,10; j=1,2).
  I diminuish by 1 both i and j, so that I can write the formula
  for computing gamma_ij in a way which is more similar to the one written in Ref [1]*/
/*---------------------------------------------*/
double Gamma_ij(int i, int j, double w){
    i-=1;
    j-=1;
    return (Gamma_jil[j][i][3]*pow(w,3) + Gamma_jil[j][i][2]*pow(w,2) + Gamma_jil[j][i][1]*w + Gamma_jil[j][i][0]);
}

/*---------------------------------------------*/
/* Coulomb log for e-i collisions as in Bobrova/Esaulov,
   neglecting factor 1/(y1^2 +Te/eps_0),
   assuming Ti = Te and assuming that we have only hydrogen*/
/*---------------------------------------------*/
double cl_ei(double ne, double kT) {
  double const q_elem6 = Q_ELEM*Q_ELEM*Q_ELEM*Q_ELEM*Q_ELEM*Q_ELEM;
  double cl;

  cl = 0.5*log(9/(4*CONST_PI) * 1/(ne*q_elem6) * (kT*kT*kT)/2);
  if (cl<1) {
    cl = 1;
  }
  return cl;
}

/*---------------------------------------------*/
/* Coulomb log for e-n collisions as in Bobrova/Esaulov,
   simplifying factor 1/(1 + (eps_0/Te)^2) to Te^2/eps_0^2,
   assuming Ti = Te and assuming that we have only hydrogen */
/*---------------------------------------------*/
double cl_en(double z, double kT){
  return 0.25*(1-z)/z * (kT*kT/(EPS_0_ESAU*EPS_0_ESAU));
}

/*---------------------------------------------
  Coulomb log for e-e collisions as in Bobrova/Esaulov,\
  neglecting factor 1/(1 +Te/eps_0),\
  assuming Ti = Te and assuming that we have only hydrogen
---------------------------------------------*/
double cl_ee(double ne, double kT) {
  double const q_elem6 = Q_ELEM*Q_ELEM*Q_ELEM*Q_ELEM*Q_ELEM*Q_ELEM;
  double cl = 0.5*log(9/(16*CONST_PI) * (kT*kT*kT)/(ne*q_elem6));
  if (cl<1) {
      cl = 1;
    }
  return cl;
}

/*------------------------------------------------
  Ratio between electron gyro frequency (pulsazione, 2*pi*frequenza)
  and electron collision freqency
---------------------------------------------------*/
double xe(double B, double f_coll_e) {
  double omega_Be = Q_ELEM*B/(CONST_me*CONST_c);
  return omega_Be/f_coll_e;
}

/*------------------------------------------------
  Parameter w
------------------------------------------------*/
double w(double cl_ee, double cl_ei, double cl_en){
  double ycl = cl_ei+cl_en;
  return cl_ee/(sqrt(2)*ycl);
}

/*-------------------------------------------------
  Electron collision frequency according to Esaulov/Bobrova
-------------------------------------------------*/
double freq_coll_e(double z, double rho, double kT) {
  double ne = ELEC_DENS(rho,z);
  double ycl, ni_e;
  double const q_elem4 = Q_ELEM*Q_ELEM*Q_ELEM*Q_ELEM;

  ycl = cl_ei(ne,kT)+cl_en(z,kT);
  ni_e = (4*sqrt(2*CONST_PI)/3) * (q_elem4)/sqrt(CONST_me) * ne/(pow(kT,1.5))*ycl;
  return ni_e;
}

/*------------------------------------------------
  Electrical resistivity according to Esaulov/Bobrova.
  If corr = False no correction using Gamma parameter is done (as if Gamma were 1)
------------------------------------------------*/
double elRes_norm(double z, double rho, double kT, int corr, double B) {
  double ne = ELEC_DENS(rho,z);
  double corr_factor, gamma5Bob;

  if (!corr) {
      gamma5Bob = 0;
  } else {
      gamma5Bob = Gamma_i(5, xe(B, freq_coll_e(z,rho,kT)), w(cl_ee(ne, kT), cl_ei(ne, kT), cl_en(z, kT)));
  }
  corr_factor = 1-gamma5Bob;
  return CONST_me * freq_coll_e(z,rho,kT) / (Q_ELEM*Q_ELEM * ne) * corr_factor;
}

/*------------------------------------------------
  Electrical resistivity as in code DUED.
  (the expression implemented here is valid ONLY FOR z<1).
  Resistivity of Spitzer, in direction parallel to magnetic field,
  (ONLY VALID FOR HYDROGEN), corrected incuding collisions e-e.
------------------------------------------------*/
double elRes_norm_DUED(double z, double rho, double kT) {
  double TkeV = ERG2KEV(kT);
  double sqrt_TkeV = TkeV*TkeV;
  double Z_ion = fmax(1.0,z);
  double Z_ion_sq = Z_ion*Z_ion;
  double ne = ELEC_DENS(rho,z);
  double elRes_norm_en, elRes_norm_ei;
  double cl_ei_el, cl_quant;
  // double Zmin = fmin(z,1.0);

  /*-----------------------------------------------------------------------*/
  /*Electrical resitivity due to collisions with neutrals according to DUED*/
  
  sqrt_TkeV = sqrt(TkeV);
  elRes_norm_en = 1.03e-14 * (1-z)/(z+1e-20) * sqrt_TkeV * (sigma_ea/1.4e-15);
  /*-----------------------------------------------------------------------*/

  /*-----------------------------------------------------------------------*/
  /*Electrical resitivity due to collisions with ions according to DUED.*/

  // Coulomb log computation
  cl_quant = 7.1 - 0.5*log(ne/1e21) + log(TkeV);
  if (TkeV < 0.01/Z_ion_sq){
      cl_ei_el = cl_quant + 2.3 + 0.5*log(TkeV/Z_ion_sq);
  } else {
      cl_ei_el = cl_quant;
  }
  cl_ei_el = fmax(cl_ei_el,1.0);

  /*Old (equivalent) version, less performing,
    (I used to call a function) */
  // cl_ei_el = cl_ei_el_DUED(ne, kT, z);

   elRes_norm_ei = 1.840e-19 * cl_ei_el /(TkeV*sqrt_TkeV) * Z_ion;
  // Old (equivalent) version, less performing
  // elRes_norm_ei = 1.840e-19 * Z_ion * cl_ei_el /(TkeV*sqrt_TkeV) * Z_ion*Zmin/(z+1e-20);  
  /*-----------------------------------------------------------------------*/

  return elRes_norm_ei + elRes_norm_en;
}
/*---------------------------------------------------
  Coulomb log for e-i according to DUED, for computation of electrical resitivity
---------------------------------------------------*/
double cl_ei_el_DUED(double ne, double kT, double z){
  double TkeV, Z_ion, cl_quant, cl;

  TkeV = ERG2KEV(kT);
  Z_ion = fmax(1.0,z);
  cl_quant = 7.1 - 0.5*log(ne/1e21) + log(TkeV);
  if (TkeV < 0.01/(Z_ion*Z_ion)){
      cl = cl_quant + 2.3 + 0.5*log(TkeV/(Z_ion*Z_ion));
  } else {
      cl = cl_quant;
  }
  cl = fmax(cl,1.0);
  return cl;
}

/*------------------------------------------------
  Thermal conductivity according to Esaulov/Bobrova CONVERTED TO work for an equation of conduction where
  temperature gradients expressed WITH TEMPERTURE DEFINED IN KELVIN (Kelvin/cm). (so I
  multiply the formula with CONST_kB)
  If corr = False no correction using Gamma parameter is done (as if Gamma were 1)
------------------------------------------------*/
double thermCond_norm(double z,double rho, double kT, int corr, double B) {
  double ne = ELEC_DENS(rho,z);
  double gamma1Bob;

  if (!corr) {
      // print("[thermCond_norm]For now I don't consider the gamma correction");
      gamma1Bob = 1;
  } else {
      gamma1Bob = Gamma_i(1, xe(B, freq_coll_e(z,rho,kT)), w(cl_ee(ne, kT), cl_ei(ne, kT), cl_en(z, kT)));
  }
  /*Notice that I */
  return CONST_kB*ne*kT/(CONST_me*freq_coll_e(z,rho,kT))*gamma1Bob;
}

/*------------------------------------------------
  Thermal conductivity of Hydrogen according to DUED. Works for an equation of conduction where
  temperature gradients expressed WITH TEMPERTURE DEFINED IN KELVIN (Kelvin/cm).
  I removed the "chie multiplier" which is (to my knowledge) always 1
  for the cases of our interest
------------------------------------------------*/
double thermCond_norm_DUED(double z, double rho, double kT) {
  double ne = ELEC_DENS(rho,z);
  double cl, deleps;
  double Z_ion = fmax(1.0,z); // This is called ZEI in dued
  double Z_ion_sq = Z_ion*Z_ion;
  double Zion_pow_089 = pow(Z_ion,0.89);
  double T=kT/CONST_kB;
  double Tsq = T*T;
  double chie;

  // I compute the Coulomb-Log (different from the one computed for resistivity)
  // Old version (just for performance): I used a function
  // cl = cl_ei_th(ne, kT, z);
  if (T <= 1.57e5*Z_ion_sq){
      cl = -18.34 + log( T*sqrt(T) / (sqrt(ne/CONST_NA) * Z_ion));
  } else {
      cl = -12.36 + log(T / sqrt(ne/CONST_NA));
  }

  deleps = 0.4*Zion_pow_089/(3.25+Zion_pow_089);

  chie = deleps * 1.96e-4 * (Tsq*sqrt(T)) / cl / Z_ion;
  // Old (equivalent) version, less performing
  // chie = deleps * (z/fmin(1.0,z)/(Z_ion_sq)) * 1.96e-4 * (Tsq*sqrt(T)) / cl;

  if (z<1) {
      chie = chie / (1 + 4e-9*(1-z)/z/Z_ion_sq *(sigma_ea/1.4e-15) * Tsq );
      // Old (equivalent) version, less performing
      // chie = chie / (1 + 4e-9*fmax(0.0,1-z)/fmin(1.0,z)/Z_ion_sq *(sigma_ea/1.4e-15) * Tsq );
  }
  return chie;
}

/*-----------------------------------------------------------------
Coulomb log for e-i according to DUED, for computation of thermal conductivity'''
-----------------------------------------------------------------*/
double cl_ei_th(double ne, double kT, double z){
  double T, Z_ion, cl;
  T = kT / CONST_kB;
  Z_ion = fmax(1,z);
  if (T <= 1.57e5*Z_ion*Z_ion){
      cl = -18.34 + log( T*sqrt(T) / (sqrt(ne/CONST_NA) * Z_ion));
  } else {
      cl = -12.36 + log(T / sqrt(ne/CONST_NA));
  }
  return cl;
}