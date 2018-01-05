#include <math.h>
#include "pluto.h"
#include "gamma_transp.h"

//debug macro
// #define DBG_FIND_CLOSEST

/*This file contains the definition of the gamma parameters
as in:
[1] N. A. Bobrova and P. V. Sasorov, Plasma Phys. Rep. 19, 409 1993
This file contains also some other quantities and parameter useful for computing
quantities like thermal conductivity, electrical resistivity.*/
/*GAMMA_ilj[i,l,j] is the parameter "capital gamma"
with subscripts i+1,j+1 and superscript l*/
double Gamma_jil[2][10][4] = {
                    {{13/4, 2, 0, 0},
                    {5/2, 0, 0, 0},
                    {212.7/56, 103.2/56, 0, 0},
                    {3/2, 0, 0, 0},
                    {3630.33/(49*16), 123.705/49, 0, 0},
                    {477/280, 0, 0, 0},
                    {5866.01/(49*16), 412.69/49, 132.52/49, 0},
                    {1.2, 1.2, 0, 0},
                    {1, 0, 0, 0},
                    {265.32/49, 440.94/49, 793.21/(49*4), 0}},
                    {{939.61/(49*16), 375.74/49, 427.68/49, 129.6/49},
                    {32079.7/(49*64), 632.025/49, 227.7/49, 0},
                    {7.161/49, 59.394/49, 127.728/49, 51.84/49},
                    {42.9675/49, 100.665/49, 70.92/49, 0},
                    {3.3201/49, 34.1064/49, 95.7888/49, 41.472/49},
                    {4.608/49, 21.744/49, 36.432/49, 0},
                    {0.31, 12.08/7, 5.76/7, 0},
                    {58.752/49, 243.253/49, 294.051/49, 109.47/49},
                    {149.4/49, 254.46/49, 465.61/(49*4), 0},
                    {576/700, 1806/700, 1068/700, 0}}
                  };

int idx_rcap = 0;
int idx_zcap = 0;
int idx_start_electr = 0;

Corr d_correction[3] = { {},{},{} };

double Gamma_i(int i, double xe, double w){
    //Returns the value of the parameter "capital gamma (xe,w)" with subscript i. See Ref [1] for info
    int J;

    if (i>=1 && i<=6){
        J=7;
    }else if (i==8 || i==9){
        J=10;
    }else{
        print1("[Ema] Gamma_i is defined only for i=1,2,3,4,5,6,8,9, no 7 or 10 or others");
    }
    return (Gamma_ij(i,1,w)*pow(xe,2) + Gamma_ij(i,2,w)) / (pow(xe,4) + Gamma_ij(J,1,w)*pow(xe,2) + pow((Gamma_ij(J,2,w)),2));
}

/*---------------------------------------------*/
/* Coulomb log for e-i collisions as in Bobrova/Esaulov,
   neglecting factor 1/(y1^2 +Te/eps_0),
   assuming Ti = Te and assuming that we have only hydrogen*/
/*---------------------------------------------*/
double cl_ei(double ne, double T) {
  double const q_elem6 = Q_ELEM*Q_ELEM*Q_ELEM*Q_ELEM*Q_ELEM*Q_ELEM;
  double cl;

  cl = 0.5*log(9/(4*CONST_PI) * 1/(ne*q_elem6)) * T*T*T/2;
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
double cl_en(double ioniz, double T){
  return 0.25*(1-ioniz)/ioniz * T*T/EPS_0_ESAU;
}

// double cl_ee(double ne, double T){
//   double const q_elem6 = Q_ELEM*Q_ELEM*Q_ELEM*Q_ELEM*Q_ELEM*Q_ELEM;
//   double cl;
//
//   cl = 0.5*log(9/(16*CONST_PI) * T*T*T/(ne*q_elem6));
//   if (cl<1) {
//     cl = 1;
//   }
//
//   return cl;
// }

//Returns the value of the parameter capital gamma (w) with subscripts i,j
/*I diminuish by 1 both i and j, so that I can write the formula
for computing gamma_ij in a way which is more similar to the one written in Ref [1]*/
double Gamma_ij(int i, int j, double w){
    i-=1;
    j-=1;
    return (Gamma_jil[j][i][3]*pow(w,3) + Gamma_jil[j][i][2]*pow(w,2) + Gamma_jil[j][i][1]*w + Gamma_jil[j][i][0]);
}

/*Finds the index of the element in vec closest to the value of v*/
int find_idx_closest(double *vec, int Nvec, double v){
  int i, i_mindiff;
  double diff;

  #ifdef DBG_FIND_CLOSEST
    print1("\nvec[0]:%g, v:%g", vec[0], v);
    print1("\nfabs(vec[0]-v):%g", fabs(vec[0]-v));
  #endif

  diff = fabs(vec[0]-v);
  for (i=1;i<Nvec;i++){
    if (fabs(vec[i]-v) < diff) {
      #ifdef DBG_FIND_CLOSEST
        print1("\nabs(vec[%d](=%e)-%e)<%e",i,vec[i],v,diff);
      #endif
      diff = fabs(vec[i]-v);
      i_mindiff = i;
    } else {
      #ifdef DBG_FIND_CLOSEST
        print1("\nfabs(vec[%d](=%e)-%e)>%e",i,vec[i],v,diff);
      #endif
    }
  }
  return i_mindiff;
}
