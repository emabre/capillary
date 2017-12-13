#include "pluto.h"
#include "gamma_transp.h"
#include "additional_const_disch.h"
// #include <math.h>

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

double cl_ei(double ne, double T) {
  /*---------------------------------------------*/
  /* Coulomb log for e-i collisions as in Bobrova/Esaulov,
     neglecting factor 1/(y1^2 +Te/eps_0),
     assuming Ti = Te and assuming that we have only hydrogen*/
  /*---------------------------------------------*/
  double const q_elem6 = Q_ELEM*Q_ELEM*Q_ELEM*Q_ELEM*Q_ELEM*Q_ELEM;
  double cl;

  cl = 0.5*log(9/(4*CONST_PI) * 1/(ne*q_elem6)) * T*T*T/2;
  if (cl<1) {
    cl = 1;
  }
  return cl;
}

double cl_en(double ioniz, double T){
  /*---------------------------------------------*/
  /* Coulomb log for e-n collisions as in Bobrova/Esaulov,
     simplifying factor 1/(1 + (eps_0/Te)^2) to Te^2/eps_0^2,
     assuming Ti = Te and assuming that we have only hydrogen */
  /*---------------------------------------------*/
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

double Gamma_ij(int i, int j, double w){
    //Returns the value of the parameter capital gamma (w) with subscripts i,j
    /*I diminuish by 1 both i and j, so that I can write the formula
    for computing gamma_ij in a way which is more similar to the one written in Ref [1]*/
    i-=1;
    j-=1;
    return (Gamma_jil[j][i][3]*pow(w,3) + Gamma_jil[j][i][2]*pow(w,2) + Gamma_jil[j][i][1]*w + Gamma_jil[j][i][0]);
}
