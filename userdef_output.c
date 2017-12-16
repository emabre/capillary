#include "pluto.h"

#define WRITE_T_MU_NE_IONIZ

/* *************************************************************** */
void ComputeUserVar (const Data *d, Grid *grid)
/*
 *
 *  PURPOSE
 *
 *    Define user-defined output variables
 *
 * EMA: la modifico per stamparmi la temperatura
 *
 ***************************************************************** */
{
  int i, j, k;
  double ***interBound;
  #ifdef WRITE_T_MU_NE_IONIZ
    double ***T, ***ioniz, ***mu, ***ne;
    double v[NVAR]; /*[Ema] I hope that NVAR as dimension is fine!*/
    int nv;
  #endif

  // export Internal boundary flags
  interBound = GetUserVar("interBound");
  DOM_LOOP(k,j,i){
    if ((int) (d->flag[k][j][i] & FLAG_INTERNAL_BOUNDARY)) {
      interBound[k][j][i] = 33.3; //Just a conventional number
    } else {
      interBound[k][j][i] = 0.0; //Just a conventional number
    }
  }

  // export T and mu ne and ionization degree
  #ifdef WRITE_T_MU_NE_IONIZ
    mu = GetUserVar("mu");
    T = GetUserVar("T");
    #if EOS==PVTE_LAW
      ioniz = GetUserVar("ioniz");
      ne = GetUserVar("ne");
    #endif

    DOM_LOOP(k,j,i){
      #if EOS==IDEAL
        mu[k][j][i] = MeanMolecularWeight(d->Vc);
        T[k][j][i] = d->Vc[PRS][k][j][i]/d->Vc[RHO][k][j][i]*KELVIN*mu[k][j][i];
      #elif EOS==PVTE_LAW
        for (nv=NVAR; nv--;) v[nv] = d->Vc[nv][k][j][i];
        if (GetPV_Temperature(v, &(T[k][j][i]) )!=0) {
          print1("ComputeUserVar:[Ema] Error computing temperature!");
        }
        GetMu(T[k][j][i], v[RHO], &(mu[k][j][i]));
        ioniz[k][j][i] = 1/mu[k][j][i] - 1;
        ne[k][j][i] = ioniz[k][j][i] * v[RHO] / CONST_mp;
      #endif
    }
  #endif


}
/* ************************************************************* */
void ChangeDumpVar ()
/*
 *
 *
 *************************************************************** */
{
  Image *image;
  // SetDumpVar("bphi", VTK_OUTPUT, YES);
  // SetDumpVar("r_coordI", VTK_OUTPUT, YES);
  // SetDumpVar("r_coordJ", VTK_OUTPUT, YES);
  // SetDumpVar("r_coordK", VTK_OUTPUT, YES);
  // SetDumpVar("rho_ema", VTK_OUTPUT, YES);
  // SetDumpVar("bx3", VTK_OUTPUT, YES);


}
