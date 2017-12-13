#include "pluto.h"

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
  int i, j, k, nv;
  double ***T, ***ioniz, ***mu, ***ne, ***knor;
  double v[NVAR]; /*[Ema] I hope that NVAR as dimension is fine!*/
  double ***phi, ***kpar;

  T = GetUserVar("T");
  mu = GetUserVar("mu");
  // knor = GetUserVar("knor");

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
    for (nv=NVAR; nv--;) v[nv] = d->Vc[nv][k][j][i];
    // I am not sure I used correctly used the grid to get x1, x2 and x3.
    // TC_kappa( v, *grid[IDIR].x, *grid[JDIR].x, *grid[KDIR].x,
    //           &(kpar[k][j][i]), &(knor[k][j][i]), &(phi[k][j][i]) );
  }
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
