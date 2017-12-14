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
  double ***press;
  int ***interBound;
  double v[NVAR]; /*[Ema] I hope that NVAR as dimension is fine!*/
  double ***phi, ***kpar;

  press = GetUserVar("press");
  interBound = GetUserVar("interBound");

  DOM_LOOP(k,j,i){
    press[k][j][i] = d->Vc[PRS][k][j][i];
    interBound[k][j][i] = d->flag[k][j][i];
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
