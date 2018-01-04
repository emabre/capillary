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
  #if MULTIPLE_GHOSTS==YES
    double ***mu_aux; /*Auxiliary variable*/
    Data *d_corrected_r, *d_corrected_z;
    double ***vx1_corr_r, ***vx1_corr_z;
    double ***vx2_corr_r, ***vx2_corr_z;
    double ***rho_corr_r, ***rho_corr_z;
    #ifdef WRITE_T_MU_NE_IONIZ
      double ***T_corr_r, ***T_corr_z;
    #endif
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
      #elif EOS==PVTE_LE_LAW
        mu = GetUserVar("mu");
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

  #if MULTIPLE_GHOSTS==YES
    *d_corrected_r = *d;
    *d_corrected_z = *d;
    ApplyMultipleGhosts(d_corrected_r, 0);
    ApplyMultipleGhosts(d_corrected_z, 1);

    // I take the corrected vx1, vx2, rho
    vx1_corr_r = GetUserVar("vx1_corr_r");
    vx1_corr_z = GetUserVar("vx1_corr_z");
    vx2_corr_r = GetUserVar("vx2_corr_r");
    vx2_corr_z = GetUserVar("vx2_corr_z");
    rho_corr_z = GetUserVar("rho_corr_z");
    rho_corr_r = GetUserVar("rho_corr_r");
    DOM_LOOP(k,j,i){
      vx1_corr_r[k][j][i] = d_corrected_r->Vc[VX1][k][j][i];
      vx1_corr_z[k][j][i] = d_corrected_z->Vc[VX1][k][j][i];
      vx2_corr_r[k][j][i] = d_corrected_r->Vc[VX2][k][j][i];
      vx2_corr_z[k][j][i] = d_corrected_z->Vc[VX2][k][j][i];
      rho_corr_r[k][j][i] = d_corrected_r->Vc[RHO][k][j][i];
      rho_corr_z[k][j][i] = d_corrected_z->Vc[RHO][k][j][i];
    }

    // I take the corrected T
    #ifdef WRITE_T_MU_NE_IONIZ
      T_corr_r= GetUserVar("T_corr_r");
      DOM_LOOP(k,j,i){
        #if EOS==IDEAL
          mu_aux[k][j][i] = MeanMolecularWeight(d_corrected_r->Vc);
          T_corr_r[k][j][i] = d_corrected_r->Vc[PRS][k][j][i]/d_corrected_r->Vc[RHO][k][j][i]*KELVIN*mu_aux[k][j][i];
        #elif EOS==PVTE_LAW
          for (nv=NVAR; nv--;) v[nv] = d_corrected_r->Vc[nv][k][j][i];
          if (GetPV_Temperature(v, &(T_corr_r[k][j][i]) )!=0) {
            print1("ComputeUserVar:[Ema] Error computing temperature!");
          }
          GetMu(T_corr_r[k][j][i], v[RHO], &(mu_aux[k][j][i]));
        #endif
      }

      T_corr_z= GetUserVar("T_corr_z");
      DOM_LOOP(k,j,i){
        #if EOS==IDEAL
          mu_aux[k][j][i] = MeanMolecularWeight(d_corrected_z->Vc);
          T_corr_z[k][j][i] = d_corrected_z->Vc[PRS][k][j][i]/d_corrected_z->Vc[RHO][k][j][i]*KELVIN*mu_aux[k][j][i];
        #elif EOS==PVTE_LAW
          for (nv=NVAR; nv--;) v[nv] = d_corrected_z->Vc[nv][k][j][i];
          if (GetPV_Temperature(v, &(T_corr_z[k][j][i]) )!=0) {
            print1("ComputeUserVar:[Ema] Error computing temperature!");
          }
          GetMu(T_corr_z[k][j][i], v[RHO], &(mu_aux[k][j][i]));
        #endif
      }
    #endif
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
