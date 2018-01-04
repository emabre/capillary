#include "pluto.h"
#include "gamma_transp.h"
#include "prototypes.h"

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
    double mu_aux; /*Auxiliary variable*/
    Data d_corrected_r, d_corrected_z;
    double ***vr_c_r, ***vr_c_z;
    double ***vz_c_r, ***vz_c_z;
    double ***rho_c_r, ***rho_c_z;
    #ifdef WRITE_T_MU_NE_IONIZ
      double ***T_c_r, ***T_c_z;
    #endif
  #endif

/******************************************************/
/*I allocate space for all the variables of the output*/
/******************************************************/
  #if MULTIPLE_GHOSTS==YES
    // I export the corrected vr, vz, rho
    vr_c_r = GetUserVar("vr_c_r");
    vr_c_z = GetUserVar("vr_c_z");
    vz_c_r = GetUserVar("vz_c_r");
    vz_c_z = GetUserVar("vz_c_z");
    rho_c_z = GetUserVar("rho_c_z");
    rho_c_r = GetUserVar("rho_c_r");
    #ifdef WRITE_T_MU_NE_IONIZ
      T_c_z= GetUserVar("T_c_z");
      T_c_r= GetUserVar("T_c_r");
    #endif
  #endif
  // export Internal boundary flags
  interBound = GetUserVar("interBound");
  #ifdef WRITE_T_MU_NE_IONIZ
    mu = GetUserVar("mu");
    T = GetUserVar("T");
    #if EOS==PVTE_LAW
      ioniz = GetUserVar("ioniz");
      ne = GetUserVar("ne");
    #endif
  #endif

/******************************************************/
/* Now I compute the variables to export*/
/******************************************************/
  DOM_LOOP(k,j,i){
    if ((int) (d->flag[k][j][i] & FLAG_INTERNAL_BOUNDARY)) {
      interBound[k][j][i] = 33.3; //Just a conventional number
    } else {
      interBound[k][j][i] = 0.0; //Just a conventional number
    }
  }

  // export T and mu ne and ionization degree
  #ifdef WRITE_T_MU_NE_IONIZ
    DOM_LOOP(k,j,i){
      #if EOS==IDEAL
        mu[k][j][i] = MeanMolecularWeight(d->Vc);
        T[k][j][i] = d->Vc[PRS][k][j][i]/d->Vc[RHO][k][j][i]*KELVIN*mu[k][j][i];
      #elif EOS==PVTE_LAW
        for (nv=NVAR; nv--;) v[nv] = d->Vc[nv][k][j][i];
        if (GetPV_Temperature(v, &(T[k][j][i]) )!=0) {
          print1("ComputeUserVar:[Ema] Error computing temperature!");
        }
        // print1("\nI just assigned %g to T[%d][%d][%d] for output",T[k][j][i], k,j,i);
        GetMu(T[k][j][i], v[RHO], &(mu[k][j][i]));
        ioniz[k][j][i] = 1/mu[k][j][i] - 1;
        ne[k][j][i] = ioniz[k][j][i] * v[RHO] / CONST_mp;
      #endif
    }
  #endif

  #if MULTIPLE_GHOSTS==YES
    d_corrected_r = *d;
    d_corrected_z = *d;
    ApplyMultipleGhosts(&d_corrected_r, 0);
    ApplyMultipleGhosts(&d_corrected_z, 1);

    DOM_LOOP(k,j,i){
      vr_c_r[k][j][i] = d_corrected_r.Vc[iVR][k][j][i];
      vr_c_z[k][j][i] = d_corrected_z.Vc[iVR][k][j][i];
      vz_c_r[k][j][i] = d_corrected_r.Vc[iVZ][k][j][i];
      vz_c_z[k][j][i] = d_corrected_z.Vc[iVZ][k][j][i];
      rho_c_r[k][j][i] = d_corrected_r.Vc[RHO][k][j][i];
      rho_c_z[k][j][i] = d_corrected_z.Vc[RHO][k][j][i];
    }

    // I take the corrected T
    #ifdef WRITE_T_MU_NE_IONIZ
      DOM_LOOP(k,j,i){
        #if EOS==IDEAL
          mu_aux = MeanMolecularWeight(d_corrected_r.Vc);
          T_c_r[k][j][i] = d_corrected_r.Vc[PRS][k][j][i]/d_corrected_r.Vc[RHO][k][j][i]*KELVIN*mu_aux;
        #elif EOS==PVTE_LAW
          for (nv=NVAR; nv--;) v[nv] = d_corrected_r.Vc[nv][k][j][i];
          if (GetPV_Temperature(v, &(T_c_r[k][j][i]) )!=0) {
            print1("ComputeUserVar:[Ema] Error computing temperature!");
          }
          GetMu(T_c_r[k][j][i], v[RHO], &(mu_aux));
        #endif
      }

      DOM_LOOP(k,j,i){
        #if EOS==IDEAL
          mu_aux = MeanMolecularWeight(d_corrected_z.Vc);
          T_c_z[k][j][i] = d_corrected_z.Vc[PRS][k][j][i]/d_corrected_z.Vc[RHO][k][j][i]*KELVIN*mu_aux;
        #elif EOS==PVTE_LAW
          for (nv=NVAR; nv--;) v[nv] = d_corrected_z.Vc[nv][k][j][i];
          if (GetPV_Temperature(v, &(T_c_z[k][j][i]) )!=0) {
            print1("ComputeUserVar:[Ema] Error computing temperature!");
          }
          GetMu(T_c_z[k][j][i], v[RHO], &(mu_aux));
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
