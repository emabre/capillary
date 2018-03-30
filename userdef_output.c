#include "pluto.h"
#include "gamma_transp.h"
#include "capillary_wall.h"
#include "current_table.h"
#include "prototypes.h"

#define WRITE_T_MU_NE_IONIZ YES
#define WRITE_J1D YES
#define WRITE_J YES

#if WRITE_J1D == YES
  void ComputeJ1DforOutput(const Data *d, Grid *grid, double ***Jz1D);
#endif

#if WRITE_J == YES
  void ComputeJforOutput(const Data *d, Grid *grid, double ***Ji, double ***Jj);
#endif

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

  #if WRITE_J1D == YES
    double ***Jz1D;
  #endif
  #if WRITE_J == YES
    double ***Jr;
    double ***Jz;
  #endif

  #if WRITE_T_MU_NE_IONIZ==YES
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
    #if WRITE_T_MU_NE_IONIZ == YES
      T_c_z= GetUserVar("T_c_z");
      T_c_r= GetUserVar("T_c_r");
    #endif
  #endif
  // export Internal boundary flags
  interBound = GetUserVar("interBound");
  #if WRITE_T_MU_NE_IONIZ == YES
    mu = GetUserVar("mu");
    T = GetUserVar("T");
    #if EOS==PVTE_LAW
      ioniz = GetUserVar("ioniz");
      ne = GetUserVar("ne");
    #endif
  #endif
  #if WRITE_J1D == YES
    Jz1D = GetUserVar("Jz1D");
  #endif
  #if WRITE_J == YES
    Jr = GetUserVar("Jr");
    Jz = GetUserVar("Jz");
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
  #if WRITE_T_MU_NE_IONIZ==YES
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
        ne[k][j][i] = ioniz[k][j][i] * (v[RHO]*UNIT_DENSITY) / CONST_mp;
      #endif
    }
  #endif

  #if MULTIPLE_GHOSTS==YES
    // memcpy(&d_corrected_r, d, sizeof(*d));
    // memcpy(&d_corrected_z, d, sizeof(*d));
    // d_corrected_r = *d;
    // d_corrected_z = *d;
    alloc_Data(&d_corrected_r);
    alloc_Data(&d_corrected_z);
    print1("\nRicordati di deallocare d_corrected_r e ..._z !");
    copy_Data_Vc(&d_corrected_z, d);
    copy_Data_Vc(&d_corrected_r, d);
    ApplyMultipleGhosts(&d_corrected_r, 0);
    ApplyMultipleGhosts(&d_corrected_z, 1);

    DOM_LOOP(k,j,i){
      vr_c_r[k][j][i] = d_corrected_r.Vc[iVR][k][j][i]*UNIT_VELOCITY;
      vr_c_z[k][j][i] = d_corrected_z.Vc[iVR][k][j][i]*UNIT_VELOCITY;
      vz_c_r[k][j][i] = d_corrected_r.Vc[iVZ][k][j][i]*UNIT_VELOCITY;
      vz_c_z[k][j][i] = d_corrected_z.Vc[iVZ][k][j][i]*UNIT_VELOCITY;
      rho_c_r[k][j][i] = d_corrected_r.Vc[RHO][k][j][i]*UNIT_DENSITY;
      rho_c_z[k][j][i] = d_corrected_z.Vc[RHO][k][j][i]*UNIT_DENSITY;
    }

    // I take the corrected T
    #if WRITE_T_MU_NE_IONIZ==YES
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
//aszascs
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
  #if WRITE_J1D == YES
    ComputeJ1DforOutput(d, grid, Jz1D);
  #endif
  #if WRITE_J == YES
    ComputeJforOutput(d, grid, Jr, Jz);
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

#if WRITE_J1D == YES
void ComputeJ1DforOutput(const Data *d, Grid *grid, double ***Jz1D){
  int i, j, k;
  // double Jstagg[NX3_TOT-1][NX2_TOT-1][NX1_TOT-1];
  double *r, ***B;
  double Jleft, Jright;
  // RBox box;
  double unit_Mfield;

  unit_Mfield = COMPUTE_UNIT_MFIELD(UNIT_VELOCITY, UNIT_DENSITY);

  r = grid[IDIR].x;
  // x = grid[0].x_glob;
  // print1("for 0, gbeg: %d, gend: %d", grid[IDIR].gbeg, grid[IDIR].gend);
  B = d->Vc[iBPHI];

  DOM_LOOP(k,j,i){
    if (i == IBEG){
      Jz1D[k][j][i] = 2/(r[i]+r[i+1]) * (B[k][j][i+1]*r[i+1]-B[k][j][i]*r[i])/(r[i+1]-r[i]);
    }
    else if (i == IEND){
      Jz1D[k][j][i] = 2/(r[i-1]+r[i]) * (B[k][j][i]*r[i]-B[k][j][i-1]*r[i-1])/(r[i]-r[i-1]);
    } else {
      Jleft = 2/(r[i-1]+r[i]) * (B[k][j][i]*r[i]-B[k][j][i-1]*r[i-1])/(r[i]-r[i-1]);
      Jright = 2/(r[i]+r[i+1]) * (B[k][j][i+1]*r[i+1]-B[k][j][i]*r[i])/(r[i+1]-r[i]);
      Jz1D[k][j][i] = 0.5*(Jleft+Jright);
    }
    // Now I make Jz1D dimensional
    Jz1D[k][j][i] *= CONST_c/(4*CONST_PI)*unit_Mfield/UNIT_LENGTH;
  }

}
#endif

#if WRITE_J == YES
#if DIMENSIONS != 2
  #error ComputeJforOutput works only in 2D
#endif
void ComputeJforOutput(const Data *d, Grid *grid, double ***Ji, double ***Jj) {
  int i,j,k;
  int dir;
  double ****J_isweep, ****J_jsweep, ****J_isweep_avg, ****J_jsweep_avg;
  Data *d_temp;
  RBox box;
  
  alloc_Data(d_temp);

  J_isweep = ARRAY_4D(3,NX3_TOT, NX2_TOT, NX1_TOT, double);
  J_jsweep = ARRAY_4D(3,NX3_TOT, NX2_TOT, NX1_TOT, double);
  J_isweep_avg = ARRAY_4D(3,NX3_TOT, NX2_TOT, NX1_TOT, double);
  J_jsweep_avg = ARRAY_4D(3,NX3_TOT, NX2_TOT, NX1_TOT, double);

 // I do a TOT_LOOP as J is defined in the same way as Vc (in therms of memory),
 // see capillary_wall.c/void alloc_Data(Data *data).
 // Copy d->B[][][] inside d_temp->B..
  TOT_LOOP(k, j, i) {
    d_temp->J[IDIR][k][j][i] = d->J[IDIR][k][j][i];
    d_temp->J[JDIR][k][j][i] = d->J[JDIR][k][j][i];
    d_temp->J[KDIR][k][j][i] = d->J[KDIR][k][j][i];
  }

  // call twice GetCurrent
  GetCurrent (d_temp, IDIR, grid);
  TOT_LOOP(k, j, i) {
    J_isweep[IDIR][k][j][i] = d->J[IDIR][k][j][i];
    J_isweep[JDIR][k][j][i] = d->J[JDIR][k][j][i];
    J_isweep[KDIR][k][j][i] = d->J[KDIR][k][j][i];
  }
  // TOT_LOOP(k, j, i) J_IDIR[k][j][i] = d_temp->J[IDIR][k][j][i];
  GetCurrent (d_temp, JDIR, grid);
  TOT_LOOP(k, j, i) {
    J_jsweep[IDIR][k][j][i] = d->J[IDIR][k][j][i];
    J_jsweep[JDIR][k][j][i] = d->J[JDIR][k][j][i];
    J_jsweep[KDIR][k][j][i] = d->J[KDIR][k][j][i];
  }
  // TOT_LOOP(k, j, i) J_JDIR[k][j][i] = d_temp->J[JDIR][k][j][i];

  // Average d_temp->J inside J
  // Average along i direction
  box.ib =       0; box.ie = NX1_TOT-1-IOFFSET;
  box.jb = JOFFSET; box.je = NX2_TOT-1-JOFFSET;
  box.kb = KOFFSET; box.ke = NX3_TOT-1-KOFFSET;
  BOX_LOOP(&box,k,j,i){
    J_isweep_avg[IDIR][k][j][i] = 0.5 * (J_isweep[IDIR][k][j][i] + J_isweep[IDIR][k][j][i+1]);
    J_isweep_avg[JDIR][k][j][i] = 0.5 * (J_isweep[JDIR][k][j][i] + J_isweep[JDIR][k][j][i+1]);
  }
  // Average along j direction
  box.ib = IOFFSET; box.ie = NX1_TOT-1-IOFFSET;
  box.jb =       0; box.je = NX2_TOT-1-JOFFSET;
  box.kb = KOFFSET; box.ke = NX3_TOT-1-KOFFSET;
  BOX_LOOP(&box,k,j,i){
    J_jsweep_avg[IDIR][k][j][i] = 0.5 * (J_jsweep[IDIR][k][j][i] + J_jsweep[IDIR][k][j][i+1]);
    J_jsweep_avg[JDIR][k][j][i] = 0.5 * (J_jsweep[JDIR][k][j][i] + J_jsweep[JDIR][k][j][i+1]);
  }
  TOT_LOOP(k, j, i) {
    Ji[k][j][i] =  0.5*(J_jsweep_avg[IDIR][k][j][i] + J_isweep_avg[IDIR][k][j][i]);
    Jj[k][j][i] =  0.5*(J_jsweep_avg[JDIR][k][j][i] + J_isweep_avg[JDIR][k][j][i]);
  }

  // Deallocate all data allocated in this function
  free_Data(&d_temp);
  FreeArray4D ((void *) J_isweep);
  FreeArray4D ((void *) J_jsweep);
  FreeArray4D ((void *) J_isweep_avg);
  FreeArray4D ((void *) J_jsweep_avg);
}
#endif