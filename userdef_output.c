#include "pluto.h"
#include "gamma_transp.h"
#include "capillary_wall.h"
#include "current_table.h"
#include "prototypes.h"

#define WRITE_T_MU_NE_IONIZ YES
#define WRITE_J NO

#if WRITE_J == YES
  void ComputeJ1DforOutput(const Data *d, Grid *grid, double ***J);
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

  #if WRITE_J == YES
    double ***J;
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
  #if WRITE_J == YES
    J = GetUserVar("J");
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
  #if WRITE_J == YES
    ComputeJ1DforOutput(d, grid, J);
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

#if WRITE_J == YES
void ComputeJ1DforOutput(const Data *d, Grid *grid, double ***J){
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
      J[k][j][i] = 2/(r[i]+r[i+1]) * (B[k][j][i+1]*r[i+1]-B[k][j][i]*r[i])/(r[i+1]-r[i]);
    }
    else if (i == IEND){
      J[k][j][i] = 2/(r[i-1]+r[i]) * (B[k][j][i]*r[i]-B[k][j][i-1]*r[i-1])/(r[i]-r[i-1]);
    } else {
      Jleft = 2/(r[i-1]+r[i]) * (B[k][j][i]*r[i]-B[k][j][i-1]*r[i-1])/(r[i]-r[i-1]);
      Jright = 2/(r[i]+r[i+1]) * (B[k][j][i+1]*r[i+1]-B[k][j][i]*r[i])/(r[i+1]-r[i]);
      J[k][j][i] = 0.5*(Jleft+Jright);
    }
    // Now I make J dimensional
    J[k][j][i] *= CONST_c/(4*CONST_PI)*unit_Mfield/UNIT_LENGTH;
  }

  // box.ib =       0; box.ie = NX1_TOT-1-IOFFSET;
  // box.jb = JOFFSET; box.je = NX2_TOT-1-JOFFSET;
  // box.kb = KOFFSET; box.ke = NX3_TOT-1-KOFFSET;
  //
  // BOX_LOOP(&box,k,j,i){
  //   Jstagg[k][j][i] = 2/(r[i]+r[i+1]) * (B[k][j][i+1]*r[i+1]-B[k][j][i+1]*r[i+1])/(r[i+1]-r[i]);
  //   // print1("\nstep : i:%d, j:%d, k:%d",i,j,k);
  // }
  //
  // KTOT_LOOP(k) JTOT_LOOP(j) {
  //   J[k][j][0] = Jstagg[k][j][0];
  //   J[k][j][NX1_TOT-1] = Jstagg[k][j][NX1_TOT-1-IOFFSET];
  // }
  // BOX_LOOP(&box,k,j,i){
  //   J[k][j][i+1] = 0.5*(Jstagg[k][j][i+1] + Jstagg[k][j][i]);
  // }

}

// -------------------------------------
// int  i, j, k;
// static int first_call = 1;
// double *dx1, *dx2, *dx3;
// double *r, *rp, *th, *thp;
// double ***Bx1, ***Bx2, ***Bx3;
// double ***Jx1, ***Jx2, ***Jx3;
// double dx1_Bx2 = 0.0, dx2_Bx1 = 0.0;
// double dx2_Bx3 = 0.0, dx3_Bx2 = 0.0;
// double dx1_Bx3 = 0.0, dx3_Bx1 = 0.0;
// double d12, d21, d13, d31, d23, d32;
// double h3;
// static double ***a23Bx3, ***a13Bx3, ***a12Bx2;
// RBox box;
//
//
// B = d->Vc[iBPHI];
// dr = grid[IDIR].dx;
//
//
// if (first_call){
//   a12Bx2 = Bx2;
//   a23Bx3 = Bx3;
//   a13Bx3 = Bx3;
//   #if GEOMETRY == CYLINDRICAL && COMPONENTS == 3
//    a13Bx3 = ARRAY_3D(NX3_MAX, NX2_MAX, NX1_MAX, double);
//   first_call = 0;
// }
//
//  r   = grid[IDIR].x; rp  = grid[IDIR].xr;
//
// /* ----------------------------------------------------
//     IDIR: Compute {Jx1, Jx2, Jx3} at i+1/2,j,k faces.
//    ---------------------------------------------------- */
//
//   box.ib =       0; box.ie = NX1_TOT-1-IOFFSET;
//   box.jb = JOFFSET; box.je = NX2_TOT-1-JOFFSET;
//   box.kb = KOFFSET; box.ke = NX3_TOT-1-KOFFSET;
//
//   BOX_LOOP(&box,k,j,i){
//
//     #if GEOMETRY == CYLINDRICAL
//      d32 = d31 = 0.0;
//      d13 = 2.0/(fabs(r[i+1])*r[i+1] - fabs(r[i])*r[i]);
//     #elif GEOMETRY == POLAR
//      d23 = d21 = 1.0/(rp[i]*dx2[j]);
//      d12 = 2.0/(fabs(r[i+1])*r[i+1] - fabs(r[i])*r[i]);
//     #elif GEOMETRY == SPHERICAL
//      h3  = rp[i]*sin(th[j]);
//
//      D_EXPAND(d12 = 1.0/(rp[i]*dx1[i]); d13 = 1.0/(rp[i]*dx1[i]);   ,
//               d21 = 1.0/(rp[i]*dx2[j]); d23 = 1.0/(h3*dx2[j]);      ,
//               d32 = 1.0/(h3*dx3[k]);    d31 = 1.0/(h3*dx3[k]);)
//     #endif
//
//     #if COMPONENTS == 3
//      dx2_Bx3 = 0.5*(CDIFF_X2(a23Bx3,k,j,i) + CDIFF_X2(a23Bx3,k,j,i+1))*d23;
//      dx3_Bx2 = 0.5*(CDIFF_X3(Bx2,k,j,i)    + CDIFF_X3(Bx2,k,j,i+1)  )*d32;
//      Jx1[k][j][i] = (dx2_Bx3 - dx3_Bx2);
//
//      dx3_Bx1 = 0.5*(CDIFF_X3(Bx1,k,j,i) + CDIFF_X3(Bx1,k,j,i+1))*d31;
//      dx1_Bx3 = FDIFF_X1(a13Bx3,k,j,i)*d13;
//      Jx2[k][j][i] = (dx3_Bx1 - dx1_Bx3);
//     #endif
//
//     dx1_Bx2 = FDIFF_X1(a12Bx2,k,j,i)*d12;
//     dx2_Bx1 = 0.5*(CDIFF_X2(Bx1,k,j,i) + CDIFF_X2(Bx1,k,j,i+1))*d21;
//     Jx3[k][j][i] = (dx1_Bx2 - dx2_Bx1);
//   }



#endif
