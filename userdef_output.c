#include "pluto.h"
#include "gamma_transp.h"
#include "capillary_wall.h"
#include "current_table.h"
#include "prototypes.h"
#include "debug_utilities.h"

#define WRITE_T_MU_NE_IONIZ YES

// This is a test
void ComputeJ1DforOutput(const Data *d, Grid *grid, double ***Jz);
void ComputeJrforOutput(const Data *d, Grid *grid, double ***Jr);

/* *************************************************************** */
void ComputeUserVar (const Data *d, Grid *grid)
/*
 *
 *  PURPOSE
 *
 *    Define user-defined output variables
 *
 *
 *
 ***************************************************************** */
{
  int i, j, k;
  char Jz_name[] = "Jz";
  char Jr_name[] = "Jr";
  char interBound_name[] = "interBound";
  char etax1_name[] = "etax1";
  char knor_name[] = "knor";

  #if WRITE_T_MU_NE_IONIZ==YES
    double ***T, ***ioniz, ***ne;
    double v[NVAR]; /*[Ema] I hope that NVAR as dimension is fine!*/
    int nv;
    double mu;
  #endif

  #if MULTIPLE_GHOSTS==YES
    Data *d_corrected_r, *d_corrected_z;
    char vr_c_r_name[] = "vr_c_r", vr_c_z_name[] = "vr_c_z";
    char vz_c_r_name[] = "vz_c_r", vz_c_z_name[] = "vz_c_z"; 
    char rho_c_r_name[] = "rho_c_r", rho_c_z_name[] = "rho_c_z";
    char Bx3_c_r_name[] = "Bx3_c_r", Bx3_c_z_name[] = "Bx3_c_z";
    double unit_Mfield;
    char T_c_z_name[] = "T_c_z", T_c_r_name[] = "T_c_r";

    unit_Mfield = COMPUTE_UNIT_MFIELD(UNIT_VELOCITY, UNIT_DENSITY);
  #endif

/******************************************************/
/*I allocate space for all the variables of the output*/
/******************************************************/
  #if MULTIPLE_GHOSTS==YES
    d_corrected_r = alloc_Data();
    d_corrected_z = alloc_Data();
    // print1("\nRicordati di deallocare d_corrected_r e ..._z !");
    copy_Data_Vc(d_corrected_z, d);
    copy_Data_Vc(d_corrected_r, d);
    ApplyMultipleGhosts(d_corrected_r, IDIR);
    ApplyMultipleGhosts(d_corrected_z, JDIR);

    // I export the corrected vr, vz, rho
    if (CheckUserVar(vr_c_r_name)) {
      double ***vr_c_r;
      vr_c_r = GetUserVar(vr_c_r_name);
      DOM_LOOP (k,j,i)  vr_c_r[k][j][i] = d_corrected_r->Vc[iVR][k][j][i]*UNIT_VELOCITY;
    }
    if (CheckUserVar(vr_c_z_name)) {
      double ***vr_c_z;
      vr_c_z = GetUserVar(vr_c_z_name);
      DOM_LOOP (k,j,i)  vr_c_z[k][j][i] = d_corrected_z->Vc[iVR][k][j][i]*UNIT_VELOCITY;
    }
    if (CheckUserVar(vz_c_r_name)) {
      double ***vz_c_r;
      vz_c_r = GetUserVar(vz_c_r_name);
      DOM_LOOP (k,j,i)  vz_c_r[k][j][i] = d_corrected_r->Vc[iVZ][k][j][i]*UNIT_VELOCITY;
    }
    if (CheckUserVar(vz_c_z_name)) {
      double ***vz_c_z;
      vz_c_z = GetUserVar(vz_c_z_name);
      DOM_LOOP (k,j,i)  vz_c_z[k][j][i] = d_corrected_z->Vc[iVZ][k][j][i]*UNIT_VELOCITY;
    }
    if (CheckUserVar(rho_c_r_name)) {
      double ***rho_c_r;
      rho_c_r = GetUserVar(rho_c_r_name);
      DOM_LOOP (k,j,i)  rho_c_r[k][j][i] = d_corrected_r->Vc[RHO][k][j][i]*UNIT_DENSITY;
    }
    if (CheckUserVar(rho_c_z_name)) {
      double ***rho_c_z;
      rho_c_z = GetUserVar(rho_c_z_name);
      DOM_LOOP (k,j,i)  rho_c_z[k][j][i] = d_corrected_z->Vc[RHO][k][j][i]*UNIT_DENSITY;
    }
    if (CheckUserVar(Bx3_c_r_name)) {
      double ***Bx3_c_r;
      Bx3_c_r = GetUserVar(Bx3_c_r_name);
      DOM_LOOP (k,j,i)  Bx3_c_r[k][j][i] = d_corrected_r->Vc[BX3][k][j][i]*unit_Mfield;
    }
    if (CheckUserVar(Bx3_c_z_name)) {
      double ***Bx3_c_z;
      Bx3_c_z = GetUserVar(Bx3_c_z_name);
      DOM_LOOP (k,j,i)  Bx3_c_z[k][j][i] = d_corrected_z->Vc[BX3][k][j][i]*unit_Mfield;
    }
    if (CheckUserVar(T_c_r_name)) {
      double ***T_c_r;
      T_c_r = GetUserVar(T_c_r_name);
      DOM_LOOP (k,j,i) {
        for (nv=NVAR; nv--;) v[nv] = d_corrected_r->Vc[nv][k][j][i];
        if (GetPV_Temperature(v, &(T_c_r[k][j][i]) )!=0) {
          print1("ComputeUserVar:[Ema] Error computing temperature!");
        }
      }
    }
    if (CheckUserVar(T_c_z_name)) {
      double ***T_c_z;
      T_c_z = GetUserVar(T_c_z_name);
      DOM_LOOP (k,j,i) {
        for (nv=NVAR; nv--;) v[nv] = d_corrected_z->Vc[nv][k][j][i];
        if (GetPV_Temperature(v, &(T_c_z[k][j][i]) )!=0) {
          print1("ComputeUserVar:[Ema] Error computing temperature!");
        }
      }
    }
  #endif

  if (CheckUserVar (Jz_name)) {
    double ***Jz;
    Jz = GetUserVar(Jz_name);
    ComputeJ1DforOutput(d, grid, Jz);
  }
  if (CheckUserVar (Jr_name)) {
    double ***Jr;
    Jr = GetUserVar(Jr_name);
    ComputeJrforOutput(d, grid, Jr);
  }

  // I Export Internal boundary flags
  if (CheckUserVar(interBound_name)) {
    double ***interBound;
    interBound = GetUserVar(interBound_name);
    DOM_LOOP(k,j,i){
      if ((int) (d->flag[k][j][i] & FLAG_INTERNAL_BOUNDARY)) {
        interBound[k][j][i] = 33.3; //Just a conventional number
      } else {
        interBound[k][j][i] = 0.0; //Just a conventional number
      }
    }
  }

  if(CheckUserVar(etax1_name)) {
    double ***etax1;
    etax1 = GetUserVar(etax1_name);
    DOM_LOOP(k,j,i) {
      for (nv=NVAR; nv--;) v[nv] = d->Vc[nv][k][j][i];
      #if RESISTIVITY != NO
        Resistive_eta( v, grid[IDIR].x_glob[i], grid[JDIR].x_glob[j], grid[KDIR].x_glob[k], NULL, &(etax1[k][j][i]));
        etax1[k][j][i] *= UNIT_ETA;
      #else
        etax1[k][j][i] = 0.0;
      #endif
    }
  }

  if (CheckUserVar(knor_name)) {
    double ***knor;
    double kpar, phi;
    knor = GetUserVar(knor_name);
    DOM_LOOP(k,j,i) {
      for (nv=NVAR; nv--;) v[nv] = d->Vc[nv][k][j][i];
      #if THERMAL_CONDUCTION != NO
        TC_kappa( v, grid[IDIR].x_glob[i], grid[JDIR].x_glob[j], grid[KDIR].x_glob[k], &kpar, &(knor[k][j][i]), &phi);
        knor[k][j][i] *= UNIT_KAPPA;
      #else 
        knor[k][j][i] = 0.0;
      #endif
    }
  }

  #if WRITE_T_MU_NE_IONIZ==YES
    T = GetUserVar("T");
    #if EOS==PVTE_LAW
      ioniz = GetUserVar("ioniz");
      ne = GetUserVar("ne");
    #endif
    DOM_LOOP(k,j,i){
      #if EOS==IDEAL
        mu = MeanMolecularWeight(d->Vc);
        T[k][j][i] = d->Vc[PRS][k][j][i]/d->Vc[RHO][k][j][i]*KELVIN*mu;
      #elif EOS==PVTE_LAW
        for (nv=NVAR; nv--;) v[nv] = d->Vc[nv][k][j][i];
        if (GetPV_Temperature(v, &(T[k][j][i]) )!=0) {
          print1("ComputeUserVar:[Ema] Error computing temperature!");
        }
        // print1("\nI just assigned %g to T[%d][%d][%d] for output",T[k][j][i], k,j,i);
        GetMu(T[k][j][i], v[RHO], &mu);
        ioniz[k][j][i] = 1/mu - 1;
        ne[k][j][i] = ioniz[k][j][i] * (v[RHO]*UNIT_DENSITY) / CONST_mp;
      #endif
    }
  #endif

  #if MULTIPLE_GHOSTS==YES
    // free_Data(&d_corrected_r);
    // free_Data(&d_corrected_z);
    free_Data(d_corrected_r);
    free_Data(d_corrected_z);
  #endif
}


/* ************************************************************* */
void ChangeDumpVar ()
/*
 *
 *
 *************************************************************** */
{
  SetDumpVar("bx1", VTK_OUTPUT, NO);
  SetDumpVar("bx2", VTK_OUTPUT, NO);
  SetDumpVar("vx3", VTK_OUTPUT, NO);
}

void ComputeJ1DforOutput(const Data *d, Grid *grid, double ***Jz){
  int i, j, k;
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
      Jz[k][j][i] = -2/(r[i]+r[i+1]) * (B[k][j][i+1]*r[i+1]-B[k][j][i]*r[i])/(r[i+1]-r[i]);
    }
    else if (i == IEND){
      Jz[k][j][i] = -2/(r[i-1]+r[i]) * (B[k][j][i]*r[i]-B[k][j][i-1]*r[i-1])/(r[i]-r[i-1]);
    } else {
      Jleft = -2/(r[i-1]+r[i]) * (B[k][j][i]*r[i]-B[k][j][i-1]*r[i-1])/(r[i]-r[i-1]);
      Jright = -2/(r[i]+r[i+1]) * (B[k][j][i+1]*r[i+1]-B[k][j][i]*r[i])/(r[i+1]-r[i]);
      Jz[k][j][i] = 0.5*(Jleft+Jright);
    }
    // Now I make Jz dimensional
    Jz[k][j][i] *= CONST_c/(4*CONST_PI)*unit_Mfield/UNIT_LENGTH;
    // I put to zero the current density where I don't need it
  }
  // Now I put to 0 the Jr where I have boundary
  // (I could not do it before, as I need some values at the internal boundary
  // to compute some values near the internal boundary)
  DOM_LOOP(k,j,i)
    if ((int) (d->flag[k][j][i] & FLAG_INTERNAL_BOUNDARY))
      Jz[k][j][i] = 0.0;
}

void ComputeJrforOutput(const Data *d, Grid *grid, double ***Jr) {
  int i, j, k;
  double *z, ***B;
  double Jleft, Jright;
  // RBox box;
  double unit_Mfield;

  unit_Mfield = COMPUTE_UNIT_MFIELD(UNIT_VELOCITY, UNIT_DENSITY);

  z = grid[JDIR].x;
  // x = grid[0].x_glob;
  // print1("for 0, gbeg: %d, gend: %d", grid[IDIR].gbeg, grid[IDIR].gend);
  B = d->Vc[iBPHI];

  DOM_LOOP(k,j,i){
    if (i == IBEG){
      Jr[k][j][i] = (B[k][j+1][i]-B[k][j][i])/(z[j+1]-z[j]);
    }
    else if (i == IEND){
      Jr[k][j][i] = (B[k][j][i]-B[k][j-1][i])/(z[j+1]-z[j]);
    } else {
      Jleft = (B[k][j][i]-B[k][j-1][i])/(z[j]-z[j-1]);
      Jright = (B[k][j+1][i]-B[k][j][i])/(z[j+1]-z[j]);
      Jr[k][j][i] = 0.5*(Jleft+Jright);
    }
    // Now I make Jr dimensional
    Jr[k][j][i] *= CONST_c/(4*CONST_PI)*unit_Mfield/UNIT_LENGTH;
  }
  // Now I put to 0 the Jr where I have boundary
  // (I could not do it before, as I need some values at the internal boundary
  // to compute some values near the internal boundary)
  DOM_LOOP(k,j,i)
    if ((int) (d->flag[k][j][i] & FLAG_INTERNAL_BOUNDARY))
      Jr[k][j][i] = 0.0;
}