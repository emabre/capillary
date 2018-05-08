#include "pluto.h"
#include "gamma_transp.h"
#include "capillary_wall.h"
#include "current_table.h"
#include "prototypes.h"

#define WRITE_T_MU_NE_IONIZ YES
#define WRITE_J1D NO
#define WRITE_J2D NO
#define WRITE_J NO

#if WRITE_J1D == YES
  void ComputeJ1DforOutput(const Data *d, Grid *grid, double ***Jz);
#endif

// This is a test
#if WRITE_J == YES
  void ComputeJ1DforOutput(const Data *d, Grid *grid, double ***Jz);
  void ComputeJrforOutput(const Data *d, Grid *grid, double ***Jr);
#endif

#if WRITE_J_OLD == YES
  void ComputeJforOutput_old(const Data *d, Grid *grid, double ***Ji, double ***Jj);
#endif

#if WRITE_J2D == YES
  void ComputeJ2DforOutput(const Data *d, Grid *grid, double ***Jr, double ***Jz, double ***Jphi);
  void GetCurrentForOutput(const Data *d, int dir, Grid *grid, double ***Jr, double ***Jz, double ***Jphi);
#endif

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
  double ***interBound;

  /*[Rob] I could exploit the (I think) globally available runtime structure to get the number and name of variables
  set for output so that I do not need to modify consistently bot the pluto.ini and the definition of these macros*/
  #if WRITE_J1D == YES || WRITE_J == YES
    double ***Jz;
  #endif
  #if WRITE_J == YES
    double ***Jr;
  #endif
  #if  WRITE_J2D == YES
    double ***Jr;
    double ***Jz;
    double ***Jphi;
  #endif

  #if WRITE_T_MU_NE_IONIZ==YES
    double ***T, ***ioniz, ***mu, ***ne;
    double v[NVAR]; /*[Ema] I hope that NVAR as dimension is fine!*/
    int nv;
  #endif

  #if MULTIPLE_GHOSTS==YES
    double mu_aux; /*Auxiliary variable*/
    //Data d_corrected_r, d_corrected_z;
    Data *d_corrected_r, *d_corrected_z;
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
  #if WRITE_J1D == YES || WRITE_J == YES
    Jz = GetUserVar("Jz");
  #endif
  #if WRITE_J2D == YES
    Jr = GetUserVar("Jr");
    Jz = GetUserVar("Jz");
    Jphi = GetUserVar("Jphi");
  #endif
  #if WRITE_J == YES
    Jr = GetUserVar("Jr");
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
    // alloc_Data(&d_corrected_r);
    // alloc_Data(&d_corrected_z);
    d_corrected_r = alloc_Data();
    d_corrected_z = alloc_Data();
    // print1("\nRicordati di deallocare d_corrected_r e ..._z !");

    // copy_Data_Vc(&d_corrected_z, d);
    // copy_Data_Vc(&d_corrected_r, d);
    // ApplyMultipleGhosts(&d_corrected_r, 0);
    // ApplyMultipleGhosts(&d_corrected_z, 1);
    copy_Data_Vc(d_corrected_z, d);
    copy_Data_Vc(d_corrected_r, d);
    ApplyMultipleGhosts(d_corrected_r, 0);
    ApplyMultipleGhosts(d_corrected_z, 1);

    DOM_LOOP(k,j,i){
      // vr_c_r[k][j][i] = d_corrected_r.Vc[iVR][k][j][i]*UNIT_VELOCITY;
      // vr_c_z[k][j][i] = d_corrected_z.Vc[iVR][k][j][i]*UNIT_VELOCITY;
      // vz_c_r[k][j][i] = d_corrected_r.Vc[iVZ][k][j][i]*UNIT_VELOCITY;
      // vz_c_z[k][j][i] = d_corrected_z.Vc[iVZ][k][j][i]*UNIT_VELOCITY;
      // rho_c_r[k][j][i] = d_corrected_r.Vc[RHO][k][j][i]*UNIT_DENSITY;
      // rho_c_z[k][j][i] = d_corrected_z.Vc[RHO][k][j][i]*UNIT_DENSITY;
      vr_c_r[k][j][i] = d_corrected_r->Vc[iVR][k][j][i]*UNIT_VELOCITY;
      vr_c_z[k][j][i] = d_corrected_z->Vc[iVR][k][j][i]*UNIT_VELOCITY;
      vz_c_r[k][j][i] = d_corrected_r->Vc[iVZ][k][j][i]*UNIT_VELOCITY;
      vz_c_z[k][j][i] = d_corrected_z->Vc[iVZ][k][j][i]*UNIT_VELOCITY;
      rho_c_r[k][j][i] = d_corrected_r->Vc[RHO][k][j][i]*UNIT_DENSITY;
      rho_c_z[k][j][i] = d_corrected_z->Vc[RHO][k][j][i]*UNIT_DENSITY;
    }

    // I take the corrected T
    #if WRITE_T_MU_NE_IONIZ==YES
      DOM_LOOP(k,j,i){
        #if EOS==IDEAL
          mu_aux = MeanMolecularWeight(d_corrected_r.Vc);
          T_c_r[k][j][i] = d_corrected_r.Vc[PRS][k][j][i]/d_corrected_r.Vc[RHO][k][j][i]*KELVIN*mu_aux;
        #elif EOS==PVTE_LAW
          // for (nv=NVAR; nv--;) v[nv] = d_corrected_r.Vc[nv][k][j][i];
          for (nv=NVAR; nv--;) v[nv] = d_corrected_r->Vc[nv][k][j][i];
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
          // for (nv=NVAR; nv--;) v[nv] = d_corrected_z.Vc[nv][k][j][i];
          for (nv=NVAR; nv--;) v[nv] = d_corrected_z->Vc[nv][k][j][i];
          if (GetPV_Temperature(v, &(T_c_z[k][j][i]) )!=0) {
            print1("ComputeUserVar:[Ema] Error computing temperature!");
          }
          GetMu(T_c_z[k][j][i], v[RHO], &(mu_aux));
        #endif
      }
    #endif
  #endif
  #if WRITE_J1D == YES || WRITE_J
    ComputeJ1DforOutput(d, grid, Jz);
  #endif
  #if WRITE_J
    ComputeJrforOutput(d, grid, Jr);
  #endif

  #if WRITE_J2D == YES
    ComputeJ2DforOutput(d, grid, Jr, Jz, Jphi);
  #endif

  #if WRITE_J_OLD == YES
    ComputeJforOutput(d, grid, Jr, Jz);
  #endif


  // free_Data(&d_corrected_r);
  // free_Data(&d_corrected_z);
  free_Data(d_corrected_r);
  free_Data(d_corrected_z);
}
/* ************************************************************* */
void ChangeDumpVar ()
/*
 *
 *
 *************************************************************** */
{
  Image *image;
}

#if WRITE_J2D

  void ComputeJ2DforOutput(const Data *d, Grid *grid, double ***Jr, double ***Jz, double ***Jphi) {
    #if DIMENSIONS != 2
      #error ComputeJ2DforOutput works only in 2D
    #endif
    #if GEOMETRY != CYLINDRICAL
      #error ComputeJ2DforOutput works only for cylindrical geometry
    #endif

    int i,j,k;
    // double J_isweep[3][NX3_TOT][NX2_TOT][NX1_TOT];
    // double J_jsweep[3][NX3_TOT][NX2_TOT][NX1_TOT];
    // double J_isweep_avg[3][NX3_TOT][NX2_TOT][NX1_TOT];
    // double J_jsweep_avg[3][NX3_TOT][NX2_TOT][NX1_TOT];
    double ****J_isweep, ****J_jsweep, ****J_isweep_avg, ****J_jsweep_avg;
    double c_area_i, c_area_ip12, c_area_im12;
    RBox box;
    double unit_Mfield;

    // print1("\n ComputeJ2DforOutput:NX3_TOT=%d, NX2_TOT=%d, NX1_TOT=%d",NX3_TOT, NX2_TOT, NX1_TOT);
    J_isweep = ARRAY_4D(3, NX3_TOT, NX2_TOT, NX1_TOT, double);
    J_jsweep = ARRAY_4D(3, NX3_TOT, NX2_TOT, NX1_TOT, double);
    J_isweep_avg = ARRAY_4D(3, NX3_TOT, NX2_TOT, NX1_TOT, double);
    J_jsweep_avg = ARRAY_4D(3, NX3_TOT, NX2_TOT, NX1_TOT, double);
    TOT_LOOP(k, j, i) {
      J_isweep[IDIR][k][j][i] = 0.0;
      J_isweep[JDIR][k][j][i] = 0.0;
      J_jsweep[IDIR][k][j][i]= 0.0;
      J_jsweep[JDIR][k][j][i]= 0.0;
      J_isweep_avg[IDIR][k][j][i]= 0.0;
      J_isweep_avg[JDIR][k][j][i]= 0.0;
      J_jsweep_avg[IDIR][k][j][i]= 0.0;
      J_jsweep_avg[JDIR][k][j][i]= 0.0;
    }

    print1("MEMENTO: \nattento a come fai la media delle J, siamo in geometria cilindrica, ricontrollare!");

    // I do a TOT_LOOP as J is defined in the same way as Vc (in therms of memory),
    // see capillary_wall.c/void alloc_Data(Data *data).

      // call twice GetCurrent
      GetCurrentForOutput (d, IDIR, grid, Jr, Jz, Jphi);
      TOT_LOOP(k, j, i) {
        J_isweep[IDIR][k][j][i] = d->J[IDIR][k][j][i];
        J_isweep[JDIR][k][j][i] = d->J[JDIR][k][j][i];
        J_isweep[KDIR][k][j][i] = d->J[KDIR][k][j][i];
      }
      GetCurrentForOutput (d, JDIR, grid, Jr, Jz, Jphi);
      TOT_LOOP(k, j, i) {
        J_jsweep[IDIR][k][j][i] = d->J[IDIR][k][j][i];
        J_jsweep[JDIR][k][j][i] = d->J[JDIR][k][j][i];
        J_jsweep[KDIR][k][j][i] = d->J[KDIR][k][j][i];
      }

      // Average along i direction
      box.ib =       0; box.ie = NX1_TOT-1-IOFFSET;
      box.jb = JOFFSET; box.je = NX2_TOT-1-JOFFSET;
      box.kb = KOFFSET; box.ke = NX3_TOT-1-KOFFSET;
      BOX_LOOP(&box,k,j,i){
        // J_isweep_avg[IDIR][k][j][i+1] = 0.5 * (J_isweep[IDIR][k][j][i] + J_isweep[IDIR][k][j][i+1]);
        // J_isweep_avg[JDIR][k][j][i+1] = 0.5 * (J_isweep[JDIR][k][j][i] + J_isweep[JDIR][k][j][i+1]);
        // Coefficient proportional to the area of competence of J(j,i)
        c_area_i = grid[IDIR].xr_glob[i+1]*grid[IDIR].xr_glob[i+1] - grid[IDIR].xl_glob[i+1]*grid[IDIR].xl_glob[i+1];
        // Coefficient proportional to the area of competence of J(j,i+1/2)
        c_area_ip12 = grid[IDIR].xr_glob[i+1]*grid[IDIR].xr_glob[i+1] - grid[IDIR].x_glob[i+1]*grid[IDIR].x_glob[i+1];
        // Coefficient proportional to the area of competence of J(j,i-1/2)
        c_area_im12 = grid[IDIR].x_glob[i+1]*grid[IDIR].x_glob[i+1] - grid[IDIR].xl_glob[i+1]*grid[IDIR].xl_glob[i+1];
        // Now I do the average
        J_isweep_avg[IDIR][k][j][i+1] = (c_area_ip12*J_isweep[IDIR][k][j][i+1]-c_area_im12*J_isweep[IDIR][k][j][i])/c_area_i;
        J_isweep_avg[JDIR][k][j][i+1] = (c_area_ip12*J_isweep[JDIR][k][j][i+1]-c_area_im12*J_isweep[JDIR][k][j][i])/c_area_i;
      }
      // Average along j direction
      box.ib = IOFFSET; box.ie = NX1_TOT-1-IOFFSET;
      box.jb =       0; box.je = NX2_TOT-1-JOFFSET;
      box.kb = KOFFSET; box.ke = NX3_TOT-1-KOFFSET;
      BOX_LOOP(&box,k,j,i){
        J_jsweep_avg[IDIR][k][j+1][i] = 0.5 * (J_jsweep[IDIR][k][j][i] + J_jsweep[IDIR][k][j+1][i]);
        J_jsweep_avg[JDIR][k][j+1][i] = 0.5 * (J_jsweep[JDIR][k][j][i] + J_jsweep[JDIR][k][j+1][i]);
      }

      // Average of the two averages
      DOM_LOOP(k, j, i) {
        Jr[k][j][i] =  0.5 * (J_jsweep_avg[IDIR][k][j][i] + J_isweep_avg[IDIR][k][j][i]);
        Jz[k][j][i] =  0.5 * (J_jsweep_avg[JDIR][k][j][i] + J_isweep_avg[JDIR][k][j][i]);
      }
      // Dimensionalization
      unit_Mfield = COMPUTE_UNIT_MFIELD(UNIT_VELOCITY, UNIT_DENSITY);
      DOM_LOOP(k, j, i) {
        Jr[k][j][i] *= CONST_c/(4*CONST_PI)*unit_Mfield/UNIT_LENGTH;
        Jz[k][j][i] *= CONST_c/(4*CONST_PI)*unit_Mfield/UNIT_LENGTH;
      }

      // I put to zero the current density where I don't need it
      DOM_LOOP(k,j,i){
      if ((int) (d->flag[k][j][i] & FLAG_INTERNAL_BOUNDARY)) {
        Jr[k][j][i] = 0.0;
        Jz[k][j][i] = 0.0;
      }
  }
  }

  /*-------------------------------------------------------------------------
  void GetCurrentForOutput (..)
  What follows has been copied from GetCurrent() removing parts which did not
  refer to 2D+CYLINDRICAL GEOMETRY + CELL CENTERED MHD,
  I also removed if cycle on the allocation of ***a23Bx3, ***a13Bx3, ***a13Bx3 and
  changed their names to ***a23Bx3_local, ***a13Bx3_local, ***a13Bx3_local (in order not
  to confuse/mess up the global variables ***a23Bx3, ***a13Bx3, ***a12Bx3)
  ----------------------------------------------------------------------------*/
  void GetCurrentForOutput (const Data *d, int dir, Grid *grid, double ***Jr, double ***Jz, double ***Jphi)
  {
    int  i, j, k;
    double *dx1, *dx2, *dx3;
    double *r, *rp;
    double ***Bx1, ***Bx2, ***Bx3;
    double ***Jx1, ***Jx2, ***Jx3;
    double dx1_Bx2 = 0.0, dx2_Bx1 = 0.0;
    double dx2_Bx3 = 0.0, dx3_Bx2 = 0.0;
    double dx1_Bx3 = 0.0, dx3_Bx1 = 0.0;
    double d12, d21, d13, d31, d23, d32;
    double ***a23Bx3_local, ***a13Bx3_local, ***a12Bx3_local;
    RBox box;

  /* ------------------------------------------------------------------
    1. Set pointer shortcuts.
        The primary magnetic field will be the cell-centered field.

        Note: in 2+1/2 dimensions, we also need Jx1 and Jx2 which
              contribute to the total energy equation but not
              to the induction equation.
              For this reason, we set  Bx3 to be equal to the \b
              cell centered value.
    ---------------------------------------------------------------- */

    EXPAND(Bx1 = d->Vc[BX1];    ,
          Bx2 = d->Vc[BX2];    ,
          Bx3 = d->Vc[BX3];)

    //[Ema] Just because I am lazy and I do not want to change the names in the whol function
    Jx1 = Jr;
    Jx2 = Jz;
    Jx3 = Jphi;

    dx1 = grid[IDIR].dx;
    dx2 = grid[JDIR].dx;
    dx3 = grid[KDIR].dx;

  /* ----------------------------------------------------------------
    2. Allocate static memory areas for the three arrays
        a13Bx3_local, a23Bx3_local, a12Bx3_local if necessary. Otherwise, they will
        just be shortcuts to the default magnetic field arrays.
    ---------------------------------------------------------------- */

      a12Bx3_local = Bx2;
      a23Bx3_local = Bx3;
      // a13Bx3_local = Bx3;
      // #if GEOMETRY == CYLINDRICAL && COMPONENTS == 3
      // cambia, possibilemnte usa un vettore allocato tipo array[N][N][N]
      a13Bx3_local = ARRAY_3D(NX3_MAX, NX2_MAX, NX1_MAX, double);
      // #endif
  /* --------------------------------------------------------------
    3. Construct goemetrical coefficients and arrays
    -------------------------------------------------------------- */

    // #if GEOMETRY == CYLINDRICAL
    r   = grid[IDIR].x; rp  = grid[IDIR].xr;
    // #if COMPONENTS == 3
      TOT_LOOP(k,j,i) a13Bx3_local[k][j][i] = Bx3[k][j][i]*fabs(r[i]);
    // #endif
    // #endif

  /* ------------------------------------------------------------- */
  /*!
    4b. For cell-centered MHD, we compute the three components of
          currents at a cell interface.
    -------------------------------------------------------------- */

    if (dir == IDIR){

    /* ----------------------------------------------------
        IDIR: Compute {Jx1, Jx2, Jx3} at i+1/2,j,k faces.
      ---------------------------------------------------- */

      box.ib =       0; box.ie = NX1_TOT-1-IOFFSET;
      box.jb = JOFFSET; box.je = NX2_TOT-1-JOFFSET;
      box.kb = KOFFSET; box.ke = NX3_TOT-1-KOFFSET;

      BOX_LOOP(&box,k,j,i){

        d12 = d13 = 1.0/dx1[i];
        d21 = d23 = 1.0/dx2[j];
        d31 = d32 = 1.0/dx3[k];

        // #if GEOMETRY == CYLINDRICAL
        d32 = d31 = 0.0;
        d13 = 2.0/(fabs(r[i+1])*r[i+1] - fabs(r[i])*r[i]);
        // #endif

        // #if COMPONENTS == 3
        dx2_Bx3 = 0.5*(CDIFF_X2(a23Bx3_local,k,j,i) + CDIFF_X2(a23Bx3_local,k,j,i+1))*d23;
        dx3_Bx2 = 0.5*(CDIFF_X3(Bx2,k,j,i)    + CDIFF_X3(Bx2,k,j,i+1)  )*d32; // Useless, it's 0
        Jx1[k][j][i] = (dx2_Bx3 - dx3_Bx2);

        dx3_Bx1 = 0.5*(CDIFF_X3(Bx1,k,j,i) + CDIFF_X3(Bx1,k,j,i+1))*d31; // Useless, it's 0
        dx1_Bx3 = FDIFF_X1(a13Bx3_local,k,j,i)*d13;
        Jx2[k][j][i] = (dx3_Bx1 - dx1_Bx3);
        // #endif

        dx1_Bx2 = FDIFF_X1(a12Bx3_local,k,j,i)*d12;
        dx2_Bx1 = 0.5*(CDIFF_X2(Bx1,k,j,i) + CDIFF_X2(Bx1,k,j,i+1))*d21;
        Jx3[k][j][i] = (dx1_Bx2 - dx2_Bx1);
      }

    }else if (dir == JDIR){

    /* ----------------------------------------------------
        JDIR: Compute {Jx1, Jx2, Jx3} at i,j+1/2,k faces.
      ---------------------------------------------------- */

      box.ib = IOFFSET; box.ie = NX1_TOT-1-IOFFSET;
      box.jb =       0; box.je = NX2_TOT-1-JOFFSET;
      box.kb = KOFFSET; box.ke = NX3_TOT-1-KOFFSET;

      BOX_LOOP(&box,k,j,i){
        d12 = d13 = 1.0/dx1[i];
        d21 = d23 = 1.0/dx2[j];
        d31 = d32 = 1.0/dx3[k];

        // #if GEOMETRY == CYLINDRICAL
        d32 = d31 = 0.0; /* axisymmetry */
        d13 = 1.0/(r[i]*dx1[i]);
        // #endif

        // #if COMPONENTS == 3
        dx2_Bx3 = FDIFF_X2(a23Bx3_local,k,j,i)*d23;
        dx3_Bx2 = 0.5*(CDIFF_X3(Bx2,k,j,i) + CDIFF_X3(Bx2,k,j+1,i))*d32;
        Jx1[k][j][i] = (dx2_Bx3 - dx3_Bx2);

        dx3_Bx1 = 0.5*(CDIFF_X3(Bx1,k,j,i)    + CDIFF_X3(Bx1,k,j+1,i))*d31;
        dx1_Bx3 = 0.5*(CDIFF_X1(a13Bx3_local,k,j,i) + CDIFF_X1(a13Bx3_local,k,j+1,i))*d13;
        Jx2[k][j][i] = (dx3_Bx1 - dx1_Bx3);
        // #endif
        dx1_Bx2 = 0.5*(CDIFF_X1(a12Bx3_local,k,j,i) + CDIFF_X1(a12Bx3_local,k,j+1,i))*d12;
        dx2_Bx1 = FDIFF_X2(Bx1,k,j,i)*d21;
        Jx3[k][j][i] = (dx1_Bx2 - dx2_Bx1);
      }

    }else if (dir == KDIR) {

    /* ----------------------------------------------------
        KDIR: Compute {Jx1, Jx2, Jx3} at i,j,k+1/2 faces.
      ---------------------------------------------------- */

      box.ib = IOFFSET; box.ie = NX1_TOT-1-IOFFSET;
      box.jb = JOFFSET; box.je = NX2_TOT-1-JOFFSET;
      box.kb =       0; box.ke = NX3_TOT-1-KOFFSET;

      BOX_LOOP(&box,k,j,i){
        d12 = d13 = 1.0/dx1[i];
        d21 = d23 = 1.0/dx2[j];
        d31 = d32 = 1.0/dx3[k];

        // #if COMPONENTS == 3
        dx2_Bx3 = 0.5*(CDIFF_X2(a23Bx3_local,k,j,i) + CDIFF_X2(a23Bx3_local,k+1,j,i))*d23;
        dx3_Bx2 = FDIFF_X3(Bx2,k,j,i)*d32;
        Jx1[k][j][i] = (dx2_Bx3 - dx3_Bx2);

        dx3_Bx1 = FDIFF_X3(Bx1,k,j,i)*d31;
        dx1_Bx3 = 0.5*(CDIFF_X1(a13Bx3_local,k,j,i) + CDIFF_X1(a13Bx3_local,k+1,j,i))*d13;
        Jx2[k][j][i] = (dx3_Bx1 - dx1_Bx3);
        // #endif
        dx1_Bx2 = 0.5*(CDIFF_X1(a12Bx3_local,k,j,i) + CDIFF_X1(a12Bx3_local,k+1,j,i))*d12;
        dx2_Bx1 = 0.5*(CDIFF_X2(Bx1,k,j,i)    + CDIFF_X2(Bx1,k+1,j,i))*d21;
        Jx3[k][j][i] = (dx1_Bx2 - dx2_Bx1);
      }
    }

    // [Ema]I free the array, as I allocate it every time
    FreeArray3D ((void *) a13Bx3_local);
  }
#endif

#if WRITE_J1D == YES || WRITE_J == YES
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
#endif
#if WRITE_J == YES
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
#endif

#if WRITE_J_OLD ==  YES
#if DIMENSIONS != 2
  #error ComputeJforOutput works only in 2D
#endif
void ComputeJforOutput_old(const Data *d, Grid *grid, double ***Ji, double ***Jj) {
  int i,j,k;
  // int dir;
  // double ****J_isweep, ****J_jsweep, ****J_isweep_avg, ****J_jsweep_avg;
  double J_isweep[3][NX3_TOT][NX2_TOT][NX1_TOT];
  double J_jsweep[3][NX3_TOT][NX2_TOT][NX1_TOT];
  double J_isweep_avg[3][NX3_TOT][NX2_TOT][NX1_TOT];
  double J_jsweep_avg[3][NX3_TOT][NX2_TOT][NX1_TOT];
  Data* d_temp;
  RBox box;

  // alloc_Data(&d_temp);
  d_temp = alloc_Data();

  // J_isweep = ARRAY_4D(3,NX3_TOT, NX2_TOT, NX1_TOT, double);
  // J_jsweep = ARRAY_4D(3,NX3_TOT, NX2_TOT, NX1_TOT, double);
  // J_isweep_avg = ARRAY_4D(3,NX3_TOT, NX2_TOT, NX1_TOT, double);
  // J_jsweep_avg = ARRAY_4D(3,NX3_TOT, NX2_TOT, NX1_TOT, double);

 // I do a TOT_LOOP as J is defined in the same way as Vc (in therms of memory),
 // see capillary_wall.c/void alloc_Data(Data *data).
 // Copy d->B[][][] inside d_temp->B..
  TOT_LOOP(k, j, i) {
    // d_temp.J[IDIR][k][j][i] = d->J[IDIR][k][j][i];
    // d_temp.J[JDIR][k][j][i] = d->J[JDIR][k][j][i];
    // d_temp.J[KDIR][k][j][i] = d->J[KDIR][k][j][i];
    d_temp->J[IDIR][k][j][i] = d->J[IDIR][k][j][i];
    d_temp->J[JDIR][k][j][i] = d->J[JDIR][k][j][i];
    d_temp->J[KDIR][k][j][i] = d->J[KDIR][k][j][i];
  }

  // call twice GetCurrent
  // GetCurrent (&d_temp, IDIR, grid);
  GetCurrent (d_temp, IDIR, grid);
  TOT_LOOP(k, j, i) {
    J_isweep[IDIR][k][j][i] = d->J[IDIR][k][j][i];
    J_isweep[JDIR][k][j][i] = d->J[JDIR][k][j][i];
    J_isweep[KDIR][k][j][i] = d->J[KDIR][k][j][i];
  }
  // GetCurrent (&d_temp, JDIR, grid);
  GetCurrent (d_temp, JDIR, grid);
  TOT_LOOP(k, j, i) {
    J_jsweep[IDIR][k][j][i] = d->J[IDIR][k][j][i];
    J_jsweep[JDIR][k][j][i] = d->J[JDIR][k][j][i];
    J_jsweep[KDIR][k][j][i] = d->J[KDIR][k][j][i];
  }

  // Average d_temp->J inside J
  // Average along i direction
  box.ib =       0; box.ie = NX1_TOT-1-IOFFSET;
  box.jb = JOFFSET; box.je = NX2_TOT-1-JOFFSET;
  box.kb = KOFFSET; box.ke = NX3_TOT-1-KOFFSET;
  BOX_LOOP(&box,k,j,i){
    J_isweep_avg[IDIR][k][j][i+1] = 0.5 * (J_isweep[IDIR][k][j][i] + J_isweep[IDIR][k][j][i+1]);
    J_isweep_avg[JDIR][k][j][i+1] = 0.5 * (J_isweep[JDIR][k][j][i] + J_isweep[JDIR][k][j][i+1]);
  }
  // Average along j direction
  box.ib = IOFFSET; box.ie = NX1_TOT-1-IOFFSET;
  box.jb =       0; box.je = NX2_TOT-1-JOFFSET;
  box.kb = KOFFSET; box.ke = NX3_TOT-1-KOFFSET;
  BOX_LOOP(&box,k,j,i){
    J_jsweep_avg[IDIR][k][j+1][i] = 0.5 * (J_jsweep[IDIR][k][j][i] + J_jsweep[IDIR][k][j+1][i]);
    J_jsweep_avg[JDIR][k][j+1][i] = 0.5 * (J_jsweep[JDIR][k][j][i] + J_jsweep[JDIR][k][j+1][i]);
  }

  // Average of the two averages
  DOM_LOOP(k, j, i) {
    Ji[k][j][i] =  0.5 * (J_jsweep_avg[IDIR][k][j][i] + J_isweep_avg[IDIR][k][j][i]);
    Jj[k][j][i] =  0.5 * (J_jsweep_avg[JDIR][k][j][i] + J_isweep_avg[JDIR][k][j][i]);
  }

  // Deallocate all data allocated in this function
  // free_Data(&d_temp);
  free_Data(d_temp);
  // FreeArray4D ((void *) J_isweep);
  // FreeArray4D ((void *) J_jsweep);
  // FreeArray4D ((void *) J_isweep_avg);
  // FreeArray4D ((void *) J_jsweep_avg);
}
#endif