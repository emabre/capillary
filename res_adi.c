#include <math.h>
#include "pluto.h"
#include "adi.h"
#include "capillary_wall.h"
#include "current_table.h"
#include "debug_utilities.h"

#define UNUSED(x) (void)(x)

// #define HARMAVG(a,b) ( 1/(0.5*( 1/(a) + 1/(b) )) )

double curr = 0;

double GetCurrADI() {
  return curr;
}

#if RESISTIVITY == ALTERNATING_DIRECTION_IMPLICIT
/****************************************************************************
Function to build the Ip,Im,Jp,Jm for the electrical resistivity problem
(**useless parameter is intentionally unused, to make this function suitable for a pointer
 which also wants that parameter)
*****************************************************************************/
void BuildIJ_Res (const Data *d, Grid *grid, Lines *lines,
                  double **Ip, double **Im, double **Jp,
                  double **Jm, double **CI, double **CJ, double **useless) {

  static int first_call=1;
  static double **protoIp, **protoIm, **protoJp, **protoJm, **protoCI, **protoCJ;
  int i,j,k;
  int Nlines, lidx, ridx;
  int nv, l;
  double eta[3]; // Electr. resistivity
  double v[NVAR];
  double ****Vc;
  double *inv_dri, *inv_dzi, *inv_dr, *inv_dz, *r_1;
  double *zL, *zR;
  double *rL, *rR;
  double *r, *z, *theta;
  UNUSED(**useless);

  /* -- set a pointer to the primitive vars array --
    I do this because it is done also in other parts of the code
    maybe it makes the program faster or just easier to write/read...*/
  Vc = d->Vc;

  theta = grid[KDIR].x;

  r = grid[IDIR].x;
  z = grid[JDIR].x;
  rL = grid[IDIR].xl;
  rR = grid[IDIR].xr;
  zL = grid[JDIR].xl;
  zR = grid[JDIR].xr;

  // I compose the grid-related part once forever (in static variables) and only update eta
  if (first_call) {
    // Name shorthands
    inv_dzi = grid[JDIR].inv_dxi; // reciprocal of cell spacing between centers (see set_geometry.c)
    inv_dri = grid[IDIR].inv_dxi;
    inv_dz = grid[JDIR].inv_dx; // reciprocal of cell spacing between interfaces (see set_geometry.c)
    inv_dr = grid[IDIR].inv_dx;
    r_1 = grid[IDIR].r_1;

    /*Allocatin of memory for the proto variables*/
    protoIp = ARRAY_2D(NX2_TOT, NX1_TOT, double);
    protoIm = ARRAY_2D(NX2_TOT, NX1_TOT, double);
    protoJp = ARRAY_2D(NX2_TOT, NX1_TOT, double);
    protoJm = ARRAY_2D(NX2_TOT, NX1_TOT, double);
    protoCI = ARRAY_2D(NX2_TOT, NX1_TOT, double);
    protoCJ = ARRAY_2D(NX2_TOT, NX1_TOT, double);

    /*[Opt] This is probably useless, it is here just for debugging purposes*/
    TOT_LOOP(k, j, i) {
      protoIp[j][i] = 0.0;
      protoIm[j][i] = 0.0;
      protoJp[j][i] = 0.0;
      protoJm[j][i] = 0.0;
      protoCI[j][i] = 0.0;
      protoCJ[j][i] = 0.0;
    }
    // I build the grid-related part of Im/p Jm/p, CI, CJ
    KDOM_LOOP(k) {
      LINES_LOOP(lines[IDIR], l, j, i) {
        /* :::: Ip :::: */
        protoIp[j][i] = inv_dr[i]*inv_dri[i]/rR[i];
        /* :::: Im :::: */
        if (rL[i]!=0.0)
          protoIm[j][i] = inv_dr[i]*inv_dri[i-1]/rL[i];
        else
          protoIm[j][i] = 1/(r[i]*r[i])/rR[i];
        /* :::: Jp :::: */
        protoJp[j][i] = inv_dz[j]*inv_dzi[j];
        /* :::: Jm :::: */
        protoJm[j][i] = inv_dz[j]*inv_dzi[j-1];
        /* :::: CI :::: */
        protoCI[j][i] = r_1[i];
        /* :::: CJ :::: */
        protoCJ[j][i] = 1.0;
      }
    }
    first_call = 0;
  }

  KDOM_LOOP(k) {
    
    /* :::: Ip and Im :::: */
    Nlines = lines[IDIR].N;
    for (l = 0; l<Nlines; l++) {
      j = lines->dom_line_idx[l];
      lidx = lines[IDIR].lidx[l];
      ridx = lines[IDIR].ridx[l];

      /* :::: Im on i=lidx :::: */
      for (nv=0; nv<NVAR; nv++)
        v[nv] = 0.5 * (Vc[nv][k][j][lidx] + Vc[nv][k][j][lidx-1]);
      Resistive_eta( v, rL[lidx], z[j], theta[k], NULL, eta);  // rL, z
      if (eta[0] != eta[1] || eta[1] != eta[2]) {
        ComplainAnisotropic(v, eta, rL[lidx], z[j], theta[k]);
        QUIT_PLUTO(1);
      }
      Im[j][lidx] = eta[0]*protoIm[j][lidx];

      /* :::: Ip and Im for internal (non boundary) interfaces :::: */
      for (i=lidx; i<ridx; i++) {
        for (nv=0; nv<NVAR; nv++)
          v[nv] = 0.5 * (Vc[nv][k][j][i] + Vc[nv][k][j][i+1]);
        Resistive_eta( v, rR[i], z[j], theta[k], NULL, eta);   // rR, z
        if (eta[0] != eta[1] || eta[1] != eta[2]) {
          ComplainAnisotropic(v, eta, rR[i], z[j], theta[k]);
          QUIT_PLUTO(1);
        }
        Ip[j][i] = eta[0]*protoIp[j][i];
        Im[j][i+1] = eta[0]*protoIm[j][i+1];
      }

      /* :::: Ip on i=ridx :::: */
      for (nv=0; nv<NVAR; nv++)
        v[nv] = 0.5 * (Vc[nv][k][j][ridx] + Vc[nv][k][j][ridx+1]);
      Resistive_eta( v, rR[ridx], z[j], theta[k], NULL, eta);   // rR, z
      if (eta[0] != eta[1] || eta[1] != eta[2]) {
        ComplainAnisotropic(v, eta, rR[ridx], z[j], theta[k]);
        QUIT_PLUTO(1);
      }
      Ip[j][ridx] = eta[0]*protoIp[j][ridx];
    }

    /* :::: Jp and Jm :::: */
    Nlines = lines[JDIR].N;
    for (l = 0; l<Nlines; l++) {
      i = lines[JDIR].dom_line_idx[l];
      lidx = lines[JDIR].lidx[l];
      ridx = lines[JDIR].ridx[l];

      /* :::: Jm on j=lidx :::: */
      for (nv=0; nv<NVAR; nv++)
        v[nv] = 0.5 * (Vc[nv][k][lidx][i] + Vc[nv][k][lidx-1][i]);
      Resistive_eta( v, r[i], zL[lidx], theta[k], NULL, eta);
      if (eta[0] != eta[1] || eta[1] != eta[2]) {
        ComplainAnisotropic(v, eta, r[i], zL[lidx], theta[k]);
        QUIT_PLUTO(1);
      }
      Jm[lidx][i] = eta[0]*protoJm[lidx][i];

      /* :::: Jp and Jm for internal (non boundary) interfaces :::: */
      for (j=lidx; j<ridx; j++) {
        for (nv=0; nv<NVAR; nv++)
          v[nv] = 0.5 * (Vc[nv][k][j][i] + Vc[nv][k][j+1][i]);
        Resistive_eta( v, r[i], zR[j], theta[k], NULL, eta);
        if (eta[0] != eta[1] || eta[1] != eta[2]) {
          ComplainAnisotropic(v, eta, r[i], zR[j], theta[k]);
          QUIT_PLUTO(1);
        }
        Jp[j][i] = eta[0]*protoJp[j][i];
        Jm[j+1][i] = eta[0]*protoJm[j+1][i];  
      }

      /* :::: Jp on j=ridx :::: */
      for (nv=0; nv<NVAR; nv++)
        v[nv] = 0.5 * (Vc[nv][k][ridx][i] + Vc[nv][k][ridx+1][i]);
      Resistive_eta( v, r[i], zR[ridx], theta[k], NULL, eta);
      if (eta[0] != eta[1] || eta[1] != eta[2]) {
        ComplainAnisotropic(v, eta, r[i], zR[ridx], theta[k]);
        QUIT_PLUTO(1);
      }
      Jp[ridx][i] = eta[0]*protoJp[ridx][i];
    }

    // I separate the computation of CI and CJ just to improve code readability,
    // I could also inglobate them in the previous cycles
    LINES_LOOP(lines[IDIR], l, j, i) {
      /* :::: CI :::: */
      CI[j][i] = protoCI[j][i];
      /* :::: CJ :::: */
      CJ[j][i] = protoCJ[j][i];
    }
  }

  #ifdef DEBUG_BUILDIJ
    printf("\n[BuildIJ_Res] Step: %ld", g_stepNumber);
    printf("\n[BuildIJ_Res] Im:");
    printmat(Im, NX2_TOT, NX1_TOT);
    printf("\n[BuildIJ_Res] Ip:");
    printmat(Ip, NX2_TOT, NX1_TOT);
    printf("\n[BuildIJ_Res] Jm:");
    printmat(Jm, NX2_TOT, NX1_TOT);
    printf("\n[BuildIJ_Res] Jp:");
    printmat(Jp, NX2_TOT, NX1_TOT);
  #endif
}

/****************************************************************************
Function to build the a matrix which contain the amount of increase of the
energy due to joule effect and magnetic field energy (flux of poynting vector due to
resistive magnetic diffusion)
[Opt]: Note that in the actual implementation this function needs to take as
input both the Hp_B and the Hm_B. But for the whole internal (i.e. boudary excluded)
only one among Hp_B and Hm_B is necessary. The other one is used to fill F
at one side (left or right) of the domain.
*****************************************************************************/
void ResEnergyIncrease(double **dUres, double** Hp_B, double** Hm_B, double **Br,
                       Grid *grid, Lines *lines,
                       int compute_inflow, double *inflow,
                       double dt, int dir){
  /* F :Power flux flowing from cell (i,j) to (i+1,j), when dir==IDIR;
        or from cell (i,j) to (i,j+1), when dir == JDIR.
  */
//  modifica questa funzinoe così: non usare più le Bcs(d altra parte non ha senso! questo è un termine sorgente, non un equazione da risolvere) per calcolare i valori al bordo,
//  usa i valori salvati *Br
  static double **F;
  double *dr, *dz;
  int i,j,l;
  int lidx, ridx;
  int Nlines = lines->N;
  int static first_call = 1;
  double *dV, *inv_dz, *r_1, *r;
  double *rL, *rR;
  Bcs *rbound, *lbound;
  double vol_lidx, vol_ridx;

  /*[Opt] Maybe I could do that it allocates static arrays with size NMAX_POINT (=max(NX1_TOT,NX2_TOT)) ?*/
  if (first_call) {
    /* I define it 2d in case I need to export it later*/
    F = ARRAY_2D(NX2_TOT, NX1_TOT, double);
    /*This is useless, it's just for debugging purposes*/
    ITOT_LOOP(i)
      JTOT_LOOP(j)
        F[j][i] = 0.0;
    first_call = 0;
  }
  /*This is useless, it's just for debugging purposes*/
  ITOT_LOOP(i)
    JTOT_LOOP(j)
      dUres[j][i] = 0.0;

  lbound = lines->lbound[BDIFF];
  rbound = lines->rbound[BDIFF];
  dr = grid[IDIR].dx;
  r_1 = grid[IDIR].r_1;
  dz = grid[JDIR].dx;
  rR = grid[IDIR].xr;
  rL = grid[IDIR].xl;

  if (dir == IDIR) {
    dV = grid[IDIR].dV;
    r = grid[IDIR].x_glob;

    for (l = 0; l<Nlines; l++) {
      j = lines->dom_line_idx[l];
      lidx = lines->lidx[l];
      ridx = lines->ridx[l];
      for (i=lidx; i<=ridx; i++){ // I start from lidx because I must treat carefully the lidx interface (at i=lidx-1), since there could be the domain axis
        // [Err] Decomment next line (original)
        F[j][i] = -Hp_B[j][i] * (Br[j][i+1] - Br[j][i])*dr[i] * 0.5*(Br[j][i+1]*r_1[i+1] + Br[j][i]*r_1[i]);
        // [Err] Delete next line (test)
        // F[j][i] = -Hp_B[j][i] * (Br[j][i+1] - Br[j][i])*dr[i] * 0.5*(Br[j][i+1] + Br[j][i])/rR[i];
      }
      /* I try to guess if the lower boundary in dir IDIR is the domain axis, if so I compute F consistently (with the usual formula
      I would get a division by zero) */
      if (lbound[l].kind == DIRICHLET && fabs(rL[lidx]) < 1e-20  && fabs(lbound[l].values[0]) < 1e-20) {
        F[j][lidx-1] = 0.0;
      } else {
        F[j][lidx-1] = -Hm_B[j][lidx] * (Br[j][lidx] - Br[j][lidx-1])*dr[lidx] * 0.5*(Br[j][lidx]*r_1[lidx] + Br[j][lidx-1]*r_1[lidx-1]);
      }
      // Build dU
      for (i=lidx; i<=ridx; i++)
        dUres[j][i] = -(rR[i]*F[j][i] - rL[i]*F[j][i-1])*dt/dV[i];

      if (compute_inflow) {
        /* --- I compute the inflow (energy entering from boundary) ---*/
        // Old
        // *inflow += F[j][lidx-1] * 2*CONST_PI*rL[lidx]*dz[j] * dt;
        // *inflow += -F[j][ridx] * 2*CONST_PI*rR[ridx]*dz[j] * dt;
        // Modified 27/11/2018
        vol_lidx = CONST_PI*(rR[lidx]*rR[lidx] - rL[lidx]*rL[lidx])*dz[j];
        vol_ridx = CONST_PI*(rR[ridx]*rR[ridx] - rL[ridx]*rL[ridx])*dz[j];
        *inflow += rL[lidx]*F[j][lidx-1]*dt/dV[lidx] * vol_lidx;
        *inflow += -rR[ridx]*F[j][ridx]*dt/dV[ridx] * vol_ridx;
      }
    }

  } else if (dir == JDIR) {
    inv_dz = grid[JDIR].inv_dx;

    for (l = 0; l<Nlines; l++) {
      i = lines->dom_line_idx[l];
      lidx = lines->lidx[l];
      ridx = lines->ridx[l];
      
      F[lidx][i] = -Hm_B[lidx][i] * (Br[lidx][i] - Br[lidx-1][i])*dz[lidx]*r_1[i]*r_1[i] * 0.5*(Br[lidx][i] + Br[lidx-1][i]);
      for (j=lidx; j<=ridx; j++){
        // [Err] Decomment next line
        F[j][i] = -Hp_B[j][i] * (Br[j+1][i] - Br[j][i])*dz[j]*r_1[i]*r_1[i] * 0.5*(Br[j+1][i] + Br[j][i]);
      }
      // Build dU
      for (j=lidx; j<=ridx; j++)
        dUres[j][i] = -(F[j][i] - F[j-1][i])*dt*inv_dz[j];
      
      if (compute_inflow) {
        /* --- I compute the inflow (energy entering from boundary) ---*/
        // Old
        // *inflow += F[lidx-1][i] * CONST_PI*(rR[i]*rR[i]-rL[i]*rL[i]) * dt;
        // *inflow += -F[ridx][i] * CONST_PI*(rR[i]*rR[i]-rL[i]*rL[i]) * dt;
        // Modified 27/11/2018
        vol_lidx = CONST_PI*(rR[i]*rR[i] - rL[i]*rL[i])*dz[lidx];
        vol_ridx = CONST_PI*(rR[i]*rR[i] - rL[i]*rL[i])*dz[ridx];
        *inflow += F[lidx-1][i]*dt*inv_dz[lidx] * vol_lidx;
        *inflow += -F[ridx][i]*dt*inv_dz[ridx] * vol_ridx;
      }
    }
  }
}

/****************************************************************************
Function to build the a matrix which contain the amount of increase of the
energy due to joule effect and magnetic field energy (flux of poynting vector due to
resistive magnetic diffusion), variant for DouglasRachford
[Opt]: Note that in the actual implementation this function needs to take as
input both the Hp_B and the Hm_B. But for the whole internal (i.e. boudary excluded)
only one among Hp_B and Hm_B is necessary. The other one is used to fill F
at one side (left or right) of the domain.
*****************************************************************************/
void ResEnergyIncreaseDR (double **dUres, double** Hp_B, double** Hm_B,
                          double **Br, double **Br_hat,
                          Grid *grid, Lines *lines, double dt, int dir){
  /* F :Power flux flowing from cell (i,j) to (i+1,j), when dir==IDIR;
        or from cell (i,j) to (i,j+1), when dir == JDIR.
  */
//  modifica questa funzinoe così: non usare più le Bcs(d altra parte non ha senso! questo è un termine sorgente, non un equazione da risolvere) per calcolare i valori al bordo,
//  usa i valori salvati *Br
  static double **F;
  double *dr, *dz;
  int i,j,l;
  int lidx, ridx;
  int Nlines = lines->N;
  int static first_call = 1;
  double *dV, *inv_dz, *r_1, *r;
  double *rL, *rR;
  Bcs *rbound, *lbound;

  /*[Opt] Maybe I could do that it allocates static arrays with size NMAX_POINT (=max(NX1_TOT,NX2_TOT)) ?*/
  if (first_call) {
    /* I define it 2d in case I need to export it later*/
    F = ARRAY_2D(NX2_TOT, NX1_TOT, double);
    /*This is useless, it's just for debugging purposes*/
    ITOT_LOOP(i)
      JTOT_LOOP(j)
        F[j][i] = 0.0;
    first_call = 0;
  }
  /*This is useless, it's just for debugging purposes*/
  ITOT_LOOP(i)
    JTOT_LOOP(j)
      dUres[j][i] = 0.0;

  lbound = lines->lbound[BDIFF];
  rbound = lines->rbound[BDIFF];
  dr = grid[IDIR].dx;
  r_1 = grid[IDIR].r_1;

  if (dir == IDIR) {
    rR = grid[IDIR].xr;
    rL = grid[IDIR].xl;
    dV = grid[IDIR].dV;
    r = grid[IDIR].x_glob;

    for (l = 0; l<Nlines; l++) {
      j = lines->dom_line_idx[l];
      lidx = lines->lidx[l];
      ridx = lines->ridx[l];

      for (i=lidx; i<=ridx; i++){ // I start from lidx because I must treat carefully the lidx interface (at i=lidx-1), since there could be the domain axis
        F[j][i] = -Hp_B[j][i] * (Br_hat[j][i+1] - Br_hat[j][i])*dr[i] * 0.5*(Br[j][i+1]*r_1[i+1] + Br[j][i]*r_1[i]);
      }
      /* I try to guess if the lower boundary in dir IDIR is the domain axis, if so I compute F consistently (with the usual formula
      I would get a division by zero) */
      if (lbound[l].kind == DIRICHLET && fabs(rL[lidx]) < 1e-20  && fabs(lbound[l].values[0]) < 1e-20) {
        F[j][lidx-1] = 0.0;
      } else {
        F[j][lidx-1] = -Hm_B[j][lidx] * (Br_hat[j][lidx] - Br_hat[j][lidx-1])*dr[lidx] * 0.5*(Br[j][lidx]*r_1[lidx] + Br[j][lidx-1]*r_1[lidx-1]);
      }

      // Build dU
      for (i=lidx; i<=ridx; i++)
        dUres[j][i] = -(rR[i]*F[j][i] - rL[i]*F[j][i-1])*dt/dV[i];
    }

  } else if (dir == JDIR) {
    dz = grid[JDIR].dx;
    inv_dz = grid[JDIR].inv_dx;

    for (l = 0; l<Nlines; l++) {
      i = lines->dom_line_idx[l];
      lidx = lines->lidx[l];
      ridx = lines->ridx[l];
      
      F[lidx][i] = -Hm_B[lidx][i] * (Br_hat[lidx][i] - Br_hat[lidx-1][i])*dz[lidx]*r_1[i]*r_1[i] * 0.5*(Br[lidx][i] + Br[lidx-1][i]);
      for (j=lidx; j<=ridx; j++){
        //[Err] ho aggiunto *r_1[i] nella formula ( e questa modifica sembra ok!)
        F[j][i] = -Hp_B[j][i] * (Br_hat[j+1][i] - Br_hat[j][i])*dz[j]*r_1[i]*r_1[i] * 0.5*(Br[j+1][i] + Br[j][i]);
      }
      
      // Build dU
      for (j=lidx; j<=ridx; j++)
        dUres[j][i] = -(F[j][i] - F[j-1][i])*dt*inv_dz[j];
    }
  }
}

// [Err] Test: decomment this function
/****************************************************************************
* Function to build the bcs of lines
* In the current implementation of this function Data *d is not used
* but I leave it there since before or later it might be needed
*****************************************************************************/
void BoundaryADI_Res(Lines lines[2], const Data *d, Grid *grid, double t, int dir) {
  int i,j,l;
  const double t_sec = t*(UNIT_LENGTH/UNIT_VELOCITY);
  double Bwall;
  double unit_Mfield;
  // [Err]
  // double L = 0.02/UNIT_LENGTH;

  // I compute the wall magnetic field
  unit_Mfield = COMPUTE_UNIT_MFIELD(UNIT_VELOCITY, UNIT_DENSITY);
  curr = current_from_time(t_sec);
  Bwall = BIOTSAV_GAUSS_A_CM(curr, RCAP)/unit_Mfield;

  if (dir == IDIR) {
    /*-----------------------------------------------*/
    /*----  Set bcs for lines in direction IDIR  ----*/
    /*-----------------------------------------------*/
    for (l=0; l<lines[IDIR].N; l++) {
      j = lines[IDIR].dom_line_idx[l];
      /* :::: Axis :::: */
      lines[IDIR].lbound[BDIFF][l].kind = DIRICHLET;
      lines[IDIR].lbound[BDIFF][l].values[0] = 0.0;
      if ( j < j_elec_start) {
        /* :::: Capillary wall (no electrode) :::: */
        lines[IDIR].rbound[BDIFF][l].kind = DIRICHLET;
        lines[IDIR].rbound[BDIFF][l].values[0] = Bwall*rcap_real;
      } else if (j >= j_elec_start && j <= j_cap_inter_end) {
        /* :::: Electrode :::: */
        // [Err] Delete next two lines
        #ifdef ELECTR_B_NEUM
          // [Err] Decomment next lines
          lines[IDIR].rbound[BDIFF][l].kind = NEUMANN_HOM;
          lines[IDIR].rbound[BDIFF][l].values[0] = 0.0;

          // [Err] remove next lines
          // lines[IDIR].rbound[BDIFF][l].kind = DIRICHLET;
          // lines[IDIR].rbound[BDIFF][l].values[0] = Bwall*rcap_real * \
          //    (1 - (grid[JDIR].x_glob[j]-(zcap_real-dzcap_real))/dzcap_real );
        #else
          lines[IDIR].rbound[BDIFF][l].kind = DIRICHLET;
          lines[IDIR].rbound[BDIFF][l].values[0] = Bwall*rcap_real * \
              (1 - (grid[JDIR].x_glob[j]-(zcap_real-dzcap_real))/dzcap_real );
        #endif

        //[Err]
        /* lines[IDIR].rbound[BDIFF][l].kind = DIRICHLET;
            if (grid[JDIR].x_glob[j]>zcap_real-dzcap_real+L) {
              lines[IDIR].rbound[BDIFF][l].values[0] = 0;
            } else {
              lines[IDIR].rbound[BDIFF][l].values[0] = Bwall*rcap_real * \
                (1 - (grid[JDIR].x_glob[j]-(zcap_real-dzcap_real))/L );
            }
        */
        // [Err] end err part
      } else {
        /* :::: Outer domain boundary :::: */
        lines[IDIR].rbound[BDIFF][l].kind = DIRICHLET;
        lines[IDIR].rbound[BDIFF][l].values[0] = 0.0;
      }
    }
  } else if (dir == JDIR) {
    /*-----------------------------------------------*/
    /*----  Set bcs for lines in direction JDIR  ----*/
    /*-----------------------------------------------*/
    for (l=0; l<lines[JDIR].N; l++) {
      i = lines[JDIR].dom_line_idx[l];
      if ( i <= i_cap_inter_end) {
        /* :::: Capillary internal (symmetry plane) :::: */
        lines[JDIR].lbound[BDIFF][l].kind = NEUMANN_HOM;
        lines[JDIR].lbound[BDIFF][l].values[0] = 0.0;
      } else {
        /* :::: Outer capillary wall :::: */
        #ifdef ELECTR_B_NEUM
          lines[JDIR].lbound[BDIFF][l].kind = NEUMANN_HOM;        
          lines[JDIR].lbound[BDIFF][l].values[0] = 0.0;
        #else
          lines[JDIR].lbound[BDIFF][l].kind = DIRICHLET;
          lines[JDIR].lbound[BDIFF][l].values[0] = 0.0;
        #endif
      }
      /* :::: Outer domain boundary :::: */
      lines[JDIR].rbound[BDIFF][l].kind = DIRICHLET;
      lines[JDIR].rbound[BDIFF][l].values[0] = 0.0;
    }
  }
}

/* *************************************************************
Function to complain that the computed eta is non isotropic
* **************************************************************/
void ComplainAnisotropic(double *v, double  *eta,
                         double r, double z, double theta) {
  double unit_Mfield;

  unit_Mfield = COMPUTE_UNIT_MFIELD(UNIT_VELOCITY, UNIT_DENSITY);

  print1("Anisotropic resistivity is not implemented in ADI!");
  print1("\nr,z,theta = %g, %g, %g", r*UNIT_LENGTH,z*UNIT_LENGTH,theta*UNIT_LENGTH);
  print1("\neta = {%g,%g,%g}", eta[1], eta[2], eta[3]);
  print1("\nv[RHO]=%g", v[RHO]*UNIT_DENSITY);
  print1("\nv[iVR]=%g", v[iVR]*UNIT_VELOCITY);
  print1("\nv[iVZ]=%g",v[iVZ]*UNIT_VELOCITY);
  print1("\nv[iVPHI]=%g", v[iVPHI]*UNIT_VELOCITY);
  print1("\nv[PRS]=%g", v[PRS]*UNIT_DENSITY*UNIT_VELOCITY*UNIT_VELOCITY);
  print1("\nv[iBR]=%g", v[iBR]*unit_Mfield);
  print1("\nv[iBZ]=%g", v[iBZ]*unit_Mfield);
  print1("\nv[iBPHI]=%g", v[iBPHI]*unit_Mfield);
}
// [Err] This function is a test, remove it later
// /****************************************************************************
// * Function to build the bcs of lines for a 1D-like problem
// *****************************************************************************/
// void BoundaryADI_Res(Lines lines[2], const Data *d, Grid *grid, double t, int dir) {
//   int i,j,l;
//   const double t_sec = t*(UNIT_LENGTH/UNIT_VELOCITY);
//   double Bwall;
//   double curr, unit_Mfield;
//   // [Err]
//   // double L = 0.02/UNIT_LENGTH;

//   // I compute the wall magnetic field
//   unit_Mfield = COMPUTE_UNIT_MFIELD(UNIT_VELOCITY, UNIT_DENSITY);
//   curr = current_from_time(t_sec);
//   Bwall = BIOTSAV_GAUSS_A_CM(curr, RCAP)/unit_Mfield;

//   if (dir == IDIR) {
//     /*-----------------------------------------------*/
//     /*----  Set bcs for lines in direction IDIR  ----*/
//     /*-----------------------------------------------*/
//     for (l=0; l<lines[IDIR].N; l++) {
//       j = lines[IDIR].dom_line_idx[l];
//       /* :::: Axis :::: */
//       lines[IDIR].lbound[BDIFF][l].kind = DIRICHLET;
//       lines[IDIR].lbound[BDIFF][l].values[0] = 0.0;
//       if ( j < j_elec_start) {
//         /* :::: Capillary wall (no electrode) :::: */
//         lines[IDIR].rbound[BDIFF][l].kind = DIRICHLET;
//         lines[IDIR].rbound[BDIFF][l].values[0] = Bwall*rcap_real;
//       } else if (j >= j_elec_start && j <= j_cap_inter_end) {
//         /* :::: Electrode :::: */
//           lines[IDIR].rbound[BDIFF][l].kind = DIRICHLET;
//           lines[IDIR].rbound[BDIFF][l].values[0] = Bwall*rcap_real;
//       } else {
//         /* :::: Outer domain boundary :::: */
//         // useless in 1D like
//         lines[IDIR].rbound[BDIFF][l].kind = DIRICHLET;
//         lines[IDIR].rbound[BDIFF][l].values[0] = 0.0;
//       }
//     }
//   } else if (dir == JDIR) {
//     /*-----------------------------------------------*/
//     /*----  Set bcs for lines in direction JDIR  ----*/
//     /*-----------------------------------------------*/
//     for (l=0; l<lines[JDIR].N; l++) {
//       i = lines[JDIR].dom_line_idx[l];
//       if ( i <= i_cap_inter_end) {
//         /* :::: Capillary internal (symmetry plane) :::: */
//         lines[JDIR].lbound[BDIFF][l].kind = NEUMANN_HOM;
//         lines[JDIR].lbound[BDIFF][l].values[0] = 0.0;
//       } else {
//         /* :::: Outer capillary wall :::: */
//         // Useless
//           lines[JDIR].lbound[BDIFF][l].kind = DIRICHLET;
//           lines[JDIR].lbound[BDIFF][l].values[0] = 0.0;
//       }
//       /* :::: Outer domain boundary :::: */
//       lines[JDIR].rbound[BDIFF][l].kind = NEUMANN_HOM;
//       lines[JDIR].rbound[BDIFF][l].values[0] = 0.0;
//     }
//   }
// }


#endif