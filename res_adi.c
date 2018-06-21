#include "pluto.h"
#include "adi.h"
#include "capillary_wall.h"
#include "current_table.h"

#define UNUSED(x) (void)(x)

#if RESISTIVITY == ALTERNATING_DIRECTION_IMPLICIT
/****************************************************************************
Function to build the Ip,Im,Jp,Jm for the electrical resistivity problem
(**useless parameter is intentionally unused, to make this function suitable for a pointer
 which also wants that parameter)
*****************************************************************************/
void BuildIJ_Res(const Data *d, Grid *grid, Lines *lines,
                   double **Ip, double **Im, double **Jp,
                   double **Jm, double **CI, double **CJ, double **useless) {
  int i,j,k;
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

  inv_dzi = grid[JDIR].inv_dxi;
  inv_dri = grid[IDIR].inv_dxi;
  inv_dz = grid[JDIR].inv_dx;
  inv_dr = grid[IDIR].inv_dx;
  r_1 = grid[IDIR].r_1;

  theta = grid[KDIR].x;

  r = grid[IDIR].x;
  z = grid[JDIR].x;
  rL = grid[IDIR].xl;
  rR = grid[IDIR].xr;
  zL = grid[JDIR].xl;
  zR = grid[JDIR].xr;

  /*[Opt] This is probably useless, it is here just for debugging purposes*/
  TOT_LOOP(k, j, i) {
    Ip[j][i] = 0.0;
    Im[j][i] = 0.0;
    Jp[j][i] = 0.0;
    Jm[j][i] = 0.0;
    CI[j][i] = 0.0;
    CJ[j][i] = 0.0;
  }
  // [Opt] I could compose the grid-related part once forever (in static variables) and only update eta
  //       in this case I could make that this function allocates Ip Im Jp Jm , instead of letting the rest of
  //       the program do that
  // lines->lidx
  KDOM_LOOP(k) {
    LINES_LOOP(lines[IDIR], l, j, i) {
      /* :::: Ip :::: */
      for (nv=0; nv<NVAR; nv++)
        v[nv] = 0.5 * (Vc[nv][k][j][i] + Vc[nv][k][j][i+1]);
      /*[Rob] I should give also J to Resistive_eta().
        I don't, but I could simply compute it by modifying easily the GetCurrent()
        function, as it only uses J and B element of the Data structure.
        Even easier: I could change the whole adi implementation and make it
        advance directly d->Vc instead of allocating vectors for T and Br */
      Resistive_eta( v, rR[i], z[j], theta[k], NULL, eta);
      if (eta[0] != eta[1] || eta[1] != eta[2]) {
        print1("Anisotropic resistivity is not implemented in ADI!");
        QUIT_PLUTO(1);
      }
      Ip[j][i] = eta[0]*inv_dr[i]*inv_dri[i]/rR[i];
      /* [Opt] Here I could use the already computed eta to compute
      also Im[1,i+1] (since it needs eta at the same interface).
      Doing so I would reduce the calls to Resistive_eta by almost a factor 1/2
      (obviously I could do the same for Im/Ip) */

      /* :::: Im :::: */
      for (nv=0; nv<NVAR; nv++)
        v[nv] = 0.5 * (Vc[nv][k][j][i] + Vc[nv][k][j][i-1]);
      Resistive_eta( v, rL[i], z[j], theta[k], NULL, eta);
      if (eta[0] != eta[1] || eta[1] != eta[2]) {
        print1("Anisotropic resistivity is not implemented in ADI!");
        QUIT_PLUTO(1);
      }
      if (rL[i]!=0.0)
        Im[j][i] = eta[0]*inv_dri[i-1]*inv_dr[i]/rL[i];
      else
        Im[j][i] = eta[0]/(r[i]*r[i])/rR[i];

      /* :::: Jp :::: */
      for (nv=0; nv<NVAR; nv++)
        v[nv] = 0.5 * (Vc[nv][k][j][i] + Vc[nv][k][j+1][i]);
      Resistive_eta( v, r[i], zR[j], theta[k], NULL, eta);
      if (eta[0] != eta[1] || eta[1] != eta[2]) {
        print1("Anisotropic resistivity is not implemented in ADI!");
        QUIT_PLUTO(1);
      }
      Jp[j][i] = eta[0]*inv_dz[j]*inv_dzi[j];

      /* :::: Jm :::: */
      for (nv=0; nv<NVAR; nv++)
        v[nv] = 0.5 * (Vc[nv][k][j][i] + Vc[nv][k][j-1][i]);
      Resistive_eta( v, r[i], zL[j], theta[k], NULL, eta);
      Jm[j][i] = eta[0]*inv_dz[j]*inv_dzi[j-1];

      /* :::: CI :::: */
      CI[j][i] = r_1[i];

      /* :::: CJ :::: */
      CJ[j][i] = 1.0;
    }
  }
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
        F[j][lidx-1] = -Hp_B[j][lidx-1] * (Br[j][lidx] - Br[j][lidx-1])*dr[lidx-1] * 0.5*(Br[j][lidx]*r_1[lidx] + Br[j][lidx-1]*r_1[lidx-1]);
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
      
      for (j=lidx-1; j<=ridx; j++){
        //[Err] ho aggiunto *r_1[i] nella formula ( e questa modifica sembra ok!)
        F[j][i] = -Hp_B[j][i] * (Br[j+1][i] - Br[j][i])*dz[j]*r_1[i]*r_1[i] * 0.5*(Br[j+1][i] + Br[j][i]);
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
  double curr, unit_Mfield;
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
        #ifdef ELECTR_NEUM
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
        #ifdef ELECTR_NEUM
          lines[JDIR].lbound[BDIFF][l].kind = NEUMANN_HOM;        
          lines[JDIR].lbound[BDIFF][l].values[0] = 0.0;
        #else
          lines[JDIR].lbound[BDIFF][l].kind = DIRICHLET;
        #endif
      }
      /* :::: Outer domain boundary :::: */
      lines[JDIR].rbound[BDIFF][l].kind = DIRICHLET;
      lines[JDIR].rbound[BDIFF][l].values[0] = 0.0;
    }
  }
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