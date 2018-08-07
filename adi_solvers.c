/*Functions and utilities for integration of parabolic (currently resistive
and thermal conduction) terms with the Alternating Direction Implicit algorithm*/

// Remarkable comments:
// [Opt] = it can be optimized (in terms of performance)
// [Err] = it is and error (usually introduced on purpose)
// [Rob] = it can/should be made more robust

#include "pluto.h"
#include "adi.h"
#include "capillary_wall.h"
#include "debug_utilities.h"
#include <time.h>
#include <stdlib.h>

/****************************************************************************
Performs an implicit update of a diffusive problem (either for B or for T).
It also applies the bcs on the ghost cells of the output matrix (**v) (useful later
for instance for ResEnergyIncrease())
*****************************************************************************/
void ImplicitUpdate (double **v, double **b, double **source,
                     double **Hp, double **Hm, double **C,
                     Lines *lines, Bcs *lbound, Bcs *rbound,
                     int compute_inflow, double *inflow, Grid *grid,
                     double dt, int dir) {
  /*[Opt] Maybe I could pass to this func. an integer which tells which bc has to be
  used inside the structure *lines, instead of passing separately the bcs (which are
  still also contained inside *lines)*/
  /*[Opt] Maybe I could use g_dir instead of passing dir, but I am afraid of
  doing caos modifiying the value of g_dir for the rest of PLUTO*/
  // int *m, *n, *i, *j; // m and n play the role of i and j (not necessarly respectively)
  // int s, dom_line_idx;
  // const int zero=0;
  int i,j;
  int Nlines = lines->N;
  int ridx, lidx, l;
  /* I allocate these as big as if I had to cover the whole domain, so that I
   don't need to reallocate at every domain line that I update */
  double *diagonal, *upper, *lower, *rhs, *x;
  double *dz, *rR, *rL;

  /*[Opt] Maybe I could do that it allocates static arrays with size NMAX_POINT (=max(NX1_TOT,NX2_TOT)) ?*/
  if (dir == IDIR) {
    diagonal = ARRAY_1D(NX1_TOT, double);
    rhs = ARRAY_1D(NX1_TOT, double);
    upper = ARRAY_1D(NX1_TOT, double);
    lower = ARRAY_1D(NX1_TOT, double);
    x = ARRAY_1D(NX1_TOT, double);
  } else if (dir == JDIR) {
    diagonal = ARRAY_1D(NX2_TOT, double);
    rhs = ARRAY_1D(NX2_TOT, double);
    upper = ARRAY_1D(NX2_TOT, double);
    lower = ARRAY_1D(NX2_TOT, double);
    x = ARRAY_1D(NX2_TOT, double);
  }
    /*[Opt] potrei sempificare il programma facendo che
  i membri di destra delle assegnazione prendono degli indici
  che puntano a j o a i a seconda del valore di dir*/
  // if (dir == IDIR) {
  //   j = &dom_line_idx;
  //   m = &zero;
  //   n = &s;
  //   i = &lidx;
  // } else if (dir == JDIR) {
  //   j = &lidx;
  //   m = &s;
  //   n = &zero;
  //   i = &dom_line_idx;
  // }
  // for (l = 0; l<Nlines; l++) {
  //   dom_line_idx = lines->dom_line_idx;
  //   lidx = lines->lidx;
  //   ridx = lines->ridx;
  //   s = lidx;
  //     [*j+*m][*i+*n];
  // }
  if (dir == IDIR) {
  /********************
  * Case direction IDIR
  *********************/
    dz = grid[JDIR].dx_glob;

    for (l = 0; l < Nlines; l++) {
      j = lines->dom_line_idx[l];
      lidx = lines->lidx[l];
      ridx = lines->ridx[l];

      upper[lidx] = -dt/C[j][lidx]*Hp[j][lidx];
      lower[ridx] = -dt/C[j][ridx]*Hm[j][ridx];
      rhs[lidx] = b[j][lidx];
      rhs[ridx] = b[j][ridx];
      for (i = lidx+1; i < ridx; i++) {
        diagonal[i] = 1 + dt/C[j][i] * (Hp[j][i]+Hm[j][i]);
        rhs[i] = b[j][i];
        upper[i] = -dt/C[j][i]*Hp[j][i];
        lower[i] = -dt/C[j][i]*Hm[j][i];
      }
      /* I include the effect of the source */
      if (source != NULL) {
        for (i = lidx; i <= ridx; i++)
          rhs[i] += source[j][i]*dt;
      }
      // I set the Bcs for left boundary
      if (lbound[l].kind == DIRICHLET){
        diagonal[lidx] = 1 + dt/C[j][lidx]*(Hp[j][lidx]+2*Hm[j][lidx]);
        rhs[lidx] += dt/C[j][lidx]*Hm[j][lidx]*2*lbound[l].values[0];
      } else if (lbound[l].kind == NEUMANN_HOM) {
        diagonal[lidx] = 1 + dt/C[j][lidx]*Hp[j][lidx];
      } else {
        print1("\n[ImplicitUpdate]Error setting left bc (in dir i), not known bc kind!");
        QUIT_PLUTO(1);
      }
      // I set the Bcs for right boundary
      if (rbound[l].kind == DIRICHLET){
        diagonal[ridx] = 1 + dt/C[j][ridx]*(2*Hp[j][ridx]+Hm[j][ridx]);
        rhs[ridx] += dt/C[j][ridx]*Hp[j][ridx]*2*rbound[l].values[0];
      } else if (rbound[l].kind == NEUMANN_HOM) {
        diagonal[ridx] = 1 + dt/C[j][ridx]*Hm[j][ridx];
      } else {
        print1("\n[ImplicitUpdate]Error setting right bc (in dir i), not known bc kind!");
        QUIT_PLUTO(1);
      }

      /*---------------------------------------------------------------------*/
      /* --- Now I solve the system --- */
      tdm_solver( x+lidx, diagonal+lidx, upper+lidx, lower+lidx+1, rhs+lidx, ridx-lidx+1);
      /*[Opt] Is this for a waste of time? maybe I could engineer better the use of tdm_solver function(or the way it is written)*/
      for (i=lidx; i<=ridx; i++)
        v[j][i] = x[i];
      
      /*---------------------------------------------------------------------*/
      /*--- I set the boundary values (ghost cells) in the solution
            [I do it now as I for the NEUMANN conditions I could't do it before solving the tridiag. system] ---*/
      // Cells near left boundary
      if (lbound[l].kind == DIRICHLET){
        // I assign the ghost value (needed by ResEnergyIncrease and maybe others..)
        v[j][lidx-1] = 2*lbound[l].values[0] - b[j][lidx];
        if (compute_inflow) {
          /*--- I compute the inflow ---*/
          *inflow += (v[j][lidx-1]-v[j][lidx]) * Hm[j][lidx] * CONST_PI*dz[j] * dt;
        }

      } else if (lbound[l].kind == NEUMANN_HOM) {
        /* I assign the ghost value (needed by ResEnergyIncrease and maybe others..).*/
        v[j][lidx-1] = v[j][lidx];
        if (compute_inflow) {
          /*--- I compute the inflow (0!!!)---*/
          *inflow += 0;
        }

      } else {
        print1("\n[ImplcitUpdate]Error setting left ghost in solution (in dir i), not known bc kind!");
        QUIT_PLUTO(1);
      }

      // Cells near right boundary
      if (rbound[l].kind == DIRICHLET){
        v[j][ridx+1] = 2*rbound[l].values[0] - v[j][ridx];
        if (compute_inflow) {
          /*--- I compute the inflow ---*/
          *inflow += (v[j][ridx+1]-v[j][ridx]) * Hp[j][ridx] * 2*CONST_PI*dz[j] * dt;
        }

      } else if (rbound[l].kind == NEUMANN_HOM) {
        v[j][ridx+1] = v[j][ridx];
        if (compute_inflow) {
          /*--- I compute the inflow (0!!!)---*/
          *inflow += 0;
        } 

      } else {
        print1("\n[ImplcitUpdate]Error setting right ghost in solution (in dir i), not known bc kind!");
        QUIT_PLUTO(1);
      }
    }
  } else if (dir == JDIR) {
    /********************
    * Case direction JDIR
    *********************/
    rR = grid[IDIR].xr_glob;
    rL = grid[IDIR].xl_glob;

    for (l = 0; l < Nlines; l++) {
      i = lines->dom_line_idx[l];
      lidx = lines->lidx[l];
      ridx = lines->ridx[l];

      upper[lidx] = -dt/C[lidx][i]*Hp[lidx][i];
      lower[ridx] = -dt/C[ridx][i]*Hm[ridx][i];
      rhs[lidx] = b[lidx][i];
      rhs[ridx] = b[ridx][i];
      for (j = lidx+1; j < ridx; j++) {
        diagonal[j] = 1 + dt/C[j][i] * (Hp[j][i]+Hm[j][i]);
        rhs[j] = b[j][i];
        upper[j] = -dt/C[j][i]*Hp[j][i];
        lower[j] = -dt/C[j][i]*Hm[j][i];
      }
      /* I include the effect of the source */
      if (source != NULL) {
        for (j = lidx; j <= ridx; j++)
          rhs[j] += source[j][i]*dt;
      }
      // I set the Bcs for left boundary
      if (lbound[l].kind == DIRICHLET){
        diagonal[lidx] = 1 + dt/C[lidx][i]*(Hp[lidx][i]+2*Hm[lidx][i]);
        rhs[lidx] += dt/C[lidx][i]*Hm[lidx][i]*2*lbound[l].values[0];
      } else if (lbound[l].kind == NEUMANN_HOM) {
        diagonal[lidx] = 1 + dt/C[lidx][i]*Hp[lidx][i];
      } else {
        print1("\n[ImplicitUpdate]Error setting left bc (in dir j), not known bc kind!");
        QUIT_PLUTO(1);
      }
      // I set the Bcs for right boundary
      if (rbound[l].kind == DIRICHLET){
        diagonal[ridx] = 1 + dt/C[ridx][i]*(2*Hp[ridx][i]+Hm[ridx][i]);
        rhs[ridx] += dt/C[ridx][i]*Hp[ridx][i]*2*rbound[l].values[0];
      } else if (rbound[l].kind == NEUMANN_HOM) {
        diagonal[ridx] = 1 + dt/C[ridx][i]*Hm[ridx][i];
      } else {
        print1("\n[ImplicitUpdate]Error setting right bcs (in dir j), not known bc kind!");
        QUIT_PLUTO(1);
      }

      /*---------------------------------------------------------------------*/
      /* --- Now I solve the system --- */
      tdm_solver( x+lidx, diagonal+lidx, upper+lidx, lower+lidx+1, rhs+lidx, ridx-lidx+1);
      for (j=lidx; j<=ridx; j++)
        v[j][i] = x[j];

      /*---------------------------------------------------------------------*/
      /*--- I set the boundary values (ghost cells) in the solution
            [I do it now as I for the NEUMANN conditions I could't do it before solving the tridiag. system] ---*/
            /*--- I set the boundary values (ghost cells) ---*/
      // Cells near left boundary
      if (lbound[l].kind == DIRICHLET){
        v[lidx-1][i] = 2*lbound[l].values[0] - v[lidx][i];
        if (compute_inflow) {
          /*--- I compute the inflow ---*/
          *inflow += (v[lidx-1][i]-v[lidx][i]) * Hm[lidx][i] * CONST_PI*(rR[i]*rR[i]-rL[i]*rL[i]) * dt;
        }

      } else if (lbound[l].kind == NEUMANN_HOM) {
        v[lidx-1][i] = v[lidx][i];
        if (compute_inflow) {
          /*--- I compute the inflow (0!!!)---*/
          *inflow += 0;
        }

      } else {
        print1("\n[ImplcitUpdate]Error setting left ghost in solution (in dir j), not known bc kind!");
        QUIT_PLUTO(1);
      }

      // Cells near right boundary
      if (rbound[l].kind == DIRICHLET) {
        v[ridx+1][i] = 2*rbound[l].values[0] - v[ridx][i];
        if (compute_inflow) {
          /*--- I compute the inflow ---*/
          *inflow += (v[ridx+1][i]-v[ridx][i]) * Hp[ridx][i] * CONST_PI*(rR[i]*rR[i]-rL[i]*rL[i]) * dt;
        }
      } else if (rbound[l].kind == NEUMANN_HOM) {
        v[ridx+1][i] = v[ridx][i];
        if (compute_inflow) {
          /*--- I compute the inflow (0!!!)---*/
          *inflow += 0;
        }
      } else {
        print1("\n[ImplcitUpdate]Error setting right ghost in solution (in dir j), not known bc kind!");
        QUIT_PLUTO(1);
      }
    }
  } else {
    print1("[ImplicitUpdate] Unimplemented choice for 'dir'!");
    QUIT_PLUTO(1);
  }

  FreeArray1D(diagonal);
  FreeArray1D(x);
  FreeArray1D(upper);
  FreeArray1D(lower);
  FreeArray1D(rhs);
}
//
/****************************************************************************
Performs an explicit update of a diffusive problem (either for B or for T).
It also applies the bcs on the ghost cells of the input matrix (**b) (useful later
for instance for ResEnergyIncrease())
*****************************************************************************/
void ExplicitUpdate (double **v, double **b, double **source,
                     double **Hp, double **Hm, double **C,
                     Lines *lines, Bcs *lbound, Bcs *rbound, double dt,
                    //  int compute_inflow, double *inflow, Grid *grid,
                     int dir) {
  int i,j,l;
  int ridx, lidx;
  int Nlines = lines->N;
  static double **rhs;
  static int first_call = 1;

  if (first_call) {
    rhs = ARRAY_2D(NX2_TOT, NX1_TOT, double);
    first_call = 0;
  }

  if (dir == IDIR) {
    /********************
    * Case direction IDIR
    *********************/
    for (l = 0; l < Nlines; l++) {
      j = lines->dom_line_idx[l];
      lidx = lines->lidx[l];
      ridx = lines->ridx[l];

      if (source != NULL) {
        for (i = lidx; i <= ridx; i++)
          rhs[j][i] = b[j][i] + source[j][i]*dt;
      } else {
        /*[Opt] Maybe I could assign directly the address. BE CAREFUL:
        if I assign the address, then I have to recover the old address of rhs (by saving temporarly the old rhs address
        inside another variable), otherwise
        at the next call of this function I will write over the memory of the old b*/
        for (i = lidx; i <= ridx; i++)
          rhs[j][i] = b[j][i];
      }

      /*--- I set the boundary values (ghost cells) ---*/
      // Cells near left boundary
      if (lbound[l].kind == DIRICHLET){
        // I assign the ghost value (needed by ResEnergyIncrease and maybe others..)
        // [Err] decomment next line
        b[j][lidx-1] = 2*lbound[l].values[0] - b[j][lidx];
        //[Err] Experimental: 2nd order accurate bc for cell centered FD (as explained in L.Chen draft on FDM)
        // b[j][lidx-1] = 1/3*b[j][lidx+1] + 8/3*lbound[l].values[0] - 2*b[j][lidx];
      } else if (lbound[l].kind == NEUMANN_HOM) {
        /* I assign the ghost value (needed by ResEnergyIncrease and maybe others..).*/
        b[j][lidx-1] = b[j][lidx];
      } else {
        print1("\n[ExplicitUpdate]Error setting left bc (in dir i), not known bc kind!");
        QUIT_PLUTO(1);
      }
      // Cells near right boundary
      if (rbound[l].kind == DIRICHLET){
        // [Err] decomment next line
        b[j][ridx+1] = 2*rbound[l].values[0] - b[j][ridx];
        //[Err] Experimental: 2nd order accurate bc for cell centered FD (as explained in L.Chen draft on FDM)
        // b[j][ridx+1] = 1/3*b[j][ridx-1] + 8/3*rbound[l].values[0] - 2*b[j][ridx];
      } else if (rbound[l].kind == NEUMANN_HOM) {
        b[j][ridx+1] = b[j][ridx];
      } else {
        print1("\n[ExplicitUpdate]Error setting right bc (in dir i), not known bc kind!");
        QUIT_PLUTO(1);
      }

      /*--- Actual update ---*/
      for (i = lidx; i <= ridx; i++)
        v[j][i] = rhs[j][i] + dt/C[j][i] * (b[j][i+1]*Hp[j][i] - b[j][i]*(Hp[j][i]+Hm[j][i]) + b[j][i-1]*Hm[j][i]);
    }
  } else if (dir == JDIR) {
    /********************
    * Case direction JDIR
    *********************/
    for (l = 0; l < Nlines; l++) {
      i = lines->dom_line_idx[l];
      lidx = lines->lidx[l];
      ridx = lines->ridx[l];

      if (source != NULL) {
        for (j = lidx; j <= ridx; j++)
          rhs[j][i] = b[j][i] + source[j][i]*dt;
      } else {
        /*[Opt] Maybe I could assign directly the address. BE CAREFUL:
        if I assign the address, then I have to recover the old address of rhs (by saving temporarly the old rhs address
        inside another variable), otherwise
        at the next call of this function I will write over the memory of the old b*/
        for (j = lidx; j <= ridx; j++)
          rhs[j][i] = b[j][i];
      }

      /*--- I set the boundary values (ghost cells) ---*/
      // Cells near left boundary
      if (lbound[l].kind == DIRICHLET){
        // [Err] decomment next line
        b[lidx-1][i] = 2*lbound[l].values[0] - b[lidx][i];
        //[Err] Experimental: 2nd order accurate bc for cell centered FD (as explained in L.Chen draft on FDM)
        // b[lidx-1][i] = 1/3*b[lidx+1][i] + 8/3*lbound[l].values[0] - 2*b[lidx][i];
      } else if (lbound[l].kind == NEUMANN_HOM) {
        b[lidx-1][i] = b[lidx][i];
      } else {
        print1("\n[ExplicitUpdate]Error setting left bc (in dir j), not known bc kind!");
        QUIT_PLUTO(1);
      }
      // Cells near right boundary
      if (rbound[l].kind == DIRICHLET){
        // [Err] decomment next line
        b[ridx+1][i] = 2*rbound[l].values[0] - b[ridx][i];
        //[Err] Experimental: 2nd order accurate bc for cell centered FD (as explained in L.Chen draft on FDM)
        // b[ridx+1][i] = 1/3*b[ridx-1][i] + 8/3*rbound[l].values[0] - 2*b[ridx][i];
      } else if (rbound[l].kind == NEUMANN_HOM) {
        b[ridx+1][i] = b[ridx][i];
      } else {
        print1("\n[ExplicitUpdate]Error setting right bc (in dir j), not known bc kind!");
        QUIT_PLUTO(1);
      }

      /*--- Actual update ---*/
      for (j = lidx; j <= ridx; j++){
        v[j][i] = rhs[j][i] + dt/C[j][i] * (b[j+1][i]*Hp[j][i] - b[j][i]*(Hp[j][i]+Hm[j][i]) + b[j-1][i]*Hm[j][i]);
        // print1("v[%d][%d]=%e\n", j,i,v[j][i]);
      }
    }
  } else {
    print1("[ImplicitUpdate] Unimplemented choice for 'dir'!");
    QUIT_PLUTO(1);
  }
}

/************************************************************
 * Solve a linear system made by a tridiagonal matrix.
 * See  https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
 * for info (accessed on 24/3/2017)
 *
 * N: the size of the arrays x, diagonal, right_hand_side.
 * x: solution
 * *********************************************************/
void tdm_solver(double *x, double const *diagonal, double const *upper,
                double const *lower, double const *right_hand_side, int const N) {
  /*[Opt] Is it really needed to define two new arrays?
          Maybe I can make that it uses directly the diagonal, upper,
          lower arrays, modifying them, the only problem is that I
          don't like the idea as a principle, that the function
          modifies its actual input with no apparent reason*/
  double up[N-1];
  double rhs[N];
  int i;

  /*[Opt] Maybe there is another way to copy the arrays... more clever, shorter, faster..*/
  for (i=0;i<N;i++)
    rhs[i] = right_hand_side[i];
  for (i=0;i<N-1;i++)
    up[i] = upper[i];

  up[0] = upper[0]/diagonal[0];
  for (i=1; i<N-1; i++)
    up[i] = up[i] / (diagonal[i]-lower[i-1]*up[i-1]);

  rhs[0] = rhs[0]/diagonal[0];
  for (i=1; i<N; i++)
    rhs[i] = (rhs[i] - lower[i-1]*rhs[i-1]) / (diagonal[i] - lower[i-1]*up[i-1]);

  x[N-1] = rhs[N-1];
  for (i=N-2; i>-1; i--)
    x[i] = rhs[i] - up[i]*x[i+1];
}

/* ***********************************************************
 * Peachman-Rachford ADI method
 * 
 * input: diff = BDIFF or TDIFF
 *        int order = FIRST_IDIR or FIRST_JDIR: tells whether the order of the directions
 *                     must be IDIR, JDIR (FIRST_IDIR) or JDIR, IDIR (FIRST_JDIR).
 *        **dEdT: may point to NULL in case diff == BDIFF
 *        **dUres: will not be updated if diff != BDIFF
 *                 (it is my opinion that for the sake of clarity it is better not to
 *                  "merge" in one single variable the quantities **dEdT and **dUres,
 *                  even though I never need both of them at the same time)
 * ***********************************************************/
void PeacemanRachford(double **v_new, double **v_old,
                      double **dUres, double **dEdT,
                      const Data *d, Grid *grid,
                      Lines *lines, int diff, int order,
                      double dt, double t0) {

    static double **v_aux; // auxiliary solution vector
    static double **Ip, **Im, **CI, **Jp, **Jm, **CJ;
    static int first_call = 1;
    double **H1p, **H1m, **H2p, **H2m, **C1, **C2;
    // void (*BoundaryADI) (Lines, const Data, Grid, double);
    BoundaryADI *ApplyBCs;
    BuildIJ *MakeIJ;
    int dir1, dir2;
    #if (JOULE_EFFECT_AND_MAG_ENG)
      static double **dUres_aux; // auxiliary vector containing a contribution to ohmic heating
      int l,i,j;
    #endif
    #if (JOULE_EFFECT_AND_MAG_ENG && !POW_INSIDE_ADI)
      static double **Br_avg, **dUres_aux1;
    #endif

    if (first_call) {
      v_aux = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      #if (JOULE_EFFECT_AND_MAG_ENG)
        dUres_aux = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      #endif
      #if (JOULE_EFFECT_AND_MAG_ENG && !POW_INSIDE_ADI)
        Br_avg = ARRAY_2D(NX2_TOT, NX1_TOT, double);
        dUres_aux1 = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      #endif
      Ip = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      Im = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      Jp = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      Jm = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      CI = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      CJ = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      first_call = 0;
    }

    /* Set the direction order*/
    if (order == FIRST_IDIR) {
      H1p = Ip;     H1m = Im;
      H2p = Jp;     H2m = Jm;
      C1 = CI;      C2 = CJ;
      dir1 = IDIR;  dir2 = JDIR;
    } else if (order == FIRST_JDIR) {
      H1p = Jp;     H1m = Jm;
      H2p = Ip;     H2m = Im;
      C1 = CJ;      C2 = CI;
      dir1 = JDIR;  dir2 = IDIR;
    }

    switch (diff) {
      #if RESISTIVITY==ALTERNATING_DIRECTION_IMPLICIT
        case BDIFF:
          ApplyBCs = BoundaryADI_Res;
          MakeIJ = BuildIJ_Res;
          break;
      #endif
      #if THERMAL_CONDUCTION==ALTERNATING_DIRECTION_IMPLICIT
        case TDIFF:
          ApplyBCs = BoundaryADI_TC;
          MakeIJ = BuildIJ_TC;
          break;
      #endif
      default:
        print1("\n[PeachmanRachford]Wrong setting for diffusion (diff) problem");
        QUIT_PLUTO(1);
        break;
    }

    ApplyBCs(lines, d, grid, t0, dir1);
    MakeIJ(d, grid, lines, Ip, Im, Jp, Jm, CI, CJ, dEdT);
    
    /**********************************
     (a.1) Explicit update sweeping DIR1
    **********************************/
    ExplicitUpdate (v_aux, v_old, NULL, H1p, H1m, C1, &lines[dir1],
                    lines[dir1].lbound[diff], lines[dir1].rbound[diff], 0.5*dt, dir1);
    #if (JOULE_EFFECT_AND_MAG_ENG && POW_INSIDE_ADI)
      if (diff == BDIFF) {
        // [Err] Decomment next line
        // [Opt] You could modify and make that the ResEnergyEncrease automatically updates a Ures variable,
        //       instead of doing it a line later 
        ResEnergyIncrease(dUres_aux, H1p, H1m, v_old, grid, &lines[dir1],
                          EN_CONS_CHECK, &en_res_in,
                          0.5*dt, dir1);
        LINES_LOOP(lines[IDIR], l, j, i)
          dUres[j][i] = dUres_aux[j][i];
      }
    #endif

    /**********************************
     (a.2) Implicit update sweeping DIR2
    **********************************/
    ApplyBCs(lines, d, grid, t0 + dt*0.5, dir2);
    ImplicitUpdate (v_new, v_aux, NULL, H2p, H2m, C2, &lines[dir2],
                      lines[dir2].lbound[diff], lines[dir2].rbound[diff],
                      (diff == TDIFF) && EN_CONS_CHECK, &en_tc_in, grid,
                      0.5*dt, dir2);
    #if (JOULE_EFFECT_AND_MAG_ENG && POW_INSIDE_ADI)
      if (diff == BDIFF) {
        // [Err] Decomment next line
        ResEnergyIncrease(dUres_aux, H2p, H2m, v_new, grid, &lines[dir2],
                          EN_CONS_CHECK, &en_res_in,
                          0.5*dt, dir2);
        LINES_LOOP(lines[IDIR], l, j, i)
          dUres[j][i] += dUres_aux[j][i];
      }
    #endif

    /**********************************
     (b.1) Explicit update sweeping DIR2
    **********************************/
    ExplicitUpdate (v_aux, v_new, NULL, H2p, H2m, C2, &lines[dir2],
                    lines[dir2].lbound[diff], lines[dir2].rbound[diff], 0.5*dt, dir2);
    #if (JOULE_EFFECT_AND_MAG_ENG && POW_INSIDE_ADI)
      if (diff == BDIFF) {
        /* [Opt]: I could inglobate this call to ResEnergyIncrease in the previous one by using dt_res_reduced instead of 0.5*dt_res_reduced
           (but in this way it is more readable)*/
        // [Err] Decomment next line       
        ResEnergyIncrease(dUres_aux, H2p, H2m, v_new, grid, &lines[dir2],
                          EN_CONS_CHECK, &en_res_in,
                          0.5*dt, dir2);
        LINES_LOOP(lines[IDIR], l, j, i)
          dUres[j][i] += dUres_aux[j][i];
      }    
    #endif

    /**********************************
     (b.2) Implicit update sweeping DIR1
    **********************************/
    ApplyBCs(lines, d, grid, t0 + dt, dir1);
    ImplicitUpdate (v_new, v_aux, NULL, H1p, H1m, C1, &lines[dir1],
                      lines[dir1].lbound[diff], lines[dir1].rbound[diff],
                      (diff == TDIFF) && EN_CONS_CHECK, &en_tc_in, grid,
                      0.5*dt, dir1);
    #if (JOULE_EFFECT_AND_MAG_ENG && POW_INSIDE_ADI)
      if (diff == BDIFF) {
        // [Err] Decomment next line
        ResEnergyIncrease(dUres_aux, H1p, H1m, v_new, grid, &lines[dir1],
                          EN_CONS_CHECK, &en_res_in,
                          0.5*dt, dir1);
        LINES_LOOP(lines[IDIR], l, j, i)
          dUres[j][i] += dUres_aux[j][i];
      }
    #endif
    #if (JOULE_EFFECT_AND_MAG_ENG && !POW_INSIDE_ADI)
      if (diff == BDIFF) {
        // I average Br
        LINES_LOOP_EXTENDED(lines[IDIR], l, j, i) {
          Br_avg[j][i] = sqrt((v_new[j][i]*v_new[j][i] + v_new[j][i]*v_old[j][i] + v_old[j][i]*v_old[j][i])/3);
        }
        // I compute the fluxes of poynting vector
        ResEnergyIncrease(dUres_aux, H2p, H2m, Br_avg, grid, &lines[dir2],
                          EN_CONS_CHECK, &en_res_in,
                          dt, dir2);
        ResEnergyIncrease(dUres_aux1, H1p, H1m, Br_avg, grid, &lines[dir1],
                          EN_CONS_CHECK, &en_res_in,
                          dt, dir1);
        // I update the increase in energy
        LINES_LOOP(lines[IDIR], l, j, i) {
          dUres[j][i] = dUres_aux[j][i];
          dUres[j][i] += dUres_aux1[j][i];
        }
      }
    #endif
}

/* ***********************************************************
 * Modified Peachman-Rachford ADI method (I have no clue whether this
 * is docuemnted in literature and how accurate it is. I hope it is fine
 * actually I use it since I have seen that in improves
 * the stability properties of the diffusion of the magnetic field
 * and computation of dUres).
 * 
 * input: diff = BDIFF or TDIFF
 *        int order = FIRST_IDIR or FIRST_JDIR: tells whether the order of the directions
 *                     must be IDIR, JDIR (FIRST_IDIR) or JDIR, IDIR (FIRST_JDIR).
 *        **dEdT: may point to NULL in case diff == BDIFF
 *        **dUres: will not be updated if diff != BDIFF
 *                 (it is my opinion that for the sake of clarity it is better not to
 *                  "merge" in one single variable the quantities **dEdT and **dUres,
 *                  even though I never need both of them at the same time)
 * ***********************************************************/
void PeacemanRachfordMod(double **v_new, double **v_old,
                      double **dUres, double **dEdT,
                      const Data *d, Grid *grid,
                      Lines *lines, int diff, int order,
                      double dt, double t0, double fract) {

    static double **v_aux; // auxiliary solution vector
    static double **Ip, **Im, **CI, **Jp, **Jm, **CJ;
    static int first_call = 1;
    double **H1p, **H1m, **H2p, **H2m, **C1, **C2;
    // void (*BoundaryADI) (Lines, const Data, Grid, double);
    BoundaryADI *ApplyBCs;
    BuildIJ *MakeIJ;
    int dir1, dir2;
    #if (JOULE_EFFECT_AND_MAG_ENG)
      static double **dUres_aux; // auxiliary vector containing a contribution to ohmic heating
      int l,i,j;
    #endif
    #if (JOULE_EFFECT_AND_MAG_ENG && !POW_INSIDE_ADI)
      static double **Br_avg, **dUres_aux1;
    #endif

    if (first_call) {
      v_aux = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      #if (JOULE_EFFECT_AND_MAG_ENG)
        dUres_aux = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      #endif
      #if (JOULE_EFFECT_AND_MAG_ENG && !POW_INSIDE_ADI)
        Br_avg = ARRAY_2D(NX2_TOT, NX1_TOT, double);
        dUres_aux1 = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      #endif
      Ip = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      Im = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      Jp = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      Jm = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      CI = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      CJ = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      first_call = 0;
    }

    /* Set the direction order*/
    if (order == FIRST_IDIR) {
      H1p = Ip;     H1m = Im;
      H2p = Jp;     H2m = Jm;
      C1 = CI;      C2 = CJ;
      dir1 = IDIR;  dir2 = JDIR;
    } else if (order == FIRST_JDIR) {
      H1p = Jp;     H1m = Jm;
      H2p = Ip;     H2m = Im;
      C1 = CJ;      C2 = CI;
      dir1 = JDIR;  dir2 = IDIR;
    }

    switch (diff) {
      #if RESISTIVITY==ALTERNATING_DIRECTION_IMPLICIT
        case BDIFF:
          ApplyBCs = BoundaryADI_Res;
          MakeIJ = BuildIJ_Res;
          break;
      #endif
      #if THERMAL_CONDUCTION==ALTERNATING_DIRECTION_IMPLICIT
        case TDIFF:
          ApplyBCs = BoundaryADI_TC;
          MakeIJ = BuildIJ_TC;
          break;
      #endif
      default:
        print1("\n[PeachmanRachford]Wrong setting for diffusion (diff) problem");
        QUIT_PLUTO(1);
        break;
    }

    ApplyBCs(lines, d, grid, t0, dir1);
    MakeIJ(d, grid, lines, Ip, Im, Jp, Jm, CI, CJ, dEdT);
    
    /**********************************
     (a.1) Explicit update sweeping DIR1
    **********************************/
    ExplicitUpdate (v_aux, v_old, NULL, H1p, H1m, C1, &lines[dir1],
                    lines[dir1].lbound[diff], lines[dir1].rbound[diff], fract*dt, dir1);
    #if (JOULE_EFFECT_AND_MAG_ENG && POW_INSIDE_ADI)
      if (diff == BDIFF) {
        // [Err] Decomment next line
        // [Opt] You could modify and make that the ResEnergyEncrease automatically updates a Ures variable,
        //       instead of doing it a line later 
        ResEnergyIncrease(dUres_aux, H1p, H1m, v_old, grid, &lines[dir1],
                          EN_CONS_CHECK, &en_res_in,
                          fract*dt, dir1);
        LINES_LOOP(lines[IDIR], l, j, i)
          dUres[j][i] = dUres_aux[j][i];
      }
    #endif
    #ifdef DEBUG_EMA
      printf("\nafter expl dir1:\n"); 
      printf("\nv_old\n");
      printmat(v_old, NX2_TOT, NX1_TOT);
      printf("\nv_aux\n");
      printmat(v_aux, NX2_TOT, NX1_TOT);
      #if (JOULE_EFFECT_AND_MAG_ENG && POW_INSIDE_ADI)
        printf("\ndUres_aux\n");
        printmat(dUres_aux, NX2_TOT, NX1_TOT);
      #endif
    #endif

    /**********************************
     (a.2) Implicit update sweeping DIR2
    **********************************/
    ApplyBCs(lines, d, grid, t0 + dt*(1-fract), dir2);
    ImplicitUpdate (v_new, v_aux, NULL, H2p, H2m, C2, &lines[dir2],
                      lines[dir2].lbound[diff], lines[dir2].rbound[diff],
                      (diff == TDIFF) && EN_CONS_CHECK, &en_tc_in, grid,
                      (1-fract)*dt, dir2);
    #if (JOULE_EFFECT_AND_MAG_ENG && POW_INSIDE_ADI)
      if (diff == BDIFF) {
        // [Err] Decomment next line
        ResEnergyIncrease(dUres_aux, H2p, H2m, v_new, grid, &lines[dir2],
                          EN_CONS_CHECK, &en_res_in,
                          (1-fract)*dt, dir2);
        LINES_LOOP(lines[IDIR], l, j, i)
          dUres[j][i] += dUres_aux[j][i];
      }
    #endif
    #ifdef DEBUG_EMA
      printf("\nafter impl dir2:\n");
      printf("\nv_new\n");
      printmat(v_new, NX2_TOT, NX1_TOT);
      #if (JOULE_EFFECT_AND_MAG_ENG && POW_INSIDE_ADI)
        printf("\ndUres_aux\n");
        printmat(dUres_aux, NX2_TOT, NX1_TOT);
      #endif
    #endif

    /**********************************
     (b.1) Explicit update sweeping DIR2
    **********************************/
    ExplicitUpdate (v_aux, v_new, NULL, H2p, H2m, C2, &lines[dir2],
                    lines[dir2].lbound[diff], lines[dir2].rbound[diff], fract*dt, dir2);
    #if (JOULE_EFFECT_AND_MAG_ENG && POW_INSIDE_ADI)
      if (diff == BDIFF) {
        /* [Opt]: I could inglobate this call to ResEnergyIncrease in the previous one by using dt_res_reduced instead of 0.5*dt_res_reduced
           (but in this way it is more readable)*/
        // [Err] Decomment next line       
        ResEnergyIncrease(dUres_aux, H2p, H2m, v_new, grid, &lines[dir2],
                          EN_CONS_CHECK, &en_res_in,
                          fract*dt, dir2);
        LINES_LOOP(lines[IDIR], l, j, i)
          dUres[j][i] += dUres_aux[j][i];
      }    
    #endif
    #ifdef DEBUG_EMA
      printf("\nafter expl dir2:\n");
      printf("\nv_new\n");
      printmat(v_new, NX2_TOT, NX1_TOT);
      printf("\nv_aux\n");
      printmat(v_aux, NX2_TOT, NX1_TOT);
      #if (JOULE_EFFECT_AND_MAG_ENG && POW_INSIDE_ADI)
        printf("\ndUres_aux\n");
        printmat(dUres_aux, NX2_TOT, NX1_TOT);
      #endif
    #endif

    /**********************************
     (b.2) Implicit update sweeping DIR1
    **********************************/
    ApplyBCs(lines, d, grid, t0 + dt, dir1);
    ImplicitUpdate (v_new, v_aux, NULL, H1p, H1m, C1, &lines[dir1],
                      lines[dir1].lbound[diff], lines[dir1].rbound[diff],
                      (diff == TDIFF) && EN_CONS_CHECK, &en_tc_in, grid,
                      (1-fract)*dt, dir1);
    #if (JOULE_EFFECT_AND_MAG_ENG && POW_INSIDE_ADI)
      if (diff == BDIFF) {
        // [Err] Decomment next line
        ResEnergyIncrease(dUres_aux, H1p, H1m, v_new, grid, &lines[dir1],
                          EN_CONS_CHECK, &en_res_in,
                          (1-fract)*dt, dir1);
        LINES_LOOP(lines[IDIR], l, j, i)
          dUres[j][i] += dUres_aux[j][i];
      }
    #endif
    #if (JOULE_EFFECT_AND_MAG_ENG && !POW_INSIDE_ADI)
      if (diff == BDIFF) {
        // I average Br
        LINES_LOOP_EXTENDED(lines[IDIR], l, j, i) {
          Br_avg[j][i] = sqrt((v_new[j][i]*v_new[j][i] + v_new[j][i]*v_old[j][i] + v_old[j][i]*v_old[j][i])/3);
        }
        // I compute the fluxes of poynting vector
        ResEnergyIncrease(dUres_aux, H2p, H2m, Br_avg, grid, &lines[dir2],
                          EN_CONS_CHECK, &en_res_in,
                          dt, dir2);
        ResEnergyIncrease(dUres_aux1, H1p, H1m, Br_avg, grid, &lines[dir1],
                          EN_CONS_CHECK, &en_res_in,
                          dt, dir1);
        // I update the increase in energy
        LINES_LOOP(lines[IDIR], l, j, i) {
          dUres[j][i] = dUres_aux[j][i];
          dUres[j][i] += dUres_aux1[j][i];
        }
      }
    #endif
    #ifdef DEBUG_EMA
      printf("\nafter impl dir1:\n");
      printf("\nv_new\n");
      printmat(v_new, NX2_TOT, NX1_TOT);
      #if (JOULE_EFFECT_AND_MAG_ENG && POW_INSIDE_ADI)
        printf("\ndUres_aux\n");
        printmat(dUres_aux, NX2_TOT, NX1_TOT);
      #endif
    #endif
}

/* ***********************************************************
 * Douglas-Rachford ADI method, absolutely NON OPTIMIZED FOR SPEED!!!
 * 
 * input: diff = BDIFF or TDIFF
 *        int order = FIRST_IDIR or FIRST_JDIR: tells whether the order of the directions
 *                     must be IDIR, JDIR (FIRST_IDIR) or JDIR, IDIR (FIRST_JDIR).
 *        **dEdT: may point to NULL in case diff == BDIFF
 *        **dUres: will not be updated if diff != BDIFF
 *                 (it is my opinion that for the sake of clarity it is better not to
 *                  "merge" in one single variable the quantities **dEdT and **dUres,
 *                  even though I never need both of them at the same time)
 * ***********************************************************/
void DouglasRachford(double **v_new, double **v_old,
                      double **dUres, double **dEdT,
                      const Data *d, Grid *grid,
                      Lines *lines, int diff, int order,
                      double dt, double t0) {

    static double **v_aux, **v_hat; // auxiliary solution vector
    static double **Ip, **Im, **CI, **Jp, **Jm, **CJ;
    static int first_call = 1;
    double **H1p, **H1m, **H2p, **H2m, **C1, **C2;
    // void (*BoundaryADI) (Lines, const Data, Grid, double);
    BoundaryADI *ApplyBCs;
    BuildIJ *MakeIJ;
    int dir1, dir2;
    int l,i,j;
    #if (JOULE_EFFECT_AND_MAG_ENG)
      static double **dUres_aux; // auxiliary vector containing a contribution to ohmic heating
    #endif
    #if (JOULE_EFFECT_AND_MAG_ENG && !POW_INSIDE_ADI)
      static double **Br_avg, **dUres_aux1;
    #endif

    print1("\nAttenzione al calcolo dell'energia che entra dai bordi per conduzione/elettromagnetica:\n");
    print1("\npotrebbe essere che sia sbagliata per come ho implmentato lo schema D-R (e per l'uso di variabili globali)\n");

    if (first_call) {
      v_aux = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      v_hat = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      #if (JOULE_EFFECT_AND_MAG_ENG)
        dUres_aux = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      #endif
      #if (JOULE_EFFECT_AND_MAG_ENG && !POW_INSIDE_ADI)
        Br_avg = ARRAY_2D(NX2_TOT, NX1_TOT, double);
        dUres_aux1 = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      #endif
      Ip = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      Im = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      Jp = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      Jm = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      CI = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      CJ = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      first_call = 0;
    }

    /* Set the direction order*/
    if (order == FIRST_IDIR) {
      H1p = Ip;     H1m = Im;
      H2p = Jp;     H2m = Jm;
      C1 = CI;      C2 = CJ;
      dir1 = IDIR;  dir2 = JDIR;
    } else if (order == FIRST_JDIR) {
      H1p = Jp;     H1m = Jm;
      H2p = Ip;     H2m = Im;
      C1 = CJ;      C2 = CI;
      dir1 = JDIR;  dir2 = IDIR;
    }

    switch (diff) {
      #if RESISTIVITY==ALTERNATING_DIRECTION_IMPLICIT
        case BDIFF:
          ApplyBCs = BoundaryADI_Res;
          MakeIJ = BuildIJ_Res;
          break;
      #endif
      #if THERMAL_CONDUCTION==ALTERNATING_DIRECTION_IMPLICIT
        case TDIFF:
          ApplyBCs = BoundaryADI_TC;
          MakeIJ = BuildIJ_TC;
          break;
      #endif
      default:
        print1("\n[PeachmanRachford]Wrong setting for diffusion (diff) problem");
        QUIT_PLUTO(1);
        break;
    }

    ApplyBCs(lines, d, grid, t0, dir1);
    MakeIJ(d, grid, lines, Ip, Im, Jp, Jm, CI, CJ, dEdT);
    
    /**********************************
     (a.1) Explicit update sweeping DIR1
    **********************************/
    ExplicitUpdate (v_aux, v_old, NULL, H1p, H1m, C1, &lines[dir1],
                    lines[dir1].lbound[diff], lines[dir1].rbound[diff], dt, dir1);
    /**********************************
     (a.2) Implicit update sweeping DIR2
    **********************************/
    ApplyBCs(lines, d, grid, t0 + dt, dir2);
    // I compute phi^ (and save it in v_hat)
    ImplicitUpdate (v_hat, v_aux, NULL, H2p, H2m, C2, &lines[dir2],
                      lines[dir2].lbound[diff], lines[dir2].rbound[diff],
                      (diff == TDIFF) && EN_CONS_CHECK, &en_tc_in, grid,
                      dt, dir2);
    /**********************************
     (b.1) Explicit update sweeping DIR2
    **********************************/
    // I compute phi~ (and save it in v_aux)
    ExplicitUpdate (v_aux, v_hat, NULL, H2p, H2m, C2, &lines[dir2],
                    lines[dir2].lbound[diff], lines[dir2].rbound[diff], dt, dir2);
    /**********************************
     (b.1_bis) This step is in order not to rewrite the updating routines,
      it is not convenient from a performance point of view
    **********************************/
    // I compute phi** (and save it in v_aux)
    LINES_LOOP(lines[IDIR], l, j, i)
      v_aux[j][i] += v_old[j][i] - v_hat[j][i];
    /**********************************
     (b.2) Implicit update sweeping DIR1
    **********************************/
    ApplyBCs(lines, d, grid, t0 + dt, dir1);
    ImplicitUpdate (v_new, v_aux, NULL, H1p, H1m, C1, &lines[dir1],
                    lines[dir1].lbound[diff], lines[dir1].rbound[diff],
                    (diff == TDIFF) && EN_CONS_CHECK, &en_tc_in, grid,
                    dt, dir1);
    //[Opt] Questa è una porcheria, avanzo esplicitamente per dt=0 solo per dare le bc in dir2 a v_new
    ExplicitUpdate (v_new, v_new, NULL, H2p, H2m, C2, &lines[dir2],
                    lines[dir2].lbound[diff], lines[dir2].rbound[diff], 0.0, dir2);
    #if (JOULE_EFFECT_AND_MAG_ENG && POW_INSIDE_ADI)
      if (diff == BDIFF) {
        // For one advancement of the energy I don't use DouglasRachf variant, since I don't need it!
        ResEnergyIncrease(dUres_aux, H1p, H1m, v_new, grid, &lines[dir1],
                          EN_CONS_CHECK, &en_res_in,
                          dt, dir1);
        LINES_LOOP(lines[IDIR], l, j, i)
          dUres[j][i] = dUres_aux[j][i];
        ResEnergyIncrease_DouglasRachford(dUres_aux, H2p, H2m, v_new, v_hat, grid, &lines[dir2], dt, dir2);
        LINES_LOOP(lines[IDIR], l, j, i)
          dUres[j][i] += dUres_aux[j][i];
      }
    #endif
    #if (JOULE_EFFECT_AND_MAG_ENG && !POW_INSIDE_ADI)
      #error "!POW_INSIDE_ADI" is not implemented in DouglasRachford
    #endif
}

/* ***********************************************************
 * Fractional-Theta (ADI) method
 * (see Glowinski, R.: Splitting methods for the numerical solution of the incompressible Navier-Stokes equations, 1985)
 * 
 * input: diff = BDIFF or TDIFF
 *        int order = FIRST_IDIR or FIRST_JDIR: tells whether the order of the directions
 *                     must be IDIR, JDIR (FIRST_IDIR) or JDIR, IDIR (FIRST_JDIR).
 *        **dEdT: may point to NULL in case diff == BDIFF
 *        **dUres: will not be updated if diff != BDIFF
 *                 (it is my opinion that for the sake of clarity it is better not to
 *                  "merge" in one single variable the quantities **dEdT and **dUres,
 *                  even though I never need both of them at the same time)
 * ***********************************************************/
void FractionalTheta(double **v_new, double **v_old,
                     double **dUres, double **dEdT,
                     const Data *d, Grid *grid,
                     Lines *lines, int diff, int order,
                     double dt, double t0, double theta) {
    
    static double **v_aux; // auxiliary solution vector
    static double **Ip, **Im, **CI, **Jp, **Jm, **CJ;
    static int first_call = 1;
    double **H1p, **H1m, **H2p, **H2m, **C1, **C2;
    // void (*BoundaryADI) (Lines, const Data, Grid, double);
    BoundaryADI *ApplyBCs;
    BuildIJ *MakeIJ;
    int dir1, dir2;
    #if (JOULE_EFFECT_AND_MAG_ENG && POW_INSIDE_ADI)
      static double **dUres_aux; // auxiliary vector containing a contribution to ohmic heating
      int l,i,j;
    #endif

    // Check theta correctness
    if (theta<0.0 || theta > 0.5){
      print1("theta=%d, outisde range ]0,0.5[");
      QUIT_PLUTO(1);
    }

    if (first_call) {
      v_aux = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      #if (JOULE_EFFECT_AND_MAG_ENG && POW_INSIDE_ADI)
        dUres_aux = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      #endif
      Ip = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      Im = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      Jp = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      Jm = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      CI = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      CJ = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      first_call = 0;
    }

    /* Set the direction order*/
    if (order == FIRST_IDIR) {
      H1p = Ip;     H1m = Im;
      H2p = Jp;     H2m = Jm;
      C1 = CI;      C2 = CJ;
      dir1 = IDIR;  dir2 = JDIR;
    } else if (order == FIRST_JDIR) {
      H1p = Jp;     H1m = Jm;
      H2p = Ip;     H2m = Im;
      C1 = CJ;      C2 = CI;
      dir1 = JDIR;  dir2 = IDIR;
    }

    switch (diff) {
      #if RESISTIVITY==ALTERNATING_DIRECTION_IMPLICIT
        case BDIFF:
          ApplyBCs = BoundaryADI_Res;
          MakeIJ = BuildIJ_Res;
          break;
      #endif
      #if THERMAL_CONDUCTION==ALTERNATING_DIRECTION_IMPLICIT
        case TDIFF:
          ApplyBCs = BoundaryADI_TC;
          MakeIJ = BuildIJ_TC;
          break;
      #endif
      default:
        print1("\n[FractionalTheta]Wrong setting for diffusion (diff) problem");
        QUIT_PLUTO(1);
        break;
    }

    ApplyBCs(lines, d, grid, t0, dir1);
    MakeIJ(d, grid, lines, Ip, Im, Jp, Jm, CI, CJ, dEdT);
    
    /**********************************
     (a.1) Explicit update sweeping DIR1
    **********************************/
    ExplicitUpdate (v_aux, v_old, NULL, H1p, H1m, C1, &lines[dir1],
                    lines[dir1].lbound[diff], lines[dir1].rbound[diff], theta*dt, dir1);
    #if (JOULE_EFFECT_AND_MAG_ENG && POW_INSIDE_ADI)
      if (diff == BDIFF) {
        // [Err] Decomment next line
        // [Opt] You could modify and make that the ResEnergyEncrease automatically updates a Ures variable,
        //       instead of doing it a line later 
        ResEnergyIncrease(dUres_aux, H1p, H1m, v_old, grid, &lines[dir1],
                          EN_CONS_CHECK, &en_res_in,
                          theta*dt, dir1);
        LINES_LOOP(lines[IDIR], l, j, i)
          dUres[j][i] = dUres_aux[j][i];
      }
    #endif

    /**********************************
     (a.2) Implicit update sweeping DIR2
    **********************************/
    ApplyBCs(lines, d, grid, t0 + theta*dt, dir2);
    ImplicitUpdate (v_new, v_aux, NULL, H2p, H2m, C2, &lines[dir2],
                      lines[dir2].lbound[diff], lines[dir2].rbound[diff],
                      (diff == TDIFF) && EN_CONS_CHECK, &en_tc_in, grid,
                      theta*dt, dir2);
    #if (JOULE_EFFECT_AND_MAG_ENG && POW_INSIDE_ADI)
      if (diff == BDIFF) {
        // [Err] Decomment next line
        ResEnergyIncrease(dUres_aux, H2p, H2m, v_new, grid, &lines[dir2],
                          EN_CONS_CHECK, &en_res_in,
                          theta*dt, dir2);
        LINES_LOOP(lines[IDIR], l, j, i)
          dUres[j][i] += dUres_aux[j][i];
      }
    #endif

    /**********************************
     (b.1) Explicit update sweeping DIR2
    **********************************/
    ExplicitUpdate (v_aux, v_new, NULL, H2p, H2m, C2, &lines[dir2],
                    lines[dir2].lbound[diff], lines[dir2].rbound[diff], (1-2*theta)*dt, dir2);
    #if (JOULE_EFFECT_AND_MAG_ENG && POW_INSIDE_ADI)
      if (diff == BDIFF) {
        /* [Opt]: I could inglobate this call to ResEnergyIncrease in the previous one by using dt_res_reduced instead of 0.5*dt_res_reduced
           (but in this way it is more readable)*/
        // [Err] Decomment next line       
        ResEnergyIncrease(dUres_aux, H2p, H2m, v_new, grid, &lines[dir2],
                          EN_CONS_CHECK, &en_res_in,
                          (1-2*theta)*dt, dir2);
        LINES_LOOP(lines[IDIR], l, j, i)
          dUres[j][i] += dUres_aux[j][i];
      }    
    #endif

    /**********************************
     (b.2) Implicit update sweeping DIR1
    **********************************/
    ApplyBCs(lines, d, grid, t0 + (1-theta)*dt, dir1);
    ImplicitUpdate (v_new, v_aux, NULL, H1p, H1m, C1, &lines[dir1],
                    lines[dir1].lbound[diff], lines[dir1].rbound[diff],
                    (diff == TDIFF) && EN_CONS_CHECK, &en_tc_in, grid,
                    (1-2*theta)*dt, dir1);
    #if (JOULE_EFFECT_AND_MAG_ENG && POW_INSIDE_ADI)
      if (diff == BDIFF) {
        // [Err] Decomment next line
        ResEnergyIncrease(dUres_aux, H1p, H1m, v_new, grid, &lines[dir1],
                          EN_CONS_CHECK, &en_res_in,
                          (1-2*theta)*dt, dir1);
        LINES_LOOP(lines[IDIR], l, j, i)
          dUres[j][i] += dUres_aux[j][i];
      }
    #endif

    /**********************************
     (c.1) Explicit update sweeping DIR1
    **********************************/
    ExplicitUpdate (v_aux, v_new, NULL, H1p, H1m, C1, &lines[dir1],
                    lines[dir1].lbound[diff], lines[dir1].rbound[diff], theta*dt, dir1);
    #if (JOULE_EFFECT_AND_MAG_ENG && POW_INSIDE_ADI)
      if (diff == BDIFF) {
        // [Err] Decomment next line
        // [Opt] You could modify and make that the ResEnergyEncrease automatically updates a Ures variable,
        //       instead of doing it a line later 
        ResEnergyIncrease(dUres_aux, H1p, H1m, v_new, grid, &lines[dir1],
                          EN_CONS_CHECK, &en_res_in,
                          theta*dt, dir1);
        LINES_LOOP(lines[IDIR], l, j, i)
          dUres[j][i] += dUres_aux[j][i];
      }
    #endif

    /**********************************
     (c.2) Implicit update sweeping DIR2
    **********************************/
    ApplyBCs(lines, d, grid, t0 + dt, dir2);
    ImplicitUpdate (v_new, v_aux, NULL, H2p, H2m, C2, &lines[dir2],
                      lines[dir2].lbound[diff], lines[dir2].rbound[diff],
                      (diff == TDIFF) && EN_CONS_CHECK, &en_tc_in, grid,
                      theta*dt, dir2);
    #if (JOULE_EFFECT_AND_MAG_ENG && POW_INSIDE_ADI)
      if (diff == BDIFF) {
        // [Err] Decomment next line
        ResEnergyIncrease(dUres_aux, H2p, H2m, v_new, grid, &lines[dir2],
                          EN_CONS_CHECK, &en_res_in,
                          theta*dt, dir2);
        LINES_LOOP(lines[IDIR], l, j, i)
          dUres[j][i] += dUres_aux[j][i];
      }
    #endif
}

/* ***********************************************************
 * SplitImplicit method.
 * 
 * input: diff = BDIFF or TDIFF
 *        int order = FIRST_IDIR or FIRST_JDIR: tells whether the order of the directions
 *                     must be IDIR, JDIR (FIRST_IDIR) or JDIR, IDIR (FIRST_JDIR).
 *        **dEdT: may point to NULL in case diff == BDIFF
 *        **dUres: will not be updated if diff != BDIFF
 *                 (it is my opinion that for the sake of clarity it is better not to
 *                  "merge" in one single variable the quantities **dEdT and **dUres,
 *                  even though I never need both of them at the same time)
 * ***********************************************************/
void SplitImplicit(double **v_new, double **v_old,
                      double **dUres, double **dEdT,
                      const Data *d, Grid *grid,
                      Lines *lines, int diff, int order,
                      double dt, double t0) {

    static double **v_aux; // auxiliary solution vector
    static double **Ip, **Im, **CI, **Jp, **Jm, **CJ;
    static int first_call = 1;
    double **H1p, **H1m, **H2p, **H2m, **C1, **C2;
    // void (*BoundaryADI) (Lines, const Data, Grid, double);
    BoundaryADI *ApplyBCs;
    BuildIJ *MakeIJ;
    int dir1, dir2;
    #if (JOULE_EFFECT_AND_MAG_ENG)
      static double **dUres_aux; // auxiliary vector containing a contribution to ohmic heating
      int l,i,j;
    #endif

    if (first_call) {
      v_aux = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      #if (JOULE_EFFECT_AND_MAG_ENG)
        dUres_aux = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      #endif

      Ip = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      Im = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      Jp = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      Jm = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      CI = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      CJ = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      first_call = 0;
    }

    /* Set the direction order*/
    if (order == FIRST_IDIR) {
      H1p = Ip;     H1m = Im;
      H2p = Jp;     H2m = Jm;
      C1 = CI;      C2 = CJ;
      dir1 = IDIR;  dir2 = JDIR;
    } else if (order == FIRST_JDIR) {
      H1p = Jp;     H1m = Jm;
      H2p = Ip;     H2m = Im;
      C1 = CJ;      C2 = CI;
      dir1 = JDIR;  dir2 = IDIR;
    }

    switch (diff) {
      #if RESISTIVITY==ALTERNATING_DIRECTION_IMPLICIT
        case BDIFF:
          ApplyBCs = BoundaryADI_Res;
          MakeIJ = BuildIJ_Res;
          break;
      #endif
      #if THERMAL_CONDUCTION==ALTERNATING_DIRECTION_IMPLICIT
        case TDIFF:
          ApplyBCs = BoundaryADI_TC;
          MakeIJ = BuildIJ_TC;
          break;
      #endif
      default:
        print1("\n[SplitImplicit]Wrong setting for diffusion (diff) problem");
        QUIT_PLUTO(1);
        break;
    }

    ApplyBCs(lines, d, grid, t0+dt, dir1);
    MakeIJ(d, grid, lines, Ip, Im, Jp, Jm, CI, CJ, dEdT);

    /**********************************
     (a) Implicit update sweeping DIR1
    **********************************/
    ImplicitUpdate (v_aux, v_old, NULL, H1p, H1m, C1, &lines[dir1],
                      lines[dir1].lbound[diff], lines[dir1].rbound[diff],
                      (diff == TDIFF) && EN_CONS_CHECK, &en_tc_in, grid,
                      dt, dir1);
    #if (JOULE_EFFECT_AND_MAG_ENG && POW_INSIDE_ADI)
      if (diff == BDIFF) {
        // [Err] Decomment next line
        ResEnergyIncrease(dUres_aux, H1p, H1m, v_aux, grid, &lines[dir1],
                          EN_CONS_CHECK, &en_res_in,
                          dt, dir1);
        LINES_LOOP(lines[IDIR], l, j, i)
          dUres[j][i] = dUres_aux[j][i];
      }
    #endif
    #ifdef DEBUG_EMA
      printf("\nafter impl dir1:\n");
      printf("\nv_new\n");
      printmat(v_new, NX2_TOT, NX1_TOT);
      #if (JOULE_EFFECT_AND_MAG_ENG && POW_INSIDE_ADI)
        printf("\ndUres_aux\n");
        printmat(dUres_aux, NX2_TOT, NX1_TOT);
      #endif
    #endif

    /**********************************
     (b) Implicit update sweeping DIR2
    **********************************/
    ApplyBCs(lines, d, grid, t0 + dt, dir2);
    ImplicitUpdate (v_new, v_aux, NULL, H2p, H2m, C2, &lines[dir2],
                      lines[dir2].lbound[diff], lines[dir2].rbound[diff],
                      (diff == TDIFF) && EN_CONS_CHECK, &en_tc_in, grid,
                      dt, dir2);
    #if (JOULE_EFFECT_AND_MAG_ENG && POW_INSIDE_ADI)
      if (diff == BDIFF) {
        // [Err] Decomment next line
        ResEnergyIncrease(dUres_aux, H2p, H2m, v_new, grid, &lines[dir2],
                          EN_CONS_CHECK, &en_res_in,
                          dt, dir2);
        LINES_LOOP(lines[IDIR], l, j, i)
          dUres[j][i] += dUres_aux[j][i];
      }
    #endif
    #ifdef DEBUG_EMA
      printf("\nafter impl dir1:\n");
      printf("\nv_new\n");
      printmat(v_new, NX2_TOT, NX1_TOT);
      #if (JOULE_EFFECT_AND_MAG_ENG && POW_INSIDE_ADI)
        printf("\ndUres_aux\n");
        printmat(dUres_aux, NX2_TOT, NX1_TOT);
      #endif
    #endif
}

/*******************************************************
 * COSE DA FARE, ma che sono secondarie:
 *
 * 1) Mettere un po' di cicli #if di controllo che la geometria
 *    il modello (MHD) e altro siano ok.
 * 2) Pensare alla compatabilità con STS o EXPL
 * 3) Capire se come è scritto va bene anche per griglia stretchata
 *    (è chiaro che comunque se stretcho la griglia lentamente va bene,
 *    ma rigorosamente parlando, va bene? Forse devo ragionare sullo sviluppo di taylor
 *    per trovarmi le derivate delle incognite sui punti di griglia)
 * 5) Alberto dice di provare dopo eventualemente (se vedo probelmi o se voglio migliorare accuratezzax) a far aggiornare a t+dt/2 Jmp,Imp
 * 6) fare che viene stampate nell'output di pluto anche il dt parabolico che ci sarebbe con step esplicito
 * 7) Attenzione alla griglia stretchata: forse devo cambiare la
 *    discretizzazione delle derivate se voglio usare rigorosamente una griglia stretchata
 * 8) Forse è meglio fare che ExplicitUpdate, ImplicitUpdate, ResEnergyIncrease usano il valore di bordo
 *    contenuto nell bc punto e basta, senza calcolare il valore che avrebbe una cella di ghost situata oltre l'ultima cella fisica.
 *    Le condizioni al contorno(di dirichlet) dicono il valore dell'incognita al bordo, quindi andrebbe sistemato anche il valore
 *    di Ip,Im,Jp,Jm al bordo (perchè non si differenzia per dr o dz ma per rR[ridx]-r[ridx] (o r[lidx]-rL[lidx]))
 ********************************************************/
