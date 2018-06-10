/*Functions and utilities for integration of parabolic (currently resistive
and thermal conduction) terms with the Alternating Directino Implicit algorithm*/
#include "pluto.h"
#include "adi.h"
#include "capillary_wall.h"

#if KBEG != KEND
  #error grid in k direction should only be of 1 point
#endif

#ifndef REVERSE_ADI_DIRS
  #define DIR1    IDIR
  #define DIR2    JDIR
  #define H1p_B   Ip_B
  #define H1m_B   Im_B
  #define H2p_B   Jp_B
  #define H2m_B   Jm_B
  #define H1p_T   Ip_T
  #define H1m_T   Im_T
  #define H2p_T   Jp_T
  #define H2m_T   Jm_T
  #define C1_B    CI_B
  #define C2_B    CJ_B
  #define C1_T    CI_T
  #define C2_T    CJ_T

#else
  #define DIR1    JDIR
  #define DIR2    IDIR
  #define H1p_B   Jp_B
  #define H1m_B   Jm_B
  #define H2p_B   Ip_B
  #define H2m_B   Im_B
  #define H1p_T   Jp_T
  #define H1m_T   Jm_T
  #define H2p_T   Ip_T
  #define H2m_T   Im_T
  #define C1_B    CJ_B
  #define C2_B    CI_B
  #define C1_T    CJ_T
  #define C2_T    CI_T 

#endif

// Remarkable comments:
// [Opt] = it can be optimized (in terms of performance)
// [Err] = it is and error (usually introduced on purpose)
// [Rob] = it can/should be made more robust

void ADI(const Data *d, Time_Step *Dts, Grid *grid) {
  static int first_call=1;
  int i,j,k, l;
  
  // [Err] Are you sure you want to do multiple steps inside a single adi call?
  int s;
  double t_start_sub_res, t_start_sub_tc;
  
  static Lines lines[2]; /*I define two of them as they are 1 per direction (r and z)*/
  #if RESISTIVITY == ALTERNATING_DIRECTION_IMPLICIT
    static double **Ip_B, **Im_B, **Jp_B, **Jm_B, **CI_B, **CJ_B;
    static double **Bra1, **Bra2, **Brb1, **Brb2;
    // static double **Br_avg;
    static double **Br;
    // Energy increse(due to electro-magnetics) terms
    static double **dUres_a1, **dUres_a2, **dUres_b1, **dUres_b2;
    int const adi_res_steps = NSUBS_RES_ADI;
    double dt_res_reduced;
  #endif
  #if THERMAL_CONDUCTION == ALTERNATING_DIRECTION_IMPLICIT
    static double **Ip_T, **Im_T, **Jp_T, **Jm_T, **CI_T, **CJ_T;
    static double **Ta1, **Ta2, **Tb1, **Tb2;
    static double **T;
    static double **dEdT;
    double v[NVAR]; /*[Ema] I hope that NVAR as dimension is fine!*/
    // double rhoe_old, rhoe_new;
    int nv;
    int const adi_tc_steps = NSUBS_TC_ADI;
    double dt_tc_reduced;
  #endif
  const double dt = g_dt;
  double ****Uc, ****Vc;
  double *r, *r_1;

  // Find the remarkable indexes (if they had not been found before)
  if (capillary_not_set) {
    if (SetRemarkableIdxs(grid)){
      print1("\nError while setting remarkable points!");
      QUIT_PLUTO(1);
    }
  }
  /* Some shortcuts */
  Vc = d->Vc;
  Uc = d->Uc;
  r = grid[IDIR].x_glob;
  r_1 = grid[IDIR].r_1;

  /* -------------------------------------------------------------------
  Build geometry and allocate some stuff
  ----------------------------------------------------------------------*/
  if (first_call) {
    GeometryADI(lines, grid);
    #if RESISTIVITY == ALTERNATING_DIRECTION_IMPLICIT
      Ip_B = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      Im_B = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      Jp_B = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      Jm_B = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      CI_B = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      CJ_B = ARRAY_2D(NX2_TOT, NX1_TOT, double);

      Bra1 = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      Bra2 = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      Brb1 = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      Brb2 = ARRAY_2D(NX2_TOT, NX1_TOT, double);

      Br = ARRAY_2D(NX2_TOT, NX1_TOT, double);

      // Br_avg = ARRAY_2D(NX2_TOT, NX1_TOT, double);

      dUres_a1 = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      dUres_a2 = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      dUres_b1 = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      dUres_b2 = ARRAY_2D(NX2_TOT, NX1_TOT, double);
    #endif

    #if THERMAL_CONDUCTION == ALTERNATING_DIRECTION_IMPLICIT
      Ip_T = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      Im_T = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      Jp_T = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      Jm_T = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      CI_T = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      CJ_T = ARRAY_2D(NX2_TOT, NX1_TOT, double);

      Ta1 = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      Ta2 = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      Tb1 = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      Tb2 = ARRAY_2D(NX2_TOT, NX1_TOT, double);

      T = ARRAY_2D(NX2_TOT, NX1_TOT, double);

      dEdT = ARRAY_2D(NX2_TOT, NX1_TOT, double);

      TOT_LOOP (k,j,i) {
        Ta1[j][i] = 0.0;
        Ta2[j][i] = 0.0;
        Tb1[j][i] = 0.0;
        Tb2[j][i] = 0.0;
      }
    #endif

    first_call=0;
  }

  t_start_sub_res = g_time; /*g_time è: "The current integration time."(dalla docuementazione in Doxigen)*/
  t_start_sub_tc = g_time; /*g_time è: "The current integration time."(dalla docuementazione in Doxigen)*/
  #if RESISTIVITY == ALTERNATING_DIRECTION_IMPLICIT
    dt_res_reduced = dt/adi_res_steps;
    //[Opt] Move this call inside  BuildIJ_X()
    BoundaryRes_ADI(lines, d, grid, t_start_sub_res);
  #endif
  #if THERMAL_CONDUCTION == ALTERNATING_DIRECTION_IMPLICIT
    dt_tc_reduced = dt/adi_tc_steps;
    //[Opt] Move this call inside  BuildIJ_X()
    BoundaryTC_ADI(lines, d, grid, t_start_sub_tc);
  #endif

  /* -------------------------------------------------------------------------
      Compute the conservative vector in order to start the cycle.
      This step will be useless if the data structure
      contains Vc as well as Uc (for future improvements).
      [Ema] (Comment copied from sts.c)
    --------------------------------------------------------------------------- */
  PrimToConsLines (Vc, Uc, lines);

  //[Opt] Move this call inside  BuildIJ_X()
  Boundary(d, ALL_DIR, grid);

  /* -------------------------------------------------------------------------
        Build the temperature (T) and/or the product magnetic field with radius (B*r)
      ------------------------------------------------------------------------- */

  #if THERMAL_CONDUCTION == ALTERNATING_DIRECTION_IMPLICIT
    #if EOS==IDEAL
      DOM_LOOP(k,j,i) T[j][i] = Vc[PRS][k][j][i]/Vc[RHO][k][j][i];
    #elif EOS==PVTE_LAW
      DOM_LOOP(k,j,i) {
        for (nv=NVAR; nv--;) v[nv] = Vc[nv][k][j][i];
        if (GetPV_Temperature(v, &(T[j][i]) )!=0) {
          print1("ADI:[Ema] Error computing temperature!\n");
        }
        Tb2[j][i] = T[j][i] = T[j][i] / KELVIN;
      }
    #else
      print1("ADI:[Ema] Error computing temperature, this EOS not implemented!")
    #endif
    BuildIJ_TC(d, grid, lines, Ip_T, Im_T, Jp_T, Jm_T, CI_T, CJ_T, dEdT);
  #endif
  #if RESISTIVITY == ALTERNATING_DIRECTION_IMPLICIT
    // Build a handy magnetic field matrix
    DOM_LOOP(k,j,i) {
      Brb2[j][i] = Br[j][i] = r[i]*Uc[k][j][i][BX3];
    }
    BuildIJ_Res(d, grid, lines, Ip_B, Im_B, Jp_B, Jm_B, CI_B, CJ_B);
  #endif

  #if THERMAL_CONDUCTION == ALTERNATING_DIRECTION_IMPLICIT
    for (s=0; s<adi_tc_steps; s++) {
   /* --------------------------------------------------------
      Peachman-Rachford Algorithm
     -------------------------------------------------------- */
      
      /**********************************
       (a.1) Explicit update sweeping DIR1
      **********************************/
      ExplicitUpdate (Ta1, Tb2, NULL, H1p_T, H1m_T, C1_T, &lines[DIR1],
                        lines[DIR1].lbound[TDIFF], lines[DIR1].rbound[TDIFF], 0.5*dt_tc_reduced, DIR1);
      /**********************************
       (a.2) Implicit update sweeping DIR2
      **********************************/
      BoundaryTC_ADI(lines, d, grid, t_start_sub_tc+0.5*dt_tc_reduced); // Get bcs at half step (not exaclty at t+0.5*dt_reduced)
      ImplicitUpdate (Ta2, Ta1, NULL, H2p_T, H2m_T, C2_T, &lines[DIR2],
                      lines[DIR2].lbound[TDIFF], lines[DIR2].rbound[TDIFF], 0.5*dt_tc_reduced, DIR2);
      /**********************************
        (b.1) Explicit update sweeping DIR2
      **********************************/
      ExplicitUpdate (Tb1, Ta2, NULL, H2p_T, H2m_T, C2_T, &lines[DIR2],
                      lines[DIR2].lbound[TDIFF], lines[DIR2].rbound[TDIFF], 0.5*dt_tc_reduced, DIR2);
      /**********************************
        (b.2) Implicit update sweeping DIR1
      **********************************/
      BoundaryTC_ADI(lines, d, grid, t_start_sub_tc+dt_tc_reduced);
      ImplicitUpdate (Tb2, Tb1, NULL, H1p_T, H1m_T, C1_T, &lines[DIR1],
                    lines[DIR1].lbound[TDIFF], lines[DIR1].rbound[TDIFF], 0.5*dt_tc_reduced, DIR1);
      t_start_sub_tc += dt_tc_reduced;
    }
  #endif

  #if RESISTIVITY == ALTERNATING_DIRECTION_IMPLICIT
  for (s=0; s<adi_res_steps; s++) {    
  /* --------------------------------------------------------
      Peachman-Rachford Algorithm
     -------------------------------------------------------- */
    
    /**********************************
     (a.1) Explicit update sweeping DIR1
    **********************************/
    // [Err] Remove next #if lines
    // #if (HAVE_ENERGY && JOULE_EFFECT_AND_MAG_ENG)
      // ResEnergyIncrease(dUres_a1, H1p_B, H1m_B, Br, grid, &lines[DIR1], 0.5*dt_res_reduced, DIR1);
      // ResEnergyIncrease(dUres_a2, H2p_B, H2m_B, Br, grid, &lines[DIR2], 0.5*dt_res_reduced, DIR2);
    // #endif
    ExplicitUpdate (Bra1, Brb2, NULL, H1p_B, H1m_B, C1_B, &lines[DIR1],
                    lines[DIR1].lbound[BDIFF], lines[DIR1].rbound[BDIFF], 0.5*dt_res_reduced, DIR1);
    #if (HAVE_ENERGY && JOULE_EFFECT_AND_MAG_ENG)
    // [Err] Decomment next line    
      ResEnergyIncrease(dUres_a1, H1p_B, H1m_B, Brb2, grid, &lines[DIR1], 0.5*dt_res_reduced, DIR1);
    #endif

    /**********************************
     (a.2) Implicit update sweeping DIR2
    **********************************/
    BoundaryRes_ADI(lines, d, grid, t_start_sub_res+0.5*dt_res_reduced);
    ImplicitUpdate (Bra2, Bra1, NULL, H2p_B, H2m_B, C2_B, &lines[DIR2],
                      lines[DIR2].lbound[BDIFF], lines[DIR2].rbound[BDIFF], 0.5*dt_res_reduced, DIR2);
    #if (HAVE_ENERGY && JOULE_EFFECT_AND_MAG_ENG)
    // [Err] Decomment next line    
      ResEnergyIncrease(dUres_a2, H2p_B, H2m_B, Bra2, grid, &lines[DIR2], 0.5*dt_res_reduced, DIR2);
    #endif
    // [Err] Remove next four lines
    // ResEnergyIncrease(dUres_a1, Ip_B, Im_B, Bra2, grid, &lines[IDIR], 0.5*dt_res_reduced, IDIR);
    // ResEnergyIncrease(dUres_a2, Jp_B, Jm_B, Bra2, grid, &lines[DIR2], 0.5*dt_res_reduced, DIR2);
    // ResEnergyIncrease(dUres_b1, Ip_B, Im_B, Bra2, grid, &lines[IDIR], 0.5*dt_res_reduced, IDIR);
    // ResEnergyIncrease(dUres_b2, Jp_B, Jm_B, Bra2, grid, &lines[DIR2], 0.5*dt_res_reduced, DIR2);

    /**********************************
     (b.1) Explicit update sweeping DIR2
    **********************************/
    ExplicitUpdate (Brb1, Bra2, NULL, H2p_B, H2m_B, C2_B, &lines[DIR2],
                    lines[DIR2].lbound[BDIFF], lines[DIR2].rbound[BDIFF], 0.5*dt_res_reduced, DIR2);
    #if (HAVE_ENERGY && JOULE_EFFECT_AND_MAG_ENG)
    /* [Opt]: I could inglobate this call to ResEnergyIncrease in the previous one by using dt_res_reduced instead of 0.5*dt_res_reduced
      (but in this way it is more readable)*/
    // [Err] Decomment next line       
      ResEnergyIncrease(dUres_b1, H2p_B, H2m_B, Bra2, grid, &lines[DIR2], 0.5*dt_res_reduced, DIR2);
    #endif

    /**********************************
     (b.2) Implicit update sweeping DIR1
    **********************************/
    BoundaryRes_ADI(lines, d, grid, t_start_sub_res+dt_res_reduced); // Get bcs at t+dt   
    ImplicitUpdate (Brb2, Brb1, NULL, H1p_B, H1m_B, C1_B, &lines[DIR1],
                      lines[DIR1].lbound[BDIFF], lines[DIR1].rbound[BDIFF], 0.5*dt_res_reduced, DIR1);
    #if (HAVE_ENERGY && JOULE_EFFECT_AND_MAG_ENG)
    // [Err] Decomment next line
      ResEnergyIncrease(dUres_b2, H1p_B, H1m_B, Brb2, grid, &lines[DIR1], 0.5*dt_res_reduced, DIR1);
    #endif
    // [Err] Remove next #if lines
    // #if (HAVE_ENERGY && JOULE_EFFECT_AND_MAG_ENG)
      // ResEnergyIncrease(dUres_b1, H1p_B, H1m_B, Brb2, grid, &lines[DIR1], 0.5*dt_res_reduced, DIR1);
      // ResEnergyIncrease(dUres_b2, H2p_B, H2m_B, Brb2, grid, &lines[DIR2], 0.5*dt_res_reduced, DIR2);
    // #endif
    t_start_sub_res += dt_res_reduced;
  }
  #endif

  /* ------------------------------------------------------------
  ------------------------------------------------------------
    Update data
  ------------------------------------------------------------
  ------------------------------------------------------------ */
  KDOM_LOOP(k) {
    LINES_LOOP(lines[IDIR], l, j, i) {
      #if (RESISTIVITY == ALTERNATING_DIRECTION_IMPLICIT)
        Uc[k][j][i][BX3] = Brb2[j][i]*r_1[i];

        #if (HAVE_ENERGY && JOULE_EFFECT_AND_MAG_ENG)
        // [Err] Decomment next line
          print1("C'è un errore perchè non hai sommato bene tutti i dUresX");
          Uc[k][j][i][ENG] += dUres_a1[j][i]+dUres_a2[j][i]+dUres_b1[j][i]+dUres_b2[j][i];
        #endif
      #endif

      #if (THERMAL_CONDUCTION == ALTERNATING_DIRECTION_IMPLICIT)
        // I get the int. energy from the temperature
        #if EOS==IDEAL
          #error Not implemented for ideal eos (but it is easy to add it!)
        #elif EOS==PVTE_LAW
            #ifdef TEST_ADI
              for (nv=NVAR; nv--;) v[nv] = Vc[nv][k][j][i];

              /*I think in this way the update does not conserve the energy*/
              rhoe_old = 3/2*CONST_kB*v[RHO]*UNIT_DENSITY/CONST_mp*T[j][i]*KELVIN;
              rhoe_old /= (UNIT_DENSITY*UNIT_VELOCITY*UNIT_VELOCITY);
              rhoe_new = 3/2*CONST_kB*v[RHO]*UNIT_DENSITY/CONST_mp*Tb2[j][i]*KELVIN;
              rhoe_new /= (UNIT_DENSITY*UNIT_VELOCITY*UNIT_VELOCITY);
              Uc[k][j][i][ENG] += rhoe_new-rhoe_old;

              /*I think in this way the update should conserve the energy(23052018)*/
              // Uc[k][j][i][ENG] += dEdT[j][i]*(Tb2[j][i]-T[j][i]);
            #else
              /*I think in this way the update should conserve the energy*/
              Uc[k][j][i][ENG] += dEdT[j][i]*(Tb2[j][i]-T[j][i]);
            #endif
        #else
          print1("ADI:[Ema] Error computing internal energy, this EOS not implemented!");
        #endif
      #endif
    }
  }

  /* -------------------------------------------------------------------------
    Compute back the primitive vector from the updated conservative vector.
  ------------------------------------------------------------------------- */
  ConsToPrimLines (Uc, Vc, d->flag, lines);
}

/****************************************************************************
Performs an implicit update of a diffusive problem (either for B or for T)
*****************************************************************************/
void ImplicitUpdate (double **v, double **b, double **source,
                     double **Hp, double **Hm, double **C,
                     Lines *lines, Bcs *lbound, Bcs *rbound, double dt,
                     int dir) {
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
      // Now I solve the system
      tdm_solver( x+lidx, diagonal+lidx, upper+lidx, lower+lidx+1, rhs+lidx, ridx-lidx+1);
      /*[Opt] Is this for a waste of time? maybe I could engineer better the use of tdm_solver function(or the way it is written)*/
      for (i=lidx; i<=ridx; i++)
        v[j][i] = x[i];
    }
  } else if (dir == JDIR) {
    /********************
    * Case direction JDIR
    *********************/
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
      // Now I solve the system
      tdm_solver( x+lidx, diagonal+lidx, upper+lidx, lower+lidx+1, rhs+lidx, ridx-lidx+1);
      for (j=lidx; j<=ridx; j++)
        v[j][i] = x[j];
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
Performs an explicit update of a diffusive problem (either for B or for T)
*****************************************************************************/
void ExplicitUpdate (double **v, double **b, double **source,
                     double **Hp, double **Hm, double **C,
                     Lines *lines, Bcs *lbound, Bcs *rbound, double dt,
                     int dir) {
  int i,j,l;
  int ridx, lidx;
  int Nlines = lines->N;
  double b_ghost;
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

      // Actual update
      for (i = lidx+1; i < ridx; i++)
        v[j][i] = rhs[j][i] + dt/C[j][i] * (b[j][i+1]*Hp[j][i] - b[j][i]*(Hp[j][i]+Hm[j][i]) + b[j][i-1]*Hm[j][i]);
      // Cells near left boundary
      if (lbound[l].kind == DIRICHLET){
        b_ghost = 2*lbound[l].values[0] - b[j][lidx];
        v[j][lidx] = rhs[j][lidx] + dt/C[j][lidx] * (b[j][lidx+1]*Hp[j][lidx] - b[j][lidx]*(Hp[j][lidx]+Hm[j][lidx]) + b_ghost*Hm[j][lidx]);
      } else if (lbound[l].kind == NEUMANN_HOM) {
        v[j][lidx] = rhs[j][lidx] + dt/C[j][lidx] * (b[j][lidx+1]*Hp[j][lidx] - b[j][lidx]*Hp[j][lidx]);
      } else {
        print1("\n[ExplicitUpdate]Error setting left bc (in dir i), not known bc kind!");
        QUIT_PLUTO(1);
      }
      // Cells near right boundary
      if (rbound[l].kind == DIRICHLET){
        b_ghost = 2*rbound[l].values[0] - b[j][ridx];
        v[j][ridx] = rhs[j][ridx] + dt/C[j][ridx] * (b_ghost*Hp[j][ridx] - b[j][ridx]*(Hp[j][ridx]+Hm[j][ridx]) + b[j][ridx-1]*Hm[j][ridx]);
      } else if (rbound[l].kind == NEUMANN_HOM) {
        v[j][ridx] = rhs[j][ridx] + dt/C[j][ridx] * (-b[j][ridx]*Hm[j][ridx] + b[j][ridx-1]*Hm[j][ridx]);
      } else {
        print1("\n[ExplicitUpdate]Error setting right bc (in dir i), not known bc kind!");
        QUIT_PLUTO(1);
      }
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

      // Actual update
      for (j = lidx+1; j < ridx; j++)
        v[j][i] = rhs[j][i] + dt/C[j][i] * (b[j+1][i]*Hp[j][i] - b[j][i]*(Hp[j][i]+Hm[j][i]) + b[j-1][i]*Hm[j][i]);
      // Cells near left boundary
      if (lbound[l].kind == DIRICHLET){
        b_ghost = 2*lbound[l].values[0] - b[lidx][i];
        v[lidx][i] = rhs[lidx][i] + dt/C[lidx][i] * (b[lidx+1][i]*Hp[lidx][i] - b[lidx][i]*(Hp[lidx][i]+Hm[lidx][i]) + b_ghost*Hm[lidx][i]);
      } else if (lbound[l].kind == NEUMANN_HOM) {
        v[lidx][i] = rhs[lidx][i] + dt/C[lidx][i] * (b[lidx+1][i]*Hp[lidx][i] - b[lidx][i]*Hp[lidx][i]);
      } else {
        print1("\n[ExplicitUpdate]Error setting left bc (in dir j), not known bc kind!");
        QUIT_PLUTO(1);
      }
      // Cells near right boundary
      if (rbound[l].kind == DIRICHLET){
        b_ghost = 2*rbound[l].values[0] - b[ridx][i];
        v[ridx][i] = rhs[ridx][i] + dt/C[ridx][i] * (b_ghost*Hp[ridx][i] - b[ridx][i]*(Hp[ridx][i]+Hm[ridx][i]) + b[ridx-1][i]*Hm[ridx][i]);
      } else if (rbound[l].kind == NEUMANN_HOM) {
        v[ridx][i] = rhs[ridx][i] + dt/C[ridx][i] * (-b[ridx][i]*Hm[ridx][i] + b[ridx-1][i]*Hm[ridx][i]);
      } else {
        print1("\n[ExplicitUpdate]Error setting right bc (in dir j), not known bc kind!");
        QUIT_PLUTO(1);
      }
    }
  } else {
    print1("[ImplicitUpdate] Unimplemented choice for 'dir'!");
    QUIT_PLUTO(1);
  }
}

/****************************************************************************
Function to initialize lines, I hope this kind of initialization is ok and
there are no problems with data continuity and similar things
*****************************************************************************/
void InitializeLines(Lines *lines, int N){
  int i;

  lines->dom_line_idx = ARRAY_1D(N, int);
  lines->lidx = ARRAY_1D(N, int);
  lines->ridx = ARRAY_1D(N, int);
  lines->N = N;
  for (i=0; i<NADI; i++) {
    lines->lbound[i] = ARRAY_1D(N, Bcs);
    lines->rbound[i] = ARRAY_1D(N, Bcs);
  }
}

/****************************************************************************
Function to build geometrical parameters belonging to lines
[Rob] Maybe it's better to move this to another file? the one containing init()??
*****************************************************************************/
void GeometryADI(Lines *lines, Grid *grid){
  int i,j;

  // Use the number of rows and cols internal to the domain to define adi lines
  InitializeLines(&lines[IDIR], NX2); // I have NX2 lines sweeping the i-direction
  InitializeLines(&lines[JDIR], NX1); // I have NX1 lines sweeping the j-direction

  // A couple of cross-checks:
  if (2*grid[JDIR].nghost+NX2 != NX2_TOT) {
    print1("Something wrong with the line definition/not understood how the grid[DIR] is made!");
    QUIT_PLUTO(1);
  }
  if (2*grid[IDIR].nghost+NX1 != NX1_TOT) {
    print1("Something wrong with the line definition/not understood how the grid[JDIR] is made!");
    QUIT_PLUTO(1);
  }

  // Use remarkable capillary indexes to define the lines
  for (j=0;j<lines[IDIR].N;j++)
    lines[IDIR].dom_line_idx[j] = j + grid[JDIR].nghost;
  for (j=0;j<lines[IDIR].N;j++){
    lines[IDIR].lidx[j] = grid[IDIR].nghost;
    if (lines[IDIR].dom_line_idx[j] <= j_cap_inter_end){
      lines[IDIR].ridx[j] = i_cap_inter_end;
    } else {
      lines[IDIR].ridx[j] = NX1_TOT - 1 - grid[IDIR].nghost;
    }
  }
  for (i=0;i<lines[JDIR].N;i++)
    lines[JDIR].dom_line_idx[i] = i + grid[IDIR].nghost;
  for (i=0;i<lines[JDIR].N;i++){
    lines[JDIR].ridx[i] = NX2_TOT - 1 - grid[JDIR].nghost;
    if (lines[JDIR].dom_line_idx[i] <= i_cap_inter_end){
      lines[JDIR].lidx[i] = grid[JDIR].nghost;
    } else {
      lines[JDIR].lidx[i] = j_cap_inter_end+1;
    }
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
 * Peachman-Rachford ADI method.
 * 
 * input: diff = BDIFF or TDIFF
 *        rz = 1 or 0: tells whether the order of the directions
 *                     must be r, z (1) or z, r (0).
 * ***********************************************************/
void PeacemanRachford(double **v_new, double **v_old, double **dUres,
                      const Data *d, Grid *grid,
                      double **Ip, double **Im, double **CI,
                      double **Jp, double **Jm, double **CJ,
                      Lines *lines, int diff, int rz,
                      double dt, double t0) {
    
    static double **v_aux; // auxiliary solution vector
    static int first_call = 1;
    double **H1p, **H1m, **H2p, **H2m, **C1, **C2;
    double **dUres_aux; // auxiliary vector containing a contribution to ohmic heating
    int dir1, dir2;
    // void (*BoundaryADI) (Lines, const Data, Grid, double);
    BoundaryADI *ApplyBCs;
    int k,j,i,l;

    if (first_call) {
      v_aux = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      #if (HAVE_ENERGY && JOULE_EFFECT_AND_MAG_ENG)
        dUres_aux = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      #endif
      first_call = 0;
    }

    /* Set the direction order*/
    if (rz == 1) {
      H1p = Ip;     H1m = Im;
      H2p = Jp;     H2m = Jm;
      C1 = CI;      C2 = CJ;
      dir1 = IDIR;  dir2 = JDIR;
    } else if (rz == 0) {
      H1p = Jp;     H1m = Jm;
      H2p = Ip;     H2m = Im;
      C1 = CJ;      C2 = CI;
      dir1 = JDIR;  dir2 = IDIR;
    }

    if (diff == BDIFF)
      ApplyBCs = BoundaryRes_ADI;
    else if (diff == TDIFF)
      ApplyBCs = BoundaryTC_ADI;
    else {
      print1("\n[PeachmanRachford]Wrong setting for diffusion (diff) problem");
      QUIT_PLUTO(1);
    }

    /**********************************
     (a.1) Explicit update sweeping DIR1
    **********************************/
    ApplyBCs(lines, d, grid, t0);
    // [Err] Remove next #if lines
    // #if (HAVE_ENERGY && JOULE_EFFECT_AND_MAG_ENG)
      // ResEnergyIncrease(dUres_a1, H1p, H1m, Br, grid, &lines[DIR1], 0.5*dt_res_reduced, DIR1);
      // ResEnergyIncrease(dUres_a2, H2p, H2m, Br, grid, &lines[DIR2], 0.5*dt_res_reduced, DIR2);
    // #endif
    ExplicitUpdate (v_aux, v_old, NULL, H1p, H1m, C1, &lines[dir1],
                    lines[dir1].lbound[diff], lines[dir1].rbound[diff], 0.5*dt, dir1);
    #if (HAVE_ENERGY && JOULE_EFFECT_AND_MAG_ENG)
      if (diff == BDIFF) {
        // [Err] Decomment next line
        // [Opt] You could modify and make that the ResEnergyEncrease automatically updates a Ures variable,
        //       instead of doing it a line later 
        ResEnergyIncrease(dUres_aux, H1p, H1m, v_aux, grid, &lines[dir1], 0.5*dt, dir1);
        KDOM_LOOP(k)
          LINES_LOOP(lines[IDIR], l, j, i)
            dUres[j][i] += dUres_aux[j][i];
      }
    #endif

    /**********************************
     (a.2) Implicit update sweeping DIR2
    **********************************/
    ApplyBCs(lines, d, grid, t0 + dt*0.5);
    ImplicitUpdate (v_new, v_aux, NULL, H2p, H2m, C2, &lines[dir2],
                      lines[dir2].lbound[diff], lines[dir2].rbound[diff], 0.5*dt, dir2);
    #if (HAVE_ENERGY && JOULE_EFFECT_AND_MAG_ENG)
      if (diff == BDIFF) {
        // [Err] Decomment next line
        ResEnergyIncrease(dUres_aux, H2p, H2m, v_new, grid, &lines[dir2], 0.5*dt, dir2);
        KDOM_LOOP(k)
          LINES_LOOP(lines[IDIR], l, j, i)
            dUres[j][i] += dUres_aux[j][i];
      }
    #endif
    // [Err] Remove next four lines
    // ResEnergyIncrease(dUres_a1, Ip, Im, Bra2, grid, &lines[IDIR], 0.5*dt_res_reduced, IDIR);
    // ResEnergyIncrease(dUres_a2, Jp, Jm, Bra2, grid, &lines[DIR2], 0.5*dt_res_reduced, DIR2);
    // ResEnergyIncrease(dUres_b1, Ip, Im, Bra2, grid, &lines[IDIR], 0.5*dt_res_reduced, IDIR);
    // ResEnergyIncrease(dUres_b2, Jp, Jm, Bra2, grid, &lines[DIR2], 0.5*dt_res_reduced, DIR2);

    /**********************************
     (b.1) Explicit update sweeping DIR2
    **********************************/
    ExplicitUpdate (v_aux, v_new, NULL, H2p, H2m, C2, &lines[dir2],
                    lines[dir2].lbound[diff], lines[dir2].rbound[diff], 0.5*dt, dir2);
    #if (HAVE_ENERGY && JOULE_EFFECT_AND_MAG_ENG)
      if (diff == BDIFF) {
        /* [Opt]: I could inglobate this call to ResEnergyIncrease in the previous one by using dt_res_reduced instead of 0.5*dt_res_reduced
           (but in this way it is more readable)*/
        // [Err] Decomment next line       
        ResEnergyIncrease(dUres_aux, H2p, H2m, v_aux, grid, &lines[DIR2], 0.5*dt, DIR2);
        KDOM_LOOP(k)
          LINES_LOOP(lines[IDIR], l, j, i)
            dUres[j][i] += dUres_aux[j][i];
      }    
    #endif

    /**********************************
     (b.2) Implicit update sweeping DIR1
    **********************************/
    ApplyBCs(lines, d, grid, t0 + dt);
    ImplicitUpdate (v_new, v_aux, NULL, H1p, H1m, C1, &lines[dir1],
                      lines[dir1].lbound[diff], lines[dir1].rbound[diff], 0.5*dt, dir1);
    #if (HAVE_ENERGY && JOULE_EFFECT_AND_MAG_ENG)
      if (diff == BDIFF) {
        // [Err] Decomment next line
        ResEnergyIncrease(dUres_aux, H1p, H1m, v_new, grid, &lines[dir1], 0.5*dt, dir1);
        KDOM_LOOP(k)
          LINES_LOOP(lines[IDIR], l, j, i)
            dUres[j][i] += dUres_aux[j][i];
      }
    #endif
    // [Err] Remove next #if lines
    // #if (HAVE_ENERGY && JOULE_EFFECT_AND_MAG_ENG)
      // ResEnergyIncrease(dUres_b1, H1p, H1m, Brb2, grid, &lines[DIR1], 0.5*dt_res_reduced, DIR1);
      // ResEnergyIncrease(dUres_b2, H2p, H2m, Brb2, grid, &lines[DIR2], 0.5*dt_res_reduced, DIR2);
    // #endif
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
