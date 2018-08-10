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

// I initialize the diffusion time, since it is nedded before the diffusion starts;
double t_diff = 0;

void ADI(const Data *d, Time_Step *Dts, Grid *grid) {
  static int first_call=1;
  int i,j,k, l;
  
  static Lines lines[2]; /*I define two of them as they are 1 per direction (r and z)*/
  #if RESISTIVITY == ALTERNATING_DIRECTION_IMPLICIT
    static double **Br_new, **Br_old;
    // static double **Br_avg;
    // Energy increse(due to electro-magnetics) terms
    static double **dUres;
  #endif
  #if THERMAL_CONDUCTION == ALTERNATING_DIRECTION_IMPLICIT
    static double **T_new, **T_old, **T;
    static double **dEdT;
    double v[NVAR]; /*[Ema] I hope that NVAR as dimension is fine!*/
    // double rhoe_old, rhoe_new;
    int nv;
  #endif

  int s;
  double t_start_sub;
  double dt_reduced;  
  int const adi_steps = NSUBS_ADI;
  const double dt = g_dt;
  double ****Uc, ****Vc;
  double *r, *r_1;

  #ifdef PRESUBS_RES
    #if  THERMAL_CONDUCTION == ALTERNATING_DIRECTION_IMPLICIT
      #error Thermal conduction is not compatible with PRESUBS_RES
    #endif
    double dt_test = 0.1*2.7964e-8;
    int Npresubs = PRESUBS_RES;
    double t_start_presub;
  #else
    int Npresubs = 0;
    double dt_test = 0.0;
  #endif

  #if FIRST_JDIR_THEN_IDIR == RANDOM
    if (first_call)
      srand(time(NULL));   // should only be called once
  #elif FIRST_JDIR_THEN_IDIR == AVERAGE
    #if (defined(FRACTIONAL_THETA) || defined(SPLIT_IMPLICIT))
      #error Average order is not yet implemented for FractionalTheta and SplitImplicit
    #endif
    static double **Br_new_other_order; //Additional Br result when the other order of direction is used
    static double **dUres_other_order; //Additional dUres result when the other order of direction is used
  #endif

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
      Br_new = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      Br_old = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      // Br_avg = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      dUres = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      
      // This is just for debug purposes
      TOT_LOOP (k,j,i) {
        Br_new[j][i] = 0.0;
        dUres[j][i] = 0.0;
      }
      #if FIRST_JDIR_THEN_IDIR == AVERAGE
        Br_new_other_order = ARRAY_2D(NX2_TOT, NX1_TOT, double);
        dUres_other_order = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      #endif
    #endif

    #if THERMAL_CONDUCTION == ALTERNATING_DIRECTION_IMPLICIT
      T_new = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      T = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      T_old = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      dEdT = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      TOT_LOOP (k,j,i) {
        T_new[j][i] = 0.0;
        T[j][i] = 0.0;
      }
    #endif

    first_call=0;
  }

  /* -------------------------------------------------------------------------
      Compute the conservative vector in order to start the cycle.
      This step will be useless if the data structure
      contains Vc as well as Uc (for future improvements).
      [Ema] (Comment copied from sts.c)
    --------------------------------------------------------------------------- */
  //[Err] Remove next line
  Boundary(d, ALL_DIR, grid);
  PrimToConsLines (Vc, Uc, lines);

  #if RESISTIVITY == ALTERNATING_DIRECTION_IMPLICIT
    DOM_LOOP(k,j,i) {
      Br_old[j][i] = r[i]*Uc[k][j][i][BX3];
    }
  #endif

  #ifdef PRESUBS_RES
    t_start_presub = g_time;
    for (s=0; s<Npresubs; s++) {
      PeacemanRachfordMod(Br_new, Br_old, dUres, NULL, d, grid, lines, BDIFF, ORDER, dt_test, t_start_presub, 0.5);
      /* ---- Update cons variables ---- */
      KDOM_LOOP(k)
        LINES_LOOP(lines[IDIR], l, j, i) {
          Uc[k][j][i][BX3] = Br_new[j][i]*r_1[i];

          #if (JOULE_EFFECT_AND_MAG_ENG)
            // [Err] Decomment next line
            Uc[k][j][i][ENG] += dUres[j][i];
          #endif
        }
      /* ---- Swap pointers to be ready for next cycle ---*/
      SwapDoublePointers (&Br_new, &Br_old);
      /* -------------------------------------------------------------------------
        Compute back the primitive vector from the updated conservative vector.
        ------------------------------------------------------------------------- */
      ConsToPrimLines (Uc, Vc, d->flag, lines);
      Boundary(d, ALL_DIR, grid);
      t_start_presub += dt_test;
    }
  #endif

  t_start_sub = g_time+dt_test*Npresubs; /*g_time è: "The current integration time."(dalla docuementazione in Doxigen)*/
  #ifndef COMMON_RATIO_NSUBS_ADI
    dt_reduced = (dt-dt_test*Npresubs)/adi_steps;
  #else
    dt_reduced = (dt-dt_test*Npresubs) / ((1-pow(COMMON_RATIO_NSUBS_ADI,adi_steps))/(1-COMMON_RATIO_NSUBS_ADI));
  #endif

  for (s=0; s<adi_steps; s++) {
    #ifdef DEBUG_EMA
      printf("\nNstep:%ld",g_stepNumber);
      printf("\ns:%d\n", s);
    #endif
    // [Err] Test, decomment later
    // Boundary(d, ALL_DIR, grid);

    /* ---- Build temperature vector ---- */
    // I must re-buil the temperature at every step as it depends on U[][][][ENG]
    #if THERMAL_CONDUCTION == ALTERNATING_DIRECTION_IMPLICIT
      #if EOS==PVTE_LAW
        DOM_LOOP(k,j,i) {
          for (nv=NVAR; nv--;) v[nv] = Vc[nv][k][j][i];
          if (GetPV_Temperature(v, &(T[j][i]) )!=0) {
            print1("ADI:[Ema] Error computing temperature!\n");
          }
          T_old[j][i] = T[j][i] = T[j][i] / KELVIN;
        }
      #else
        print1("ADI:[Ema] Error computing temperature, this EOS not implemented!")
      #endif

      /* ---- Avdance T with ADI ---- */
      #ifdef FRACTIONAL_THETA
        FractionalTheta(T_new, T_old, NULL, dEdT, d, grid, lines, TDIFF, ORDER, dt_reduced, t_start_sub, FRACTIONAL_THETA);
      #else
        PeacemanRachford(T_new, T_old, NULL, dEdT, d, grid, lines, TDIFF, ORDER, dt_reduced, t_start_sub);
      #endif

      /* ---- Update cons variables ---- */
      KDOM_LOOP(k)
        LINES_LOOP(lines[IDIR], l, j, i) {
          // I get the int. energy from the temperature
          #if EOS==IDEAL
            #error Not implemented for ideal eos (but it is easy to add it!)
          #elif EOS==PVTE_LAW
              #ifdef TEST_ADI
                for (nv=NVAR; nv--;) v[nv] = Vc[nv][k][j][i];

                /*I think in this way the update does not conserve the energy*/
                rhoe_old = 3/2*CONST_kB*v[RHO]*UNIT_DENSITY/CONST_mp*T[j][i]*KELVIN;
                rhoe_old /= (UNIT_DENSITY*UNIT_VELOCITY*UNIT_VELOCITY);
                rhoe_new = 3/2*CONST_kB*v[RHO]*UNIT_DENSITY/CONST_mp*T_new[j][i]*KELVIN;
                rhoe_new /= (UNIT_DENSITY*UNIT_VELOCITY*UNIT_VELOCITY);
                Uc[k][j][i][ENG] += rhoe_new-rhoe_old;

                /*I think in this way the update should conserve the energy(23052018)*/
                // Uc[k][j][i][ENG] += dEdT[j][i]*(T_new[j][i]-T[j][i]);
              #else
                /*I think in this way the update should conserve the energy*/
                Uc[k][j][i][ENG] += dEdT[j][i]*(T_new[j][i]-T_old[j][i]);
              #endif
          #else
            print1("ADI:[Ema] Error computing internal energy, this EOS not implemented!");
          #endif
        }

      /* ---- Swap pointers to be ready for next cycle ---*/
      SwapDoublePointers (&T_new, &T_old);
    #endif

    #if RESISTIVITY == ALTERNATING_DIRECTION_IMPLICIT

      // No need to re-build the magnetic field, as it does not depend on U[][][][whatever]

      /* ---- Avdance B*r with ADI ---- */
      #ifdef SPLIT_IMPLICIT
        SplitImplicit(Br_new, Br_old, dUres, NULL, d, grid, lines, BDIFF, ORDER, dt_reduced, t_start_sub);
      #elif defined(FRACTIONAL_THETA)
        FractionalTheta(Br_new, Br_old, dUres, NULL, d, grid, lines, BDIFF, ORDER, dt_reduced, t_start_sub, FRACTIONAL_THETA);
      #else
        PeacemanRachfordMod(Br_new, Br_old, dUres, NULL, d, grid, lines, BDIFF, ORDER, dt_reduced, t_start_sub, FRACT);
        // DouglasRachford(Br_new, Br_old, dUres, NULL, d, grid, lines, BDIFF, ORDER, dt_reduced, t_start_sub);
      #endif
      

      /* ---- Update cons variables ---- */
      KDOM_LOOP(k)
        LINES_LOOP(lines[IDIR], l, j, i) {
          Uc[k][j][i][BX3] = Br_new[j][i]*r_1[i];

          #if (JOULE_EFFECT_AND_MAG_ENG)
            // [Err] Decomment next line
            Uc[k][j][i][ENG] += dUres[j][i];
          #endif
        }

      /* ---- Swap pointers to be ready for next cycle ---*/
      SwapDoublePointers (&Br_new, &Br_old);
    #endif

    /* -------------------------------------------------------------------------
        Compute back the primitive vector from the updated conservative vector.
        ------------------------------------------------------------------------- */
    ConsToPrimLines (Uc, Vc, d->flag, lines);

    t_start_sub += dt_reduced;
    //[Err] Remove next line
    Boundary(d, ALL_DIR, grid);

    #ifdef COMMON_RATIO_NSUBS_ADI
    // I update the dt_reduced
      dt_reduced *= COMMON_RATIO_NSUBS_ADI;
    #endif
  }
  // Update the time where the diffusion process has arrived
  t_diff = t_start_sub;

  //[Err] Test, delete later
  if (g_stepNumber == 3) {
    DumpQuit (d, RuntimeGet(), grid);
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

/* ***********************************************************
 * Function to swap double pointers to double
 * ***********************************************************/
void SwapDoublePointers (double ***a, double ***b) {
  double **temp;
  temp = *a;
  *a = *b;
  *b = temp;
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
