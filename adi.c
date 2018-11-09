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

#if THERMAL_CONDUCTION == ALTERNATING_DIRECTION_IMPLICIT
  // Temperature, to make it available outside (by means of a function)
  static double **T_old;
#endif

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
    static double **T_new;
    static double **dEdT;
    double v[NVAR]; /*[Ema] I hope that NVAR as dimension is fine!*/
    // double rhoe_old, rhoe_new;
    int nv;
  #endif

  int s;
  double t_start_sub;
  double dt_reduced;  
  int const adi_steps = NSUBS_ADI_TOT;
  const double dt = g_dt;
  double ****Uc, ****Vc;
  double *r, *r_1;

  #if FIRST_JDIR_THEN_IDIR == RANDOM
    if (first_call)
      srand(time(NULL));   // should only be called once
  #elif FIRST_JDIR_THEN_IDIR == AVERAGE
    #if METHOD_RES==FRACTIONAL_THETA||METHOD_RES==SPLIT_IMPLICIT||METHOD_TC==FRACTIONAL_THETA||METHOD_TC==SPLIT_IMPLICIT 
      #error Average order is not yet implemented for FractionalTheta or SplitImplicit
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
      T_old = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      dEdT = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      TOT_LOOP (k,j,i) {
        T_new[j][i] = 0.0;
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
  PrimToConsLines (Vc, Uc, lines);

  #if RESISTIVITY == ALTERNATING_DIRECTION_IMPLICIT
    DOM_LOOP(k,j,i) {
      Br_old[j][i] = r[i]*Uc[k][j][i][BX3];
    }
  #endif

  t_start_sub = g_time; /*g_time è: "The current integration time."(dalla docuementazione in Doxigen)*/
  dt_reduced = dt/adi_steps;

  for (s=0; s<adi_steps; s++) {
    #ifdef DEBUG_EMA
      printf("\nNstep:%ld",g_stepNumber);
      printf("\ns:%d\n", s);
    #endif
    Boundary(d, ALL_DIR, grid);

    /* ---- Build temperature vector ---- */
    // I must re-buil the temperature at every step as it depends on U[][][][ENG]
    #if THERMAL_CONDUCTION == ALTERNATING_DIRECTION_IMPLICIT
      #if EOS==PVTE_LAW
        DOM_LOOP(k,j,i) {
          for (nv=NVAR; nv--;) v[nv] = Vc[nv][k][j][i];
          if (GetPV_Temperature(v, &(T_old[j][i]) )!=0) {
            #if WARN_ERR_COMP_TEMP
              print1("ADI:[Ema]Err.comp.temp\n");
            #endif
          }
          T_old[j][i] = T_old[j][i] / KELVIN;
        }
      #else
        print1("ADI:[Ema]Err.comp.temp, this EOS not implemented!")
      #endif

      /* ---- Avdance T with ADI ---- */
      #if METHOD_TC==SPLIT_IMPLICIT
        #error SPLIT_IMPLICIT has not yet been tested with thermal conduction
        if (NSUBS_TC!=1) {
          print1("\n[ADI] In SPLIT_IMPLICIT method only NSUBS_TC=1 is implemented");
          QUIT_PLUTO(1);
        }
      #elif METHOD_TC==FRACTIONAL_THETA
        if (NSUBS_TC!=1) {
          print1("\n[ADI] In FRACTIONAL_THETA method only NSUBS_TC=1 is implemented");
          QUIT_PLUTO(1);
        }
        FractionalTheta(T_new, T_old, NULL, dEdT, d, grid, lines, TDIFF, ORDER, dt_reduced, t_start_sub, FRACTIONAL_THETA_THETA_TC);
      #elif METHOD_TC==DOUGLAS_RACHFORD
        DouglasRachford(T_new, T_old, NULL, dEdT, d, grid, lines, TDIFF, ORDER, dt_reduced, t_start_sub, NSUBS_TC);
      #elif METHOD_TC==PEACEMAN_RACHFORD_MOD
        PeacemanRachfordMod(T_new, T_old, NULL, dEdT, d, grid, lines, TDIFF, ORDER, dt_reduced, t_start_sub, FRACT_TC, NSUBS_TC);
      #elif METHOD_TC==STRANG_LIE
        #error STRANG_LIE has not yet been tested with thermal conduction
      #elif METHOD_TC==STRANG
        #error STRANG has not yet been tested with thermal conduction
      #else
        print1("[ADI]No suitable scheme for thermal conduction has been selected!");
        QUIT_PLUTO(1);
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
                rhoe_old = 3/2*CONST_kB*v[RHO]*UNIT_DENSITY/CONST_mp*T_old[j][i]*KELVIN;
                rhoe_old /= (UNIT_DENSITY*UNIT_VELOCITY*UNIT_VELOCITY);
                rhoe_new = 3/2*CONST_kB*v[RHO]*UNIT_DENSITY/CONST_mp*T_new[j][i]*KELVIN;
                rhoe_new /= (UNIT_DENSITY*UNIT_VELOCITY*UNIT_VELOCITY);
                Uc[k][j][i][ENG] += rhoe_new-rhoe_old;

                /*I think in this way the update should conserve the energy(23052018)*/
                // Uc[k][j][i][ENG] += dEdT[j][i]*(T_new[j][i]-T_old[j][i]);
              #else
                #ifdef DEBUG_TNEGATIVE
                  if (T_new[j][i]<=0.0) {
                    print1("! ADI: T_new[%d][%d] = %g < 0", j,i,T_new[j][i] );
                  }
                #endif
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
      #if METHOD_RES==SPLIT_IMPLICIT
        SplitImplicit(Br_new, Br_old, dUres, NULL, d, grid, lines, BDIFF, ORDER, dt_reduced, t_start_sub, NSUBS_RES);
      #elif METHOD_RES==FRACTIONAL_THETA
        if (NSUBS_RES!=1) {
          print1("\n[ADI] In FRACTIONAL_THETA method only NSUBS_RES=1 is implemented");
          QUIT_PLUTO(1);
        }
        FractionalTheta(Br_new, Br_old, dUres, NULL, d, grid, lines, BDIFF, ORDER, dt_reduced, t_start_sub, FRACTIONAL_THETA_THETA_RES);
      #elif METHOD_RES==DOUGLAS_RACHFORD
        DouglasRachford(Br_new, Br_old, dUres, NULL, d, grid, lines, BDIFF, ORDER, dt_reduced, t_start_sub, NSUBS_RES);
      #elif METHOD_RES==PEACEMAN_RACHFORD_MOD
        PeacemanRachfordMod(Br_new, Br_old, dUres, NULL, d, grid, lines, BDIFF, ORDER, dt_reduced, t_start_sub, FRACT_RES, NSUBS_RES);
      #elif METHOD_RES==STRANG_LIE
        Strang_Lie (Br_new, Br_old, dUres, NULL, d, grid, lines, BDIFF, ORDER, dt_reduced, t_start_sub, NSUBS_RES);
      #elif METHOD_RES==STRANG
        Strang (Br_new, Br_old, dUres, NULL, d, grid, lines, BDIFF, ORDER, dt_reduced, t_start_sub, NSUBS_RES);      
      #else
        print1("[ADI]No suitable scheme for resistivity has been selected!");
        QUIT_PLUTO(1);
      #endif
      

      /* ---- Update cons variables ---- */
      KDOM_LOOP(k)
        LINES_LOOP(lines[IDIR], l, j, i) {
          Uc[k][j][i][BX3] = Br_new[j][i]*r_1[i];

          #if (JOULE_EFFECT_AND_MAG_ENG && (!MAG_PS_OUTSIDE_SSTEP))
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
  }

  // Update the time where the diffusion process has arrived
  t_diff = t_start_sub;
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

#if THERMAL_CONDUCTION == ALTERNATING_DIRECTION_IMPLICIT
  /* ***********************************************************
  * Function to get T_old outside this file
  * ***********************************************************/
  double GetT_old (int j, int i) {
    return T_old[j][i];
  }
#endif

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
