/*Functions and utilities for integration of parabolic (currently resistive
and thermal conduction) terms with the Alternating Directino Implicit algorithm*/
#include "pluto.h"
#include "adi.h"
#include "capillary_wall.h"

void ADI(const Data *d, Time_Step *Dts, Grid *grid) {
  static int first_call=1;
  int i,j,k, l, nv;
  static Lines lines[2]; /*I define two of them as they are 1 per direction (r and z)*/
  #if RESISTIVITY == ALTERNATING_DIRECTION_IMPLICIT
    static double **Ip_B, **Im_B, **Jp_B, **Jm_B, **CB;
    static double **Ba1, **Ba2, **Bb1, **Bb2;
    static double **B;
  #endif
  #if THERMAL_CONDUCTION == ALTERNATING_DIRECTION_IMPLICIT
    static double **Ip_T, **Im_T, **Jp_T, **Jm_T, **CT;
    static double **Ta1, **Ta2, **Tb1, **Tb2;
    static double **T;
  #endif
  static double **sourcea1, **sourcea2, **sourceb1, **sourceb2;
  // static double **dEdT;
  double dt;
  double v[NVAR]; /*[Ema] I hope that NVAR as dimension is fine!*/
  double ****V;

  // Find the remarkable indexes (if they had not been found before)
  if (capillary_not_set) {
    if (SetRemarkableIdxs(grid)){
      print1("\nError while setting remarkable points!");
      QUIT_PLUTO(1);
    }
  }
  // Build geometry and allocate some stuff
  if (first_call) {
    GeometryADI(lines, grid);
    #if RESISTIVITY == ALTERNATING_DIRECTION_IMPLICIT
      Ip_B = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      Im_B = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      Jp_B = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      Jm_B = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      CB = ARRAY_2D(NX2_TOT, NX1_TOT, double);

      Ba1 = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      Ba2 = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      Bb1 = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      Bb2 = ARRAY_2D(NX2_TOT, NX1_TOT, double);

      B = ARRAY_2D(NX2_TOT, NX1_TOT, double);
    #endif
    #if THERMAL_CONDUCTION == ALTERNATING_DIRECTION_IMPLICIT
      Ip_T = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      Im_T = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      Jp_T = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      Jm_T = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      CT = ARRAY_2D(NX2_TOT, NX1_TOT, double);

      Ta1 = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      Ta2 = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      Tb1 = ARRAY_2D(NX2_TOT, NX1_TOT, double);
      Tb2 = ARRAY_2D(NX2_TOT, NX1_TOT, double);

      T = ARRAY_2D(NX2_TOT, NX1_TOT, double);
    #endif

    sourcea1 = ARRAY_2D(NX2_TOT, NX1_TOT, double);
    sourcea2 = ARRAY_2D(NX2_TOT, NX1_TOT, double);
    sourceb1 = ARRAY_2D(NX2_TOT, NX1_TOT, double);
    sourceb2 = ARRAY_2D(NX2_TOT, NX1_TOT, double);
    
    // dEdT = ARRAY_2D(NX2_TOT, NX1_TOT, double);
    first_call=0;
  }

  // Build a handy magnetic field matrix
  DOM_LOOP(k,j,i)
    B[j][i] = d->Uc[BX3][k][j][i];

  // Build the temperature matrix
  #if EOS==IDEAL
      DOM_LOOP(k,j,i) T[j][i] = d->Vc[PRS][k][j][i]/V[RHO][k][j][i];
    #elif EOS==PVTE_LAW
      DOM_LOOP(k,j,i) {
        for (nv=NVAR; nv--;) v[nv] = d->Vc[nv][k][j][i];
        if (GetPV_Temperature(v, &(T[j][i]) )!=0) {
          print1("ADI:[Ema] Error computing temperature!\n");
        }
        T[j][i] = T[j][i] / KELVIN;
      }
    #else
      print1("ADI:[Ema] Error computing temperature, this EOS not implemented!")
    #endif

  BuildIJ(d, grid, Ip_T, Im_T, Jp_T, Jm_T, CT, BDIFF);
  BuildIJ(d, grid, Ip_B, Im_B, Jp_B, Jm_B, CB, TDIFF);
  BoundaryADI(); // Get bcs at t

  /**********************************
   (a.1) Explicit update sweeping IDIR
  **********************************/
  #if RESISTIVITY == ALTERNATING_DIRECTION_IMPLICIT
    ExplicitUpdate( Ba1, B, Ip_B, Im_B, CB, &lines[IDIR], NULL, 0.5*dt);
  #endif
  #if RESISTIVITY == ALTERNATING_DIRECTION_IMPLICIT  // Sure only if it is adi?? maybe it's ok even if it is sts or expl
    // Include eta*J^2 source term using B
    // ...
  #else
    LINES_LOOP(lines[0], l, j, i)
      sourcea1[j][i] = 0.0;
  #endif
  #if THERMAL_CONDUCTION == ALTERNATING_DIRECTION_IMPLICIT
    ExplicitUpdate( Ta1, T, Ip_T, Im_T, CT, &lines[IDIR], sourcea1, 0.5*dt);
  #endif

  /**********************************
   (a.2) Implicit update sweeping JDIR
  **********************************/
  BoundaryADI(); // Get bcs at half step (not exaclty at t+0.5*dt)
  #if RESISTIVITY == ALTERNATING_DIRECTION_IMPLICIT
    ImplicitUpdate( Ba2, Ba1, Jp_B, Jm_B, CB, &lines[JDIR], NULL, 0.5*dt);
  #endif
  #if RESISTIVITY == ALTERNATING_DIRECTION_IMPLICIT
    // Build eta*J^2 source term using Ba2
    // ...
  #else
    LINES_LOOP(lines[0], l, j, i)
      sourcea2[j][i] = 0.0;
  #endif
  #if THERMAL_CONDUCTION == ALTERNATING_DIRECTION_IMPLICIT
    ImplicitUpdate( Ta2, Ta1, Jp_T, Jm_T, CT, &lines[JDIR], sourcea2, 0.5*dt);
  #endif

  /**********************************
   (b.1) Explicit update sweeping JDIR
  **********************************/
  #if RESISTIVITY == ALTERNATING_DIRECTION_IMPLICIT
    ExplicitUpdate(Bb1, Ba2, Jp_B, Jm_B, CB, lines, NULL, 0.5*dt);
  #endif
  // /* -- This is USELESS as I already computed the source before, delete in future*/
  // #if RESISTIVITY == ALTERNATING_DIRECTION_IMPLICIT
  //   // Build eta*J^2 source term using Ba2
  //   // ...
  // #else
  //   LINES_LOOP(lines[0], l, j, i)
  //     sourceb1[j][i] = 0.0;
  // #endif
  #if THERMAL_CONDUCTION == ALTERNATING_DIRECTION_IMPLICIT
    ExplicitUpdate(Tb1, Ta2, Jp_T, Jm_T, CT, lines, sourcea1, 0.5*dt);
  #endif

  /**********************************
   (b.2) Implicit update sweeping IDIR
  **********************************/
  BoundaryADI(); // Get bcs at t+dt
  #if RESISTIVITY == ALTERNATING_DIRECTION_IMPLICIT
    ImplicitUpdate(Bb2, Bb1, Ip_B, Im_B, CB, lines, NULL, 0.5*dt);
  #endif
  #if RESISTIVITY == ALTERNATING_DIRECTION_IMPLICIT // Sure only if it is adi?? maybe it's ok even if it is sts or expl
    // Build eta*J^2 source term using Bb2
    // ...
  #else
    LINES_LOOP(lines[0], l, j, i)
      sourceb2[j][i] = 0.0;
  #endif
  #if THERMAL_CONDUCTION == ALTERNATING_DIRECTION_IMPLICIT
    ImplicitUpdate(Tb2, Tb1, Ip_T, Im_T, CT, lines, sourcea1, 0.5*dt);
  #endif

  /***********************************
   * Update data
   * *********************************/
  KDOM_LOOP(k) {
    for (j=0;j<lines[IDIR].N;j++){
      #if (RESISTIVITY == ALTERNATING_DIRECTION_IMPLICIT)
        d->Uc[k][j][i][BX3] = Bb2[j][i];
      #endif

      #if HAVE_ENERGY

        #if (THERMAL_CONDUCTION == ALTERNATING_DIRECTION_IMPLICIT)
          // I get the int. energy from the temperature
          #if EOS==IDEAL
            #error Not implemented for ideal eos (but it is easy to add it!)
          #elif EOS==PVTE_LAW
            DOM_LOOP(k,j,i) {
              for (nv=NVAR; nv--;) v[nv] = d->Vc[nv][k][j][i];
              if (d->Uc[k][j][i][ENG] = InternalEnergyFunc(v, Tb2[j][i]))
                print1("ADI:[Ema] Error computing internal energy!\n");
            }
          #else
            print1("ADI:[Ema] Error computing internal energy, this EOS not implemented!")
          #endif
        
        #elif (RESISTIVITY      == ALTERNATING_DIRECTION_IMPLICIT) && \
              (THERMAL_CONDUCTION != ALTERNATING_DIRECTION_IMPLICIT)
          // I provide the increase in int.energy from the integration of the pow.source
          /* I am just using trapezi integration of source term,
          but I could use something like cav-simps with very little effort
          (I have already sourcea, sourceb, sourcec(?)!*/
          d->Uc[k][j][i][ENG] += sourceb1[j][i]*dt
        #endif
      #endif
    }
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

/****************************************************************************
Function to build the bcs of lines
*****************************************************************************/
void BoundaryADI() {
  
}

/****************************************************************************
Function to build the Ip,Im,Jp,Jm
*****************************************************************************/
void BuildIJ(Data *d, Grid *grid, double **Ip, double **Im, double **Jp,
             double **Jm, double **C, int kind) {
  
}

/****************************************************************************
Performs an implicit update of a diffusive problem (either for B or for T)
*****************************************************************************/
void ImplicitUpdate( double **v, double **rhs, double **Hp, double **Hm, double **C,
                     Lines *lines, double **source, double dt) {

}

/****************************************************************************
Performs an explicit update of a diffusive problem (either for B or for T)
*****************************************************************************/
void ExplicitUpdate( double **v, double **rhs, double **Hp, double **Hm, double **C,
                     Lines *lines, double **source, double dt) {

}