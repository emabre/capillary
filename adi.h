#ifndef ADI_H
#define ADI_H
/* Contains some structures used for the adi method*/

#if KBEG != KEND
  #error grid in k direction should only be of 1 point
#endif

#if !HAVE_ENERGY && JOULE_EFFECT_AND_MAG_ENG
  #error You need an energy to accunt for Joule effect and mag. energy
#endif

#if (JOULE_EFFECT_AND_MAG_ENG && RESISTIVITY!=ALTERNATING_DIRECTION_IMPLICIT)
  #error Joule effect and mag. energy requires ADI for resistivity
#endif

#if (POW_INSIDE_ADI != YES)
  // POW_INSIDE_ADI != YES is not allowed now since I would involve a correct averaging of the bcs between timesteps
  // Note that averaging the bcs correctly is complicated.. how do you treat the internal corner cell(which represents 2 bcs!)?
  // and also, how do you compute the BCs of time steps where you did not use the bcs to advance the solution??
  #error Did you check carefully that you can do this? (e.g.: is the averaging of the bcs/ghost cells.. ok?)
#endif

#ifdef PRESUBS_RES
  #if THERMAL_CONDUCTION == ALTERNATING_DIRECTION_IMPLICIT
    #error PRESUBS_RES is not yet fully compatible with ADI thermal conduction
  #endif
#endif

/*********************
 * Some useful macros
 * *******************/
#define DIRICHLET    1
#define NEUMANN_HOM  2

// macro for calling RuntimeSet()
#define AFTER_SETOUTPUT 1

/* To loop on lines, the correct use is, usually:
   LINES_LOOP(lines[IDIR], l, j, i)
   or:
   LINES_LOOP(lines[JDIR], l, i, j)*/
#define LINES_LOOP(lines, l, line_idx , line_sweeper) \
  for ((line_idx)=lines.dom_line_idx[(l)=0]; (l)<lines.N; (line_idx)=lines.dom_line_idx[++(l)]) \
  for ((line_sweeper)=lines.lidx[(l)]; (line_sweeper)<=lines.ridx[(l)]; (line_sweeper)++)

// Macros which are valid choices for setting FIRST_JDIR_THEN_IDIR (in addition, YES and NO are valid too)
#define RANDOM 2
#define AVERAGE 3
#define PERMUTE 4
// Macros defining the value of the ORDER macro, thus defining the precedence of the directions in ADI
// Note: if you change the values of FIRST_IDIR and FIRST_JDIR, you should also change RANDOM_ORDER and PERMUTE_ORDER
#define FIRST_IDIR 0
#define FIRST_JDIR 1
#define RANDOM_ORDER (rand()%2)
// This permutes the order just after a step with odd number
// #define PERMUTE_ORDER ( (g_stepNumber/2)%2 )
// This permutes the order just after a step with even number
#define PERMUTE_ORDER ( ((g_stepNumber+1)/2)%2 )

#if FIRST_JDIR_THEN_IDIR == NO
  #define ORDER FIRST_IDIR
#elif FIRST_JDIR_THEN_IDIR == YES
  #define ORDER FIRST_JDIR
#elif FIRST_JDIR_THEN_IDIR == RANDOM
  #define ORDER RANDOM_ORDER
#elif FIRST_JDIR_THEN_IDIR == AVERAGE
  #ifdef ORDER
    #undef ORDER
  #endif
#elif FIRST_JDIR_THEN_IDIR == PERMUTE
  #define ORDER PERMUTE_ORDER
#else
  #error wrong choice for FIRST_JDIR_THEN_IDIR
#endif

// // For swapping arrays
// #define SWAP_DOUBLE_POINTERS
/**********************/

#if (THERMAL_CONDUCTION==ALTERNATING_DIRECTION_IMPLICIT) && \
    (RESISTIVITY==ALTERNATING_DIRECTION_IMPLICIT)
  #define TDIFF 0
  #define BDIFF 1
  #define NADI  2
#elif THERMAL_CONDUCTION==ALTERNATING_DIRECTION_IMPLICIT
  // [Err] Decomment next 2 lines
  #define TDIFF 0
  #define BDIFF 200 // On purpose a very high number which I will never use
  #define NADI  1
  // [Err] delete next 3 lines
  // #define TDIFF 0
  // #define BDIFF 1
  // #define NADI  2
#elif RESISTIVITY==ALTERNATING_DIRECTION_IMPLICIT
  // [Err] Decomment next 2 lines
  #define TDIFF 200 // On purpose a very high number which I will never use
  #define BDIFF 0
  #define NADI  1
  // [Err] delete next 3 lines
  // #define TDIFF 0
  // #define BDIFF 1
  // #define NADI  2
#else
  /* I still put 1 bc in the type definition (for ease of programming)
  but the kind won't be allocated so it's practically no problem */
  #define NADI  1
  #define BDIFF 200
  #define TDIFF 300
#endif

/******************************************/
/* I define the schemes                   */
#define FRACTIONAL_THETA      1
#define SPLIT_IMPLICIT        2
#define STRANG_LIE            3
#define DOUGLAS_RACHFORD      4
#define PEACEMAN_RACHFORD_MOD 5
#define STRANG                6
/**************************************************/
/* Consistency check of some definitions          */
#if RESISTIVITY==ALTERNATING_DIRECTION_IMPLICIT
  #if METHOD_RES==FRACTIONAL_THETA && !defined(FRACTIONAL_THETA_THETA_RES)
    #error FRACTIONAL_THETA_THETA must be defined when FRACTIONAL_THETA is used
  #endif
  #if METHOD_RES==PEACEMAN_RACHFORD_MOD && !defined(FRACT_RES)
    #error FRACT_RES must be defined when FRACTIONAL_THETA is used
  #endif
#endif
#if THERMAL_CONDUCTION==ALTERNATING_DIRECTION_IMPLICIT
  #if METHOD_TC==FRACTIONAL_THETA && !defined(FRACTIONAL_THETA_THETA_TC)
    #error FRACTIONAL_THETA_THETA must be defined when FRACTIONAL_THETA is used
  #endif
  #if METHOD_TC==PEACEMAN_RACHFORD_MOD && !defined(FRACT_TC)
    #error FRACT_TC must be defined when FRACTIONAL_THETA is used
  #endif
#endif
/***************************************************/

// Time where the diffusion process has arrived (code units)
double extern t_diff;

/* Remember that for the diffusion of B, the BC must refer to the product B*r
as that is the quantity which is advanced by the scheme*/
typedef struct BCS{
  int kind;     /**< Kind of boundary condition: 1=Dirichlet, 2=Hom.Neumann*/
  double values[2];   /**< Values necessary to define the boundary condition
  // (in bc_values[] only element [][0] is used now, in future maybe also [][1], for Robin conditions)*/
} Bcs;

typedef struct LINES{
  int *dom_line_idx;     /**< Indexes (of rows or columns) corresponding to each line*/
  Bcs *lbound[NADI],*rbound[NADI];   /**< Left and right boundary conditions */
  int *lidx, *ridx;      /**< Leftmost and rightmost indexes of the lines. */
  double N;              /**< Number of lines */
} Lines;

// I define a function pointer type, that will take the value of the right bc function
// typedef void (*BoundaryADI) (Lines lines[2], const Data *d, Grid *grid, double t);
typedef void BoundaryADI (Lines lines[2], const Data *d, Grid *grid, double t, int dir);
// I define a function pointer type, that will take the value of the right IJ builder function
typedef void BuildIJ (const Data *d, Grid *grid, Lines *lines, double **Ip, double **Im,
                      double **Jp, double **Jm, double **CI, double **CJ, double **dEdT);

void InitializeLines (Lines *, int);
void GeometryADI (Lines *lines, Grid *grid);
void BoundaryADI_Res(Lines lines[2], const Data *d, Grid *grid, double t, int dir);
void BoundaryADI_TC(Lines lines[2], const Data *d, Grid *grid, double t, int dir);

void FractionalTheta(double **v_new, double **v_old,
                     double **dUres, double **dEdT,
                     const Data *d, Grid *grid,
                     Lines *lines, int diff, int order,
                     double dt, double t0, double theta);

void SplitImplicit(double **v_new, double **v_old,
                  double **dUres, double **dEdT,
                  const Data *d, Grid *grid,
                  Lines *lines, int diff, int order,
                  double dt, double t0, int M);

void PeacemanRachfordMod(double **v_new, double **v_old,
                      double **dUres, double **dEdT,
                      const Data *d, Grid *grid,
                      Lines *lines, int diff, int order,
                      double dt, double t0, double fract, int M);

void DouglasRachford( double **v_new, double **v_old,
                      double **dUres, double **dEdT,
                      const Data *d, Grid *grid,
                      Lines *lines, int diff, int order,
                      double dt, double t0, int M);

void Strang_Lie (double **v_new, double **v_old,
                 double **dUres, double **dEdT,
                 const Data *d, Grid *grid,
                 Lines *lines, int diff, int order,
                 double dt, double t0, int M);
void Strang    (double **v_new, double **v_old,
                double **dUres, double **dEdT,
                const Data *d, Grid *grid,
                Lines *lines, int diff, int order,
                double dt, double t0, int M);

void ExplicitUpdate (double **v, double **b, double **source,
                     double **Hp, double **Hm, double **C,
                     Lines *lines, Bcs *lbound, Bcs *rbound,
                     int compute_inflow, double *inflow, Grid *grid,
                     double dt, int dir);
                     
void ExplicitUpdateDR (double **v, double **b, double **b_der, double **source,
                       double **Hp, double **Hm, double **C,
                       Lines *lines,
                       int compute_inflow, double *inflow,
                       double dt, int dir);
void ApplyBCsonGhosts(double **v, Lines *lines,
                      Bcs *lbound, Bcs *rbound,
                      int dir);
void ImplicitUpdate (double **v, double **b, double **source,
                     double **Hp, double **Hm, double **C,
                     Lines *lines, Bcs *lbound, Bcs *rbound,
                     int compute_inflow, double *inflow, Grid *grid,
                     double dt, int dir);

void tdm_solver(double *x, double const *diagonal, double *up,
                double const *lower, double *rhs, int const N);

double GetCurrADI();

double GetT_old(int j, int i);

#if THERMAL_CONDUCTION  == ALTERNATING_DIRECTION_IMPLICIT
  void BuildIJ_TC (const Data *d, Grid *grid, Lines *lines, double **Ip, double **Im,
                  double **Jp, double **Jm, double **CI, double **CJ, double **dEdT);
  #ifdef TEST_ADI
    void HeatCapacity_test(double *v, double r, double z, double theta, double *dEdT);
  #endif
#endif

#if RESISTIVITY == ALTERNATING_DIRECTION_IMPLICIT
  void BuildIJ_Res (const Data *d, Grid *grid, Lines *lines, double **Ip, double **Im,
                      double **Jp, double **Jm, double **CI, double **CJ, double **useless);
  #if (HAVE_ENERGY && JOULE_EFFECT_AND_MAG_ENG)
    void ResEnergyIncrease(double **dUres, double** Ip_B, double** Im_B, double **Br,
                            Grid *grid, Lines *lines,
                            int compute_inflow, double *inflow,
                            double dt, int dir);
    void ResEnergyIncreaseDR(double **dUres, double** Hp_B, double** Hm_B,
                                           double **Br, double **Br_hat,
                                           Grid *grid, Lines *lines, double dt, int dir);
  #endif
  void ComplainAnisotropic(double *v, double  *eta, double r, double z, double theta);
#endif

/* Stuff to do prim->cons and cons->prim conversions*/
void ConsToPrimLines (Data_Arr U, Data_Arr V, unsigned char ***flag, Lines *lines);
void PrimToConsLines (Data_Arr V, Data_Arr U, Lines *lines);

void SwapDoublePointers (double ***a, double ***b);

#endif
