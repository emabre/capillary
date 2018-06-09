#ifndef ADI_H
#define ADI_H
/* Contains some structures used for the adi method*/

// #define TEST_ADI

/*********************
 * Some useful macros
 * *******************/
#define DIRICHLET    1
#define NEUMANN_HOM  2

/* To loop on lines, the correct use is, usually:
   LINES_LOOP(lines[IDIR], l, j, i)
   or:
   LINES_LOOP(lines[JDIR], l, i, j)*/
#define LINES_LOOP(lines, l, line_idx , line_sweeper) \
  for ((line_idx)=lines.dom_line_idx[(l)=0]; (l)<lines.N; (line_idx)=lines.dom_line_idx[++(l)]) \
  for ((line_sweeper)=lines.lidx[(l)]; (line_sweeper)<=lines.ridx[(l)]; (line_sweeper)++)
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
#endif

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

void InitializeLines (Lines *, int);
void GeometryADI (Lines *lines, Grid *grid);
void BoundaryRes_ADI(Lines lines[2], const Data *d, Grid *grid, double t);
void BoundaryTC_ADI(Lines lines[2], const Data *d, Grid *grid, double t);

void ExplicitUpdate (double **v, double **rhs, double **source,
                     double **Hp, double **Hm, double **C,
                     Lines *lines, Bcs *lbound, Bcs *rbound, double dt, int const dir);
void ImplicitUpdate (double **v, double **rhs, double **source,
                     double **Hp, double **Hm, double **C,
                     Lines *lines, Bcs *lbound, Bcs *rbound, double dt, int const dir);
void tdm_solver(double *x, double const *diagonal, double const *upper,
                double const *lower, double const *right_hand_side, int const N);
                     
#if THERMAL_CONDUCTION == ALTERNATING_DIRECTION_IMPLICIT
void BuildIJ_TC (const Data *d, Grid *grid, Lines *lines, double **Ip, double **Im,
                    double **Jp, double **Jm, double **CI, double **CJ, double **dEdT);
#ifdef TEST_ADI
  void HeatCapacity_test(double *v, double r, double z, double theta, double *dEdT);
#endif
#endif

#if RESISTIVITY == ALTERNATING_DIRECTION_IMPLICIT
void BuildIJ_Res (const Data *d, Grid *grid, Lines *lines, double **Ip, double **Im,
                     double **Jp, double **Jm, double **CI, double **CJ);
  #if (HAVE_ENERGY && JOULE_EFFECT_AND_MAG_ENG)
    void ResEnergyIncrease(double **dUres, double** Ip_B, double** Im_B, double **Br,
                           Grid *grid, Lines *lines, double dt, int dir);
  #endif
#endif

/* Stuff to do prim->cons and cons->prim conversions*/
void ConsToPrimLines (Data_Arr U, Data_Arr V, unsigned char ***flag, Lines *lines);
void PrimToConsLines (Data_Arr V, Data_Arr U, Lines *lines);

#endif