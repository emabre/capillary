/* Contains some structures used for the adi method*/

#define DIRICHLET    1
#define NEUMANN_HOM  2

#if (THERMAL_CONDUCTION==ALTERNATING_DIRECTION_IMPLICIT) && \
    (RESISTIVITY==ALTERNATING_DIRECTION_IMPLICIT)
  #define TDIFF 0
  #define BDIFF 1
  #define NADI  2
#elif THERMAL_CONDUCTION==ALTERNATING_DIRECTION_IMPLICIT
  #define TDIFF 0
  #define NADI  1
#elif RESISTIVITY==ALTERNATING_DIRECTION_IMPLICIT
  #define BDIFF 0
  #define NADI  1
#else
  /* I still put 1 bc in the type definition (for ease of programming)
  but the bc_kind won't be allocated so it's practically no problem */
  #define NADI  1
#endif

typedef struct BCS{
  int bc_kind;     /**< Kind of boundary condition: 1=Dirichlet, 2=Hom.Neumann*/
  double bc_values[NADI][2];   /**< Values necessary to define the boundary condition
  // (in bc_values[] only element [][0] is used now, in future maybe also [][1], for Robin conditions)*/
} Bcs;

typedef struct LINES{
  int *dom_line_idx;     /**< Indexes (of rows or columns) corresponding to each line*/
  Bcs *lbound,*rbound;   /**< Left and right boundary conditions */
  int *lidx, *ridx;      /**< Leftmost and rightmost indexes of the lines. */
  double N;              /**< Number of lines */
} Lines;

void InitializeLines(Lines *, int);
void GeometryADI(Lines *lines, Grid *grid);