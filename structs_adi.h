/* Contains some structures used for the adi method*/

#define DIRICHLET    1
#define NEUMANN_HOM  2

typedef struct BCS{
  int bc_kind;     /**< Kind of boundary condition*/
  double bc_values[2];   /**< Values necessary to define the boundary condition*/
} Bcs;

typedef struct LINES{
  int *dom_line_idx;     /**< Indexes (of rows or columns) corresponding to each line*/
  Bcs *lbound,*rbound;   /**< Left and right boundary conditions */
  int *lidx, *ridx;      /**< Leftmost and rightmost indexes of the lines. */
  double N;              /**< Number of lines */
} Lines;

void initialize_Lines(Lines *, int);
int CountLines(Data *d, Grid *grid, int dir);