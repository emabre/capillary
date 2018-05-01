/* Contains some structures used for the adi method*/

#define DIRICHLET    1
#define NEUMANN_HOM  2

typedef struct BCS{
  int bc_kind;     /**< Kind of boundary condition: 1=Dirichlet, 2=Hom.Neumann*/
  double bc_values[2];   /**< Values necessary to define the boundary condition
  // (in bc_values[] only element 0 is used now, in future maybe also 1, for Robin conditions)*/
} Bcs;

typedef struct LINES{
  int *dom_line_idx;     /**< Indexes (of rows or columns) corresponding to each line*/
  Bcs *lbound,*rbound;   /**< Left and right boundary conditions */
  int *lidx, *ridx;      /**< Leftmost and rightmost indexes of the lines. */
  double N;              /**< Number of lines */
} Lines;

void initialize_Lines(Lines *, int);
void GeometryADI(Lines *lines, Grid *grid);