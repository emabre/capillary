#include "pluto.h"
#include "adi.h"

#if KBEG != KEND
  #error grid in k direction should only be of 1 point
#endif

void ConsToPrimLines (Data_Arr U, Data_Arr V, unsigned char ***flag, Lines *lines)
/*!
 *  Convert conservative variables \c U to
 *  an primitive variables \c V, along an array of lines ( a Lines struct)
 *  Note that <tt>[nv]</tt> is the fastest running index for \c U 
 *  while it is the slowest running index for \c V.
 *********************************************************************** */
{
  int   i, j, k, nv, err;
  int   l, Nlines;
  int   ibeg, iend;
  int   current_dir;
  static double **v;

  if (v == NULL){
    v = ARRAY_2D(NMAX_POINT, NVAR, double);
  }

/* ----------------------------------------------
    Save current sweep direction and by default,
    perform the conversion along X1 stripes
   ---------------------------------------------- */

  current_dir = g_dir; /* save current direction */
  g_dir = IDIR;
  
/* -----------------------------------------------
    Do the actual conversion
   ------------------------------------------------ */
  Nlines = lines[IDIR].N;
  KDOM_LOOP (k) {
    g_k = k;
    for (l = 0; l < Nlines; l++) {
      g_j = j = lines[IDIR].dom_line_idx[l];
      ibeg = lines[IDIR].lidx[l];
      iend = lines[IDIR].ridx[l];

      err = ConsToPrim (U[k][j], v, ibeg, iend, flag[k][j]);

      if (err) {
        print1("[ConsToPrimLines] Error converting Cons->Prim (k:%d,j:%d)", k,j);
        // QUIT_PLUTO(1);
        print1("\nI move on...\n");        
      }

      for (i = ibeg; i <= iend; i++) NVAR_LOOP(nv) V[nv][k][j][i] = v[i][nv];
    }
  }
  g_dir = current_dir; /* restore current direction */

}
/* ********************************************************************* */
void PrimToConsLines (Data_Arr V, Data_Arr U, Lines *lines)
/*!
 *  Convert a primitive variables \c V  to
 *  conservative variables \c U, along an array of lines (a Lines structs)
 *  Note that <tt>[nv]</tt> is the fastest running index for \c U 
 *  while it is the slowest running index for \c V.
 *
 *********************************************************************** */
{
  int   i, j, k, nv;
  int   ibeg, iend;
  int   l, Nlines;
  int   current_dir;
  static double **v, **u;

  if (v == NULL) {
    v = ARRAY_2D(NMAX_POINT, NVAR, double);
    u = ARRAY_2D(NMAX_POINT, NVAR, double);
  }

  current_dir = g_dir; /* save current direction */
  g_dir = IDIR;

/* -----------------------------------------------
    Do the actual conversion
   ------------------------------------------------ */
  Nlines = lines[IDIR].N;
  KDOM_LOOP (k) {
    for (l = 0; l < Nlines; l++) {
      j = lines[IDIR].dom_line_idx[l];
      ibeg = lines[IDIR].lidx[l];
      iend = lines[IDIR].ridx[l];

      for (i = ibeg; i <= iend; i++) VAR_LOOP(nv) v[i][nv] = V[nv][k][j][i];

      PrimToCons (v, U[k][j], ibeg, iend);
    }
  }
  g_dir = current_dir; /* restore current direction */

}
