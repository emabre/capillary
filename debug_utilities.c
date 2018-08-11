#include "pluto.h"
#include "capillary_wall.h"
#include "debug_utilities.h"

/***************************************************************************************
 * [Ema] Print matrix (useful for calling inside gdb) 
 ****************************************************************************************/
void printmat(double **matrix, int dim2, int dim1)
{   
    char coltop[11];
    int i, j;
    for (i = 0; i < dim1; ++i){
      sprintf(coltop, "--%d--", i );
      printf("%21s", coltop);
    }
    printf("\n");
    for (j = 0; j < dim2; ++j) {
      printf("%2d|", j);
        for (i = 0; i < dim1; ++i)
            printf("%21.15g", matrix[j][i]);
        printf(";\n");
    }
    for (i = 0; i < dim1; ++i){
      sprintf(coltop, "--%d--", i );
      printf("%21s", coltop);
    }
    printf("\n");

    // for (i = 0; i < dim1; ++i) {
    //     for (j = 0; j < dim2; ++j)
    //         fprintf(stdout, "%g ", matrix[i][j]);
    //     fprintf(stdout, "\n");
    // }

    // for (i = 0; i < dim1; ++i) {
    //     for (j = 0; j < dim2; ++j)
    //         fprintf(stderr, "%g ", matrix[i][j]);
    //     fprintf(stderr, "\n");
    // }
}

/***************************************************************************************
 * [Ema] Print matrix from 4D vector choosing dimensions to print (useful for calling inside gdb)
 * whichX (X=0,1,2,3) : Integer which tells whether the dimension X has to be printed (whichX = -1)
 *          or kept fixed (whichX = value to keep)
 ****************************************************************************************/
void printmat4d(double ****matrix, int dim2, int dim1, int which0, int which1, int which2, int which3 ) {
    char coltop[11];
    int i, j;
    
    for (i = 0; i < dim1; ++i){
      sprintf(coltop, "--%d--", i );
      printf("%21s", coltop);
    }
      
    printf("\n");
    for (j = 0; j < dim2; ++j) {
      printf("%2d|", j);
      for (i = 0; i < dim1; ++i) {
        if        (which0==-1 && which1==-1) {
          printf("%21.15g", matrix[j][i][which2][which3]);
        } else if (which0==-1 && which2==-1) {
          printf("%21.15g", matrix[j][which1][i][which3]);
        } else if (which0==-1 && which3==-1) {
          printf("%21.15g", matrix[j][which1][which2][i]);
        } else if (which1==-1 && which2==-1) {
          printf("%21.15g", matrix[which0][j][i][which3]);
        } else if (which1==-1 && which3==-1) {
          printf("%21.15g", matrix[which0][j][which2][i]);
        } else if (which2==-1 && which3==-1) {
          printf("%21.15g", matrix[which0][which1][j][i]);
        }
      }
      printf(";\n");
    }

    for (i = 0; i < dim1; ++i){
      sprintf(coltop, "--%d--", i );
      printf("%21s", coltop);
    }
    printf("\n");
}


/************************************************************************************
 * printbox: to print the details of a RBox struct
 * **********************************************************************************/
void printbox(RBox box, char *info) {

  printf("\n---------------------------------------");
  printf("\nBox info: %s", info);

  printf("\nib: %d", box.ib);
  printf("\tie: %d", box.ie);
  printf("\njb: %d", box.jb);
  printf("\tje: %d", box.je);
  printf("\nkb: %d", box.kb);
  printf("\tke: %d", box.ke);
  printf("\nvpos: %d ", box.vpos);
  printf("(CENTER=%d, ", CENTER);
  printf("X1FACE=%d, ", X1FACE);
  printf("X2FACE=%d, ", X2FACE);
  printf("X3FACE=%d)", X3FACE);
  printf("\n---------------------------------------");

  /* I don't print di,dj,dk as they are
     automatically set by the ::BOX_LOOP macro. */
}

/************************************************************************************
 * print_d_corr: to print the details of a RBox struct.
 * whichX (X=0,1,2) : Integer which tells whether the dimension X has to be printed (whichX = -1)
 *          or kept fixed (whichX = value to keep)
 * **********************************************************************************/
void printcorr(Corr corr, char *info) {
  int nv;
  int pp;

  printf("\n--------------------------------------");
  printf("\nCorr info: %s", info);

  printf("\nNpoints: %d", corr.Npoints);

  for (pp=0; pp<corr.Npoints; pp++) {
    printf("\npoint %d; k=%d, j=%d, i=%d",pp, corr.k[pp], corr.j[pp], corr.i[pp]);
    printf("\n Corr Vc vars: (nv, Vc[nv])\n");
    for (nv=0; nv<NVAR; nv++) {
      printf("(%d, %g); ", nv, corr.Vc[nv][pp]);
    }
  }

  printf("\n--------------------------------------");
}

/************************************************************************************
 * DumpQuit: function to perform output of data at the very moment when it is called
 *           and then it quits pluto (to modify it removing the quit of pluto you
 *           have to figure out how to deal with the file names and so on. In this
 *           function I spoil the Runtime *ini structure, thus THE PLUTO PROGRAM MUST NOT CONTINUE)
 *    (written copy-modifying PLUTO's CheckForOutput() )
 * **********************************************************************************/
void DumpQuit (const Data *d, Runtime *ini, Grid *grid) {
  static int first_call = 1;
  int  n;
  Output *output;
  // static time_t clock_beg[MAX_OUTPUT_TYPES], clock_end;

    print1("\n[DumpQuit] I write data (if any) and quit\n");

  // /* -- on first execution initialize
  //     current beginning time for all output types -- */

  // if (first_call){
  //   for (n = 0; n < MAX_OUTPUT_TYPES; n++) time(clock_beg + n);
  // }

  // /* -- get current time -- */
  // time(&clock_end);

  /* -------------------------------------------------------
          start main loop on outputs
    ------------------------------------------------------- */
  for (n = 0; n < MAX_OUTPUT_TYPES; n++){
    output = ini->output + n;
   // I modify the nfile field to disinguish this output from all the others 
    output->nfile = 9998;
    
    if (output->dt > 0.0 || output->dn > 0){
      WriteData(d, output, grid);
    }
  }

  first_call = 0;

  QUIT_PLUTO(1);
}
