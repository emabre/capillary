#include "debug_utilities.h"
#include "pluto.h"
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
        printf("|\n");
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
 * which[X] : Integer which tells whether the dimension X has to be printed (which[X] = -1)
 *          or kept fixed (which[X] = value to keep)
 ****************************************************************************************/
void printmat4d(double ****matrix, int dim2, int dim1, int which[4]) {
    int i, j;

    for (i = 0; i < dim1; ++i)
      printf("%11s", "---");
      
    printf("\n");
    for (j = 0; j < dim2; ++j) {
      printf("|");
      for (i = 0; i < dim1; ++i) {
        if        (which[0]==-1 && which[1]==-1) {
          printf("%11.5g", matrix[j][i][which[2]][which[3]]);
        } else if (which[0]==-1 && which[2]==-1) {
          printf("%11.5g", matrix[j][which[1]][i][which[3]]);
        } else if (which[0]==-1 && which[3]==-1) {
          printf("%11.5g", matrix[j][which[1]][which[2]][i]);
        } else if (which[1]==-1 && which[2]==-1) {
          printf("%11.5g", matrix[which[0]][j][i][which[3]]);
        } else if (which[1]==-1 && which[3]==-1) {
          printf("%11.5g", matrix[which[0]][j][which[2]][i]);
        } else if (which[2]==-1 && which[3]==-1) {
          printf("%11.5g", matrix[which[0]][which[1]][j][i]);
        }
      }
      printf("|\n");
    }

    for (i = 0; i < dim1; ++i)
      printf("%11s", "---");

    printf("\n");
}