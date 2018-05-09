#include "debug_utilities.h"
#include "pluto.h"
/***************************************************************************************
 * [Ema] Print matrix (useful for calling inside gdb) 
 ****************************************************************************************/
void printmat(double **matrix, int dim2, int dim1)
{
    int i, j;
    for (i = 0; i < dim1; ++i)
      printf("%11s", "---");
    printf("\n");
    for (j = 0; j < dim2; ++j) {
      printf("|");
        for (i = 0; i < dim1; ++i)
            printf("%11.5g", matrix[j][i]);
        printf("|\n");
    }
    for (i = 0; i < dim1; ++i)
      printf("%11s", "---");
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