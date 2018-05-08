#include "debug_utilities.h"
#include "pluto.h"
/***************************************************************************************
 * [Ema] Print matrix (useful for calling inside gdb) 
 ****************************************************************************************/
void printMat(double **matrix, int dim1, int dim2)
{
    int i, j;
    for (i = 0; i < dim1; ++i)
    {
        for (j = 0; j < dim2; ++j)
            print1("%d ", matrix[i][j]);
        print1("\n");
    }
}