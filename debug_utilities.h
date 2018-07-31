#ifndef DEBUG_UTILITIES
#define DEBUG_UTILITIES

void printmat(double **matrix, int dim1, int dim2);
void printmat4d(double ****matrix, int dim2, int dim1, int which0, int which1, int which2, int which3 );

void printbox(RBox box, char *info);

#endif