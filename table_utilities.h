#ifndef TABLE_UTILITIES_H
#define TABLE_UTILITIES_H

int ReadASCIITableSettings(const char* table_finame, double *xmin, double *xmax, int *Nx, double *ymin, double *ymax, int *Ny);
int ReadASCIITableMatrix(const char* table_finame, double **f, int Nx, int Ny);
void ReprintTable(Table2D *tab, const char *tabname);

#endif