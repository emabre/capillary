#ifndef TABLE_UTILITIES_H
#define TABLE_UTILITIES_H

int ReadASCIITableSettings(const char* table_finame, double *xmax, double *xmin, int *Nx, double *ymax, double *ymin, int *Ny);
int ReadASCIITableMatrix(const char* table_finame, double **f, int Nx, int Ny);
void ReprintTable(Table2D *tab, const char *tabname);

#endif