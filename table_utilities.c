#include "pluto.h"
// #include<math.h>

/* Function to read only the first lines (containing settings) of a file containing ascii table of some quantity*/
int ReadASCIITableSettings(const char* table_finame, int *logspacing, double *xmin, double *xmax, int *Nx, double *ymin, double *ymax, int *Ny) {
  FILE *table;

  table = fopen(table_finame, "r");

  // Read the first 2 comment lines
  fscanf(table,"%*[^\n] \n");
  fscanf(table,"%*[^\n] \n");

  // Reading settings
  fscanf(table, "%d\n", logspacing);
  fscanf(table, "%lf", xmin); // xmax
  fscanf(table, "%lf", xmax); // xmin
  fscanf(table, "%d", Nx); // xmin
  fscanf(table, "%lf", ymin); // xmin
  fscanf(table, "%lf", ymax); // xmin
  fscanf(table, "%d", Ny); // xmin

  fclose(table);
  return 0;
}

/*Function to read the part of a file containing an ascii table of some quanity and store it in a matrix (f)*/
int ReadASCIITableMatrix(const char* table_finame, double **f, int Nx, int Ny) {
  FILE *table;
  // char headers[300];
  int i,j;

  table = fopen(table_finame, "r");

  // Skip settings part
  for (i=0;i<8;i++) {
    fscanf(table, "%*[^\n]\n");
    // fscanf(table, "%[^\n]\n", headers);
  }

  for (i=0; i<Nx; i++) {
    for (j=0; j<Ny; j++) {
      fscanf(table, "%lf", &(f[j][i]));
      // printf("\nJust read: %g", f[i][j]);
    }
  }

  fclose(table);
  return 0;
}

/* For testing purposes */
void ReprintTable(Table2D *tab, const char *tabname) {
  int i,j;

  print1("\n-----------------------------------------------");
  print1("\nTable '%s' is:", tabname);

  print1("\nmax T=%e", pow(10, tab->lnxmin));
  print1("\nmin T=%e",pow(10, tab->lnxmax));
  print1("\nN T=%d",tab->nx);
  print1("\nmin rho=%e",pow(10, tab->lnymin));
  print1("\nmax rho=%e",pow(10, tab->lnymax));
  print1("\nN rho=%d",tab->ny);

  print1("\n%3s%10s", "*","*");
  for (j = 0; j <tab->ny; j++){
    // sprintf(coltop, "%d", j );
    // print1("%10s ", coltop);
    // sprintf(coltop, "%d", j );
    print1("%10d", j);
  }
  print1("\n%3s%10s", "*","*");
  for (j = 0; j <tab->ny; j++){
    print1("%10.3e", tab->y[j]);
  }
  print1("\n");
  for(i=0; i<tab->nx; i++) {
    print1("%3d", i);
    print1("%10.3e", tab->x[i]);
    for (j=0; j<tab->ny; j++) {
      print1("%10.3e", tab->f[j][i]);
    }
    print1("\n");
  }
  print1("-----------------------------------------------\n");

}