#ifndef GAMMA_TRANSP_H
#define GAMMA_TRANSP_H

#define MULTIPLE_GHOSTS YES

/* ********************************************************************* */
/*! [Ema]The Corr structure contains the correction to the solution 3D array
 to apply when advancing different directions (useful to have multiple ghost cells
located on the same internal cell)
   ********************************************************************* */
typedef struct CORR{
  double **Vc;  /**< The main three-index data array used for cell-centered
                       primitive variables. The index order is
                       <tt>Vc[nv][pp]</tt> where \c nv gives the variable
                       index while \c pp \c is the index counting the cell where
                       one employes the correction (Corr struct resembles a sparse
                     matrix).*/
  int *i; /* k,j and i arrays are the
  locations of the cells to correct in the x_3,
  x_2 and x_1 direction. Which means that, for each pp you have to correct
  location i[pp],j[pp],k[pp] by overwriting Vs[nv][pp]
  to d-Vs[nv][k][j][i] for every nv*/
  int *j;
  int *k;
  int Npoints; // number of points to be corrected (pp will be always <=Npoints-1)
} Corr;


#define Q_ELEM 4.8032e-10 /*Elementary (electron) charge in statcoulomb*/
/*Epsilon 0 constant in Esaulov's paper Plasma Phys control Fusion 43 571 */
#define EPS_0_ESAU 5.81e-11 /* in erg*/

double Gamma_jil[2][10][4];
int extern idx_rcap, idx_zcap, idx_start_electr;

// For correcting the internal boundary and putting a double ghost on the wall
// corner cell at the capillary exit
Corr extern d_correction[3];

double Gamma_i(int i, double xe, double w);
double Gamma_ij(int i, int j, double w);
double cl_ei(double ne, double T);
double cl_en(double ioniz, double T);
// double cl_ee(double ne, double T);

int find_idx_closest(double *vec, int Nvec, double v);

void alloc_Data(Data *data);
void copy_Data_Vc(Data *d_target, const Data *d_source);

#endif
