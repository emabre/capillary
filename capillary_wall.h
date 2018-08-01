#ifndef CAPILLARY_WALL_H
#define CAPILLARY_WALL_H

#define MULTIPLE_GHOSTS YES

/* Macros to refer to the indexes of RBox rbox_center_capWall[]
   depending on the capillary wall(boundary) region you need
   IF YOU ADD A NEW ONE: increase the size of RBox rbox_center_capWall[]
   inside capillary_wall.c*/
#define CAP_WALL_INTERNAL         0
#define CAP_WALL_EXTERNAL         1
#define CAP_WALL_CORNER_INTERNAL  2
#define CAP_WALL_CORNER_EXTERNAL  3

/* Variables defined for keeping information on the geometry of the capillary:
*/
double extern const zcap, dzcap, rcap;
// Actual values used inside the simulation for zcap, rcap, dzcap;
double extern zcap_real, rcap_real, dzcap_real;
int extern capillary_not_set;
int extern i_cap_inter_end, j_cap_inter_end, j_elec_start;

/* Variables defined for computation of energy conservation:
  en_res_loss : energy lost by resistivity (resistive part of poynting flux
                through the boundary)
  en_cond_loss : energy lost by conduction through boundary
  en_Bvloss : energy lost due to term B(B*v) through boundary
              in equation 6.4 of userguide
  en_advloss: energy lost by advection of the total energy rhough boundary*/
double extern en_cond_loss, en_res_loss, en_adv_loss, en_Bvloss;

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

//int extern idx_rcap, idx_zcap, idx_start_electr;

// For correcting the internal boundary and putting a double ghost on the wall
// corner cell at the capillary exit
Corr extern d_correction[3];

int SetRemarkableIdxs(Grid *grid);
int FindIdxClosest(double *vec, int Nvec, double v);

// void alloc_Data(Data *data);
Data* alloc_Data();
void copy_Data_Vc(Data *d_target, const Data *d_source);
void free_Data(Data *data);

// RBox *GetRBoxCap(int side, int vpos);
void SetRBox_capWall(int Nghost);
void ReflectiveBoundCap (double ****q, int nv, int s, int side, int vpos);
void ZeroBoundCap (double ****q, int nv, int s, int side, int vpos);
void SetNotEvolvedVar (int nv);

#endif
