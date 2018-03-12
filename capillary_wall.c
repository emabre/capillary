#include <math.h>
#include "pluto.h"
#include "capillary_wall.h"

//debug macro
// #define DBG_FIND_CLOSEST

int idx_rcap = 0;
int idx_zcap = 0;
int idx_start_electr = 0;

Corr d_correction[3] = { {},{},{} };

/*Finds the index of the element in vec closest to the value of v*/
int find_idx_closest(double *vec, int Nvec, double v){
  int i, i_mindiff;
  double diff;

  #ifdef DBG_FIND_CLOSEST
    print1("\nvec[0]:%g, v:%g", vec[0], v);
    print1("\nfabs(vec[0]-v):%g", fabs(vec[0]-v));
  #endif

  diff = fabs(vec[0]-v);
  for (i=1;i<Nvec;i++){
    if (fabs(vec[i]-v) < diff) {
      #ifdef DBG_FIND_CLOSEST
        print1("\nabs(vec[%d](=%e)-%e)<%e",i,vec[i],v,diff);
      #endif
      diff = fabs(vec[i]-v);
      i_mindiff = i;
    } else {
      #ifdef DBG_FIND_CLOSEST
        print1("\nfabs(vec[%d](=%e)-%e)>%e",i,vec[i],v,diff);
      #endif
    }
  }
  return i_mindiff;
}

/*Allocate space for a Data element*/
void alloc_Data(Data *data) {
  // print1 ("\n> Memory allocation\n");
  data->Vc = ARRAY_4D(NVAR, NX3_TOT, NX2_TOT, NX1_TOT, double);
  data->Uc = ARRAY_4D(NX3_TOT, NX2_TOT, NX1_TOT, NVAR, double);

  #ifdef STAGGERED_MHD
   data->Vs = ARRAY_1D(DIMENSIONS, double ***);
   D_EXPAND(
     data->Vs[BX1s] = ArrayBox( 0, NX3_TOT-1, 0, NX2_TOT-1,-1, NX1_TOT-1); ,
     data->Vs[BX2s] = ArrayBox( 0, NX3_TOT-1,-1, NX2_TOT-1, 0, NX1_TOT-1); ,
     data->Vs[BX3s] = ArrayBox(-1, NX3_TOT-1, 0, NX2_TOT-1, 0, NX1_TOT-1);)
  #endif

  #if UPDATE_VECTOR_POTENTIAL == YES
   D_EXPAND(                                                  ,
     data->Ax3 = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);  ,
     data->Ax1 = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
     data->Ax2 = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
   )
  #endif

  #if RESISTIVITY != NO
   data->J = ARRAY_4D(3,NX3_TOT, NX2_TOT, NX1_TOT, double);
  #endif

  data->flag = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, unsigned char);
}

/*Copies the Vc field of a Data structure*/
void copy_Data_Vc(Data *d_target, const Data *d_source) {
  int i,j,k,nv;
  TOT_LOOP (k,j,i) {
    VAR_LOOP (nv) d_target->Vc[nv][k][j][i] = d_source->Vc[nv][k][j][i];
  }
}
