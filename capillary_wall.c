#include <math.h>
#include "pluto.h"
#include "capillary_wall.h"
#include "debug_utilities.h"

//debug macro
// #define DBG_FIND_CLOSEST

// Box useful for setting bcs internal to the domain
static RBox rbox_center_capWall[2];
static RBox rbox_center_capCorn[1];

double const zcap = ZCAP/UNIT_LENGTH;
double const dzcap = DZCAP/UNIT_LENGTH;
double const rcap = RCAP/UNIT_LENGTH;
// Actual values used inside the simulation for zcap, rcap, dzcap;
double zcap_real, rcap_real, dzcap_real;
int capillary_not_set = 1;
int i_cap_inter_end, j_cap_inter_end, j_elec_start;

double en_tc_in = 0;
double en_adv_in = 0;
double en_res_in = 0;

Corr d_correction[3] = { {},{},{} };

int not_allocated_d_correction = 1;

/******************************************************************/
/* Sets the remarkable indexes of the grid for the capillary*/
/******************************************************************/
int SetRemarkableIdxs(Grid *grid){

  /* Capillary:
                                     j=j_elec_start (first cell belonging to electrode)
                                      :     j=j_cap_inter_end (ghost)
               r                      :      |
               ^    |                 :      *
               |    |       wall      :      *
                    |                 v      *
i=i_cap_inter_end+1 |****(ghosts)*****o*******|
i=i_cap_inter_end   |                        j=j_cap_inter_end+1 (first outside, not ghost)
                    |
                    |
i=0                 o-------------------------------->(axis)
                   j=0           -> z
    // I should not change the grid size exacly on the capillary end!
    */

  /* I find the indexes of the cells closest to the capillary bounds*/
  i_cap_inter_end = IBEG + FindIdxClosest(&(grid[IDIR].xr_glob[IBEG]), IEND-IBEG+1, rcap);
  j_cap_inter_end = JBEG + FindIdxClosest(&(grid[JDIR].xr_glob[JBEG]), JEND-JBEG+1, zcap);
  j_elec_start = JBEG + FindIdxClosest(&(grid[JDIR].xl_glob[JBEG]), JEND-JBEG+1, zcap-dzcap);

  if (j_elec_start > j_cap_inter_end) {
    print1("\n[SetRemarkableIdxs]Electrode appears to start after end of capillary! Quitting.");
    QUIT_PLUTO(1);
  }

  rcap_real = grid[IDIR].xr_glob[i_cap_inter_end];
  zcap_real = grid[JDIR].xr_glob[j_cap_inter_end];
  dzcap_real = grid[JDIR].xr_glob[j_cap_inter_end]-grid[JDIR].xl_glob[j_elec_start];

  print1("\n  -------------------------------------------------------------------------");
  print1("\n  Indexes of remarkable internal bounary points:");
  print1("\n  i_cap_inter_end: \t%d", i_cap_inter_end);
  print1("\n  j_cap_inter_end: \t%d", j_cap_inter_end);
  print1("\n  j_elec_start:    \t%d\n", j_elec_start);
  print1("\n  Remarkable points:");
  print1("\n  Capillary radius,      set: %g; \tactual: %g \t(cm)", RCAP, rcap_real*UNIT_LENGTH);
  print1("\n  Capillary half length, set: %g; \tactual: %g \t(cm)", ZCAP, zcap_real*UNIT_LENGTH);
  print1("\n  Electrode length,      set: %g; \tactual: %g \t(cm)", DZCAP, dzcap_real*UNIT_LENGTH);
  print1("\n  ( electrode actual start: z=%g; \t(cm) )",(zcap_real-dzcap_real)*UNIT_LENGTH);
  print1("\n");
  print1("\n  Just so that you know:");
  print1("\n  NX3_TOT=%d, NX2_TOT=%d, NX1_TOT=%d",NX3_TOT, NX2_TOT, NX1_TOT);
  print1("\n  IBEG=%d, IEND=%d, JBEG=%d, JEND=%d, KBEG=%d, KEND=%d", IBEG, IEND, JBEG, JEND, KBEG, KEND);
  print1("\n  ---------------------------------------------------------------------------");
  print1("\n");
  
  capillary_not_set = 0;
  return 0;
}

/******************************************************************/
/*Finds the index of the element in vec closest to the value of v*/
/******************************************************************/
int FindIdxClosest(double *vec, int Nvec, double v){
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

/*--------------------------------------------------------------------
 alloc_Data(Data *data) : Allocate space for a Data element         
--------------------------------------------------------------------*/
Data* alloc_Data() {
  Data* newdata = (Data *)malloc(sizeof(Data));
  /*******************************************************************
   *  What follows has been copied from initialize.c, lines: 448-472 *
   *******************************************************************/
  // print1 ("\n> Memory allocation\n");
  newdata->Vc = ARRAY_4D(NVAR, NX3_TOT, NX2_TOT, NX1_TOT, double);
  newdata->Uc = ARRAY_4D(NX3_TOT, NX2_TOT, NX1_TOT, NVAR, double);

  #ifdef STAGGERED_MHD
   newdata->Vs = ARRAY_1D(DIMENSIONS, double ***);
   D_EXPAND(
     newdata->Vs[BX1s] = ArrayBox( 0, NX3_TOT-1, 0, NX2_TOT-1,-1, NX1_TOT-1); ,
     newdata->Vs[BX2s] = ArrayBox( 0, NX3_TOT-1,-1, NX2_TOT-1, 0, NX1_TOT-1); ,
     newdata->Vs[BX3s] = ArrayBox(-1, NX3_TOT-1, 0, NX2_TOT-1, 0, NX1_TOT-1);)
  #endif

  #if UPDATE_VECTOR_POTENTIAL == YES
   D_EXPAND(                                                  ,
     newdata->Ax3 = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);  ,
     newdata->Ax1 = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
     newdata->Ax2 = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
   )
  #endif

  #if RESISTIVITY != NO
   newdata->J = ARRAY_4D(3,NX3_TOT, NX2_TOT, NX1_TOT, double);
  #endif

  newdata->flag = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, unsigned char);

  return newdata;
}

/*--------------------------------------------------------------------
 Copies the Vc field of a Data structure
--------------------------------------------------------------------*/
void copy_Data_Vc(Data *d_target, const Data *d_source) {
  int i,j,k,nv;
  TOT_LOOP (k,j,i) {
    VAR_LOOP (nv) d_target->Vc[nv][k][j][i] = d_source->Vc[nv][k][j][i];
  }
}

/*--------------------------------------------------------------------
 free_Data(Data *data) Free the memory of a Data element         
--------------------------------------------------------------------*/
void free_Data(Data *data) {
  // print1 ("\n> Memory allocation\n");
  FreeArray4D ((void *) data->Vc);
  FreeArray4D ((void *) data->Uc);

  #ifdef STAGGERED_MHD
    #error deallocation of Vs is not yet implemented 
  #endif

  #if UPDATE_VECTOR_POTENTIAL == YES
    #error deallocation of Ax1, Ax2, Ax3 is not yet implemented
  #endif

  #if RESISTIVITY != NO
    FreeArray4D ((void *) data->J);
  #endif

  FreeArray3D ((void *) data->flag);

  free((Data *) data);
}

void SetRBox_capWall(int Nghost) {
  int s;
  /* ---------------------------------------------------
    0. set CAP_WALL_INTERNAL grid index ranges
   --------------------------------------------------- */
  
  s = CAP_WALL_INTERNAL;

  rbox_center_capWall[s].vpos = CENTER;

  rbox_center_capWall[s].ib = i_cap_inter_end + 1;
  rbox_center_capWall[s].ie = i_cap_inter_end + 1 + Nghost - 1;
  
  rbox_center_capWall[s].jb = Nghost;
  rbox_center_capWall[s].je = j_cap_inter_end - Nghost;
  
  rbox_center_capWall[s].kb = 0;
  rbox_center_capWall[s].ke = NX3_TOT-1;

  /* ---------------------------------------------------
    2. set CAP_WALL_EXTERNAL grid index ranges
   --------------------------------------------------- */

  s = CAP_WALL_EXTERNAL;

  rbox_center_capWall[s].vpos = CENTER;

  rbox_center_capWall[s].ib = i_cap_inter_end + 1 + Nghost;
  rbox_center_capWall[s].ie = NX1_TOT - 1 - Nghost;

  rbox_center_capWall[s].jb = j_cap_inter_end - Nghost + 1;
  rbox_center_capWall[s].je = j_cap_inter_end;

  rbox_center_capWall[s].kb = 0;
  rbox_center_capWall[s].ke = NX3_TOT-1;

  /* ---------------------------------------------------
    3. set grid index ranges for capillary
       corner correction (rbox_center_capCorn),
   --------------------------------------------------- */
  s = 0;
  
  rbox_center_capCorn[s].vpos = CENTER;

  rbox_center_capCorn[s].ib = i_cap_inter_end + 1;
  rbox_center_capCorn[s].ie = i_cap_inter_end + 1 + Nghost - 1;

  rbox_center_capCorn[s].jb = j_cap_inter_end - Nghost + 1;
  rbox_center_capCorn[s].je = j_cap_inter_end;

  rbox_center_capCorn[s].kb = 0;
  rbox_center_capCorn[s].ke = NX3_TOT-1;

  #ifdef DEBUG_BCS
    printbox(rbox_center_capWall[CAP_WALL_INTERNAL], "rbox_center_capWall[CAP_WALL_INTERNAL]");
    printbox(rbox_center_capWall[CAP_WALL_EXTERNAL], "rbox_center_capWall[CAP_WALL_EXTERNAL]");
    printbox(rbox_center_capCorn[0], "rbox_center_capCorn[0]");
  #endif
}

/* ********************************************************************* */
RBox *GetRBoxCap(int side, int vpos)
/*!
 *  Returns a pointer to a local static RBox 
 *
 *  \param[in]  side  the region of the computational domain where 
 *                    the box is required.
 *                    Possible values :
 *                    CAP_WALL_INTERNAL       
 *                    CAP_WALL_EXTERNAL       
 *                    CAP_WALL_CORNER_INTERNAL
 *                    CAP_WALL_CORNER_EXTERNAL
 * 
 *  \param[in]  vpos  the variable position inside the cell:
 *                    CENTER is the only avaiable for now.
 *
 *********************************************************************** */
{
  if (vpos != CENTER) {
    print1("\n[GetRBoxCap] Only vpos == CENTER is implemented!");
  }
  if      (side == CAP_WALL_INTERNAL) 
    return &(rbox_center_capWall[CAP_WALL_INTERNAL]);
  else if (side == CAP_WALL_EXTERNAL)
    return &(rbox_center_capWall[CAP_WALL_EXTERNAL]);
  else if (side == CAP_WALL_CORNER_INTERNAL)
    return &(rbox_center_capCorn[0]);
  else if (side == CAP_WALL_CORNER_EXTERNAL)
    return &(rbox_center_capCorn[0]);
  else
    return NULL;
}

void ReflectiveBoundCap (double ****q, int nv, int s, int side, int vpos)
/*!
 * [Created by Ema]
 * Make symmetric (s = 1) or anti-symmetric (s=-1) profiles.
 * The sign is set by the FlipSign() function. 
 *
 * \param [in,out] q   a 3D flow quantity
 * \param [in]    nv    kind of variable to apply the bc, e.g.: RHO, iBPHI, PRS...
 * \param [in] s   an integer taking only the values +1 (symmetric 
 *                 profile) or -1 (antisymmetric profile)
 *   
 *********************************************************************** */
{
  int   i, j, k, pp;
  RBox *box = GetRBoxCap(side, vpos);

  if (not_allocated_d_correction && (side == CAP_WALL_CORNER_INTERNAL || side == CAP_WALL_CORNER_EXTERNAL)) {
    // I allocate memory for d_correction
    // [Rob] I could define a function to do this in two lines
    // First I count how many points are needed
    pp = 0;
    BOX_LOOP(box, k, j, i) pp++;

    d_correction[IDIR].Npoints = pp;
    d_correction[IDIR].i = ARRAY_1D(d_correction[IDIR].Npoints, int);
    d_correction[IDIR].j = ARRAY_1D(d_correction[IDIR].Npoints, int);
    d_correction[IDIR].k = ARRAY_1D(d_correction[IDIR].Npoints, int);
    d_correction[IDIR].Vc = ARRAY_2D( NVAR, d_correction[IDIR].Npoints, double);

    d_correction[JDIR].Npoints = pp;
    d_correction[JDIR].i = ARRAY_1D(d_correction[JDIR].Npoints, int);
    d_correction[JDIR].j = ARRAY_1D(d_correction[JDIR].Npoints, int);
    d_correction[JDIR].k = ARRAY_1D(d_correction[JDIR].Npoints, int);
    d_correction[JDIR].Vc = ARRAY_2D( NVAR, d_correction[JDIR].Npoints, double);

    not_allocated_d_correction = 0;
  }

  if (side == CAP_WALL_INTERNAL) {
    /* [Ema] Values are simply reflected across the boundary
          (depending on "s", with sign changed or not!),
          even if the ghost cells are more than one,
          e.g.: with 2 ghosts cells per side:
                q[nv][k][j][IEND+1] = q[nv][k][j][IEND]
                q[nv][k][j][IEND+2] = q[nv][k][j][IEND-1]
          remember: IEND is the last cell index inside the real domain.
    */
    BOX_LOOP(box,k,j,i) q[nv][k][j][i] = s*q[nv][k][j][2*i_cap_inter_end-i+1];

  } else if (side == CAP_WALL_EXTERNAL){  
    BOX_LOOP(box,k,j,i) q[nv][k][j][i] = s*q[nv][k][2*(j_cap_inter_end+1)-j-1][i];

  } else if (side == CAP_WALL_CORNER_INTERNAL) {
    pp = 0;
    BOX_LOOP(box,k,j,i) {
      d_correction[IDIR].Vc[nv][pp] = s*q[nv][k][j][2*i_cap_inter_end-i+1];
      d_correction[IDIR].i[pp] = i;
      d_correction[IDIR].j[pp] = j;
      d_correction[IDIR].k[pp] = k;
      pp++;
    }
    if (pp != d_correction[IDIR].Npoints) {
      print1("[ReflectiveBoundCap] Not all the correction cells(IDIR) have been filled!");
      QUIT_PLUTO(1);
    }

  } else if ( side == CAP_WALL_CORNER_EXTERNAL) {
    pp = 0;
    BOX_LOOP(box,k,j,i) {
      d_correction[JDIR].Vc[nv][pp] = s*q[nv][k][2*(j_cap_inter_end+1)-j-1][i];
      d_correction[JDIR].i[pp] = i;
      d_correction[JDIR].j[pp] = j;
      d_correction[JDIR].k[pp] = k;
      pp++;
    }
    if (pp != d_correction[JDIR].Npoints) {
      print1("[ReflectiveBoundCap] Not all the correction cells(JDIR) have been filled!");
      QUIT_PLUTO(1);
    }

  } else {
    print1("\n[ReflectiveBoundCap] Wrong choice for 'side'");
    QUIT_PLUTO(1);
  }
}

/***********************************************************
 * Do all the necessary stuff for variables which are not evolved
 * in time (at the moment of writing: simply set to 0 the corrections d_correction)
 * *********************************************************/
void SetNotEvolvedVar (int nv) {
  ZeroBoundCap (NULL, nv, 0, CAP_WALL_CORNER_INTERNAL, CENTER);
  ZeroBoundCap (NULL, nv, 0, CAP_WALL_CORNER_EXTERNAL, CENTER);
}

/*************************************************************
 * Set to 0 the correction for a certain variable in a certain direction
 * ***********************************************************/
void ZeroBoundCap (double ****q, int nv, int s, int side, int vpos)
/*!
 * [Created by Ema]
 * Set a boundary(ghost cells defined by some box or the corrections d_correction.Vc[][] ...) to value 0.
 *
 * \param [in,out] q   a 3D flow quantity
 * \param [in]    nv    kind of variable to apply the bc, e.g.: RHO, iBPHI, PRS...
 *   
 *********************************************************************** */
{
  int   i, j, k, pp;
  RBox *box = GetRBoxCap(side, vpos);

  if (not_allocated_d_correction && (side == CAP_WALL_CORNER_INTERNAL || side == CAP_WALL_CORNER_EXTERNAL)) {
    // I allocate memory for d_correction
    // [Rob] I could define a function to do this in two lines
    // First I count how many points are needed
    pp = 0;
    BOX_LOOP(box, k, j, i) pp++;

    d_correction[IDIR].Npoints = pp;
    d_correction[IDIR].i = ARRAY_1D(d_correction[IDIR].Npoints, int);
    d_correction[IDIR].j = ARRAY_1D(d_correction[IDIR].Npoints, int);
    d_correction[IDIR].k = ARRAY_1D(d_correction[IDIR].Npoints, int);
    d_correction[IDIR].Vc = ARRAY_2D( NVAR, d_correction[IDIR].Npoints, double);

    d_correction[JDIR].Npoints = pp;
    d_correction[JDIR].i = ARRAY_1D(d_correction[JDIR].Npoints, int);
    d_correction[JDIR].j = ARRAY_1D(d_correction[JDIR].Npoints, int);
    d_correction[JDIR].k = ARRAY_1D(d_correction[JDIR].Npoints, int);
    d_correction[JDIR].Vc = ARRAY_2D( NVAR, d_correction[JDIR].Npoints, double);

    not_allocated_d_correction = 0;
  }

  if (side == CAP_WALL_INTERNAL) {
    BOX_LOOP(box,k,j,i) q[nv][k][j][i] = 0;

  } else if (side == CAP_WALL_EXTERNAL){  
    BOX_LOOP(box,k,j,i) q[nv][k][j][i] = 0;

  } else if (side == CAP_WALL_CORNER_INTERNAL) {
    pp = 0;
    BOX_LOOP(box,k,j,i) {
      d_correction[IDIR].Vc[nv][pp] = 0;
      d_correction[IDIR].i[pp] = i;
      d_correction[IDIR].j[pp] = j;
      d_correction[IDIR].k[pp] = k;
      pp++;
    }
    if (pp != d_correction[IDIR].Npoints) {
      print1("[ReflectiveBoundCap] Not all the correction cells(IDIR) have been filled!");
      QUIT_PLUTO(1);
    }

  } else if ( side == CAP_WALL_CORNER_EXTERNAL) {
    pp = 0;
    BOX_LOOP(box,k,j,i) {
      d_correction[JDIR].Vc[nv][pp] = 0;
      d_correction[JDIR].i[pp] = i;
      d_correction[JDIR].j[pp] = j;
      d_correction[JDIR].k[pp] = k;
      pp++;
    }
    if (pp != d_correction[JDIR].Npoints) {
      print1("[ReflectiveBoundCap] Not all the correction cells(JDIR) have been filled!");
      QUIT_PLUTO(1);
    }

  } else {
    print1("\n[ReflectiveBoundCap] Wrong choice for 'side'");
    QUIT_PLUTO(1);
  }
}

/******************************************************
 * IsOutCone : Checks whether a certain point is
 * outside from a cone (departing from the corner of the capillary)
 * with a defined angle
 * ***************************************************/
int IsOutCone(double angle, double r, double z) {

  if (atan((r-rcap_real) / (z-zcap_real)) >= angle && r >= rcap_real) {
    return 1;
  } else {
    return 0;
  }
}

#if MULTIPLE_GHOSTS == YES
  /***********************************************
  * Author :  Ema
  * date : 03/01/18
  * Purpose: Apply multiple ghost cells in internal boundary,
  *          which means overwrite the present Data *d in certain points with
  *          values which depends on the integration direction
  *
  ***********************************************/
  void ApplyMultipleGhosts(const Data *d, int direction) {
    int nv, pp, k,j,i;
    for (pp = 0; pp < d_correction[direction].Npoints; pp++){
      i = d_correction[direction].i[pp];
      j = d_correction[direction].j[pp];
      k = d_correction[direction].k[pp];
      VAR_LOOP(nv) d->Vc[nv][k][j][i] = d_correction[direction].Vc[nv][pp];
    }
  }
#endif