/* ///////////////////////////////////////////////////////////////////// */
/* fatto da Ema partendo dall'esempio di Field diffusion, 2/8/2017
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"
#include "gamma_transp.h"

#define BWALL 2.758e-1
// #define BWALL 1.1*2.758e-1 // old : 0.1*2.758e-1
#define T0 5000.0
#define TWALL 5000.0
#define RCAP 0.03
#define ZCAP 3.0 /*the capillary is long 2*ZCAP and wide 2*RCAP */
/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *
 *********************************************************************** */
{
  double T,mu;/*Temperature in K and mean particle weight*/

  #if GEOMETRY != CYLINDRICAL
   #error geometry not valid
  #endif
  #if PHYSICS != MHD
   #error physics not valid
  #endif

  if (x1<(RCAP/UNIT_LENGTH)) { //in cyl coords x1 is r, x2 is z
    us[iBPHI] = BWALL*x1/(RCAP/UNIT_LENGTH);
  } else {
    us[iBPHI] = BWALL;
  }
    us[iBZ] = us[iBR] = 0.0;

    us[iVPHI] = us[iVZ] = us[iVR] = 0.0;
    us[RHO] = 10.0;

    T = T0;
    #if EOS==IDEAL
        mu = MeanMolecularWeight(us);
    #elif EOS==PVTE_LAW
        GetMu(T, us[RHO], &mu);
    #endif
    us[PRS] = us[RHO]*T / (KELVIN*mu); /*for the usage of macro "KELVIN" see page 45 of the manual*/

}

/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/*
 *
 *********************************************************************** */
{}


/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid)
/*
 *
 *********************************************************************** */
{
  int  i, j, k;
  double t;
  int  vsign[NVAR]; /*vector containing signs which will be set by Flipsign*/
  double T,mu;/*Temperature in K and mean particle weight, for the usage of macro "KELVIN" see page 45 of the manual*/

  if (idx_rcap==0 && idx_zcap==0) {
    print1("inidici prima di algoritmo ricerca bordi interni\n");
    print1("idx_rcap: %d, idx_zcap: %d\n",idx_rcap,idx_zcap);
    /* I find the indexes of the cells closest to the capillary bounds*/
    idx_rcap = find_idx_closest(grid[0].x_glob, grid[0].gend-grid[0].gbeg+1, RCAP/UNIT_LENGTH);
    //print1("grid[0].gbeg:%d, grid[0].gend:%d\n",grid[0].gbeg,grid[0].gend );
    //print1("grid[0].x_glob:%g,grid[0].gend-grid[0].gbeg:%d,rcap%g\n",grid[0].x_glob,grid[0].gend-grid[0].gbeg,RCAP);
    KTOT_LOOP(k) JTOT_LOOP(j) {
      for(i=idx_rcap+1;i<NX1_TOT; i++) d->flag[k][j][i] |= FLAG_INTERNAL_BOUNDARY;
    }
  }
 /*[Ema] g_time è: "The current integration time."(dalla docuementazione in Doxigen) */
  t = g_time; /*at the moment unused*/

  if (side == X1_END){
    if (box->vpos == CENTER) {
    /**********************************************
     side r = rmax
     **********************************************/
    // Setting the Magnetic field
      BOX_LOOP(box,k,j,i){
        d->Vc[iBPHI][k][j][i] = BWALL;
        d->Vc[iBZ][k][j][i] = 0.0;
        d->Vc[iBR][k][j][i] = 0.0;
      }

      // Setting v and rho
      FlipSign (X1_END, REFLECTIVE, vsign); // forse va inizializzato vsign??
      ReflectiveBound (d->Vc[RHO], vsign[RHO], X1_END, CENTER);
      ReflectiveBound (d->Vc[iVZ], vsign[iVZ], X1_END, CENTER);
      ReflectiveBound (d->Vc[iVR], vsign[iVR], X1_END, CENTER);
      ReflectiveBound (d->Vc[iVPHI], vsign[iVPHI], X1_END, CENTER);
      //ReflectiveBound (d->Vc[PRS], vsign[PRS], X1_END, CENTER);

      // Setting T (which means pressure!)
      T = TWALL;
      BOX_LOOP(box,k,j,i){
        #if EOS==IDEAL
            mu = MeanMolecularWeight(d->Vc);
        #elif EOS==PVTE_LAW
            GetMu(T, d->Vc[RHO][k][j][i], &mu);
        #endif
          d->Vc[PRS][k][j][i] = d->Vc[RHO][k][j][i]*T / (KELVIN*mu);
      }
    } else {
        print1("[Ema]UserDefBoundary: Not setting BCs!!!!\n");
    }
  } else if (side == 0) {
    /**********************************
    Internal Boundary
    ***********************************/
    // I make reflective boundary using as ghost the cell idx_rcap+1
    T = TWALL;
    KTOT_LOOP(k) JTOT_LOOP(j) {
      d->Vc[RHO][k][j][idx_rcap+1] = d->Vc[RHO][k][j][idx_rcap];
      d->Vc[iVR][k][j][idx_rcap+1] = -(d->Vc[iVR][k][j][idx_rcap]);
      d->Vc[iVZ][k][j][idx_rcap+1] = d->Vc[iVZ][k][j][idx_rcap];
      d->Vc[PRS][k][j][idx_rcap] = d->Vc[RHO][k][j][idx_rcap]*T / (KELVIN*mu);
      d->flag[k][j][idx_rcap+1] |= FLAG_INTERNAL_BOUNDARY;
    }
  }
}
