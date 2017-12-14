/* ///////////////////////////////////////////////////////////////////// */
/* fatto da Ema partendo dall'esempio di Field diffusion, 2/8/2017
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"
#include "gamma_transp.h"

//#define BWALL 2.758e-1
//#define BWALL 1.1*2.758e-1 // old : 0.1*2.758e-1
#define BWALL 0.0
#define T0 5000.0
#define TWALL 5000.0
//void BoundValues (double *v, double x1, double x2, double x3, double t);

/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *
 *********************************************************************** */
{
  double T,Tout,mu,muout;/*Temperature in K and mean particle weight*/

  #if GEOMETRY != CYLINDRICAL
   #error geometry not valid
  #endif
  #if PHYSICS != MHD
   #error physics not valid
  #endif

  if (x1<rcap) { //in cyl coords x1 is r, x2 is z
    us[iBPHI] = BWALL*x1/100.0;
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
  } else {
    us[iBPHI] = 0.0;
    us[iBZ] = us[iBR] = 0.0;
    us[iVPHI] = us[iVZ] = us[iVR] = 0.0;
    us[RHO] = 1.0;

    Tout = T/2.0;
    #if EOS==IDEAL
        mu = MeanMolecularWeight(us);
    #elif EOS==PVTE_LAW
        GetMu(Tout, us[RHO], &muout);
    #endif
    us[PRS] = us[RHO]*Tout / (KELVIN*muout);
  }
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
  double *r_coord, *z_coord;
  double t;
  int  vsign[NVAR]; /*vector containing signs which will be set by Flipsign*/
  double T,mu;/*Temperature in K and mean particle weight, for the usage of macro "KELVIN" see page 45 of the manual*/

  if (idx_rcap==0 && idx_zcap==0) {
    print1("inidici prima di algoritmo ricerca bordi interni\n");
    print1("idx_rcap: %d, idx_zcap: %d\n",idx_rcap,idx_zcap);
    /* I find the indexes of the cells closest to the capillary bounds*/
    idx_rcap = find_idx_closest(grid[0].x_glob,grid[0].gend-grid[0].gbeg,rcap);
    //print1("grid[0].gbeg:%d, grid[0].gend:%d\n",grid[0].gbeg,grid[0].gend );
    //print1("grid[0].x_glob:%g,grid[0].gend-grid[0].gbeg:%d,rcap%g\n",grid[0].x_glob,grid[0].gend-grid[0].gbeg,rcap);
  }
 /*[Ema] g_time è: "The current integration time."(dalla docuementazione in Doxigen) */
  t = g_time; /*at the moment unused*/

  if (side == X1_END){  /* side r=Rcap*/
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){
        d->Vc[iBPHI][k][j][i] = BWALL;
        d->Vc[iBZ][k][j][i] = 0.0;
        d->Vc[iBR][k][j][i] = 0.0;
      }
      FlipSign (X1_END, REFLECTIVE, vsign); // forse va inizializzato vsign??
      ReflectiveBound (d->Vc[RHO], vsign[RHO], X1_END, CENTER);
      ReflectiveBound (d->Vc[iVZ], vsign[iVZ], X1_END, CENTER);
      ReflectiveBound (d->Vc[iVR], vsign[iVR], X1_END, CENTER);
      ReflectiveBound (d->Vc[iVPHI], vsign[iVPHI], X1_END, CENTER);
      //ReflectiveBound (d->Vc[PRS], vsign[PRS], X1_END, CENTER);
      BOX_LOOP(box,k,j,i){
        #if EOS==IDEAL
            mu = MeanMolecularWeight(d->Vc);
        #elif EOS==PVTE_LAW
            GetMu(T, d->Vc[RHO][k][j][i], &mu);
        #endif
          d->Vc[PRS][k][j][i] = d->Vc[RHO][k][j][i]*(TWALL) / (KELVIN*mu);
      }
    } else {
        print1("[Ema]UserDefBoundary: Not setting BCs!!!!\n");
    }
  } else if (side == 0) { // this is to set internal boundary
    // I make reflective boundary using as ghost the cell idx_rcap+1
    TOT_LOOP(k,j,i){
      if (i==idx_rcap+1) {
        d->Vc[RHO][k][j][i] = d->Vc[RHO][k][j][i-1];
        d->Vc[iVR][k][j][i] = -(d->Vc[iVR][k][j][i-1]);
        d->Vc[iVZ][k][j][i] = d->Vc[iVZ][k][j][i-1];
        d->Vc[iVPHI][k][j][i] = d->Vc[RHO][k][j][i]*(TWALL) / (KELVIN*mu);
        d->flag[k][j][i] |= FLAG_INTERNAL_BOUNDARY;
      }
    }
  }

}
