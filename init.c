/* ///////////////////////////////////////////////////////////////////// */
/* fatto da Ema partendo dall'esempio di Field diffusion, 2/8/2017
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

//#define BWALL 2.758e-1
#define BWALL 1.1*2.758e-1 // old : 0.1*2.758e-1
#define T0 5000.0
#define TWALL 5000.0
#define RLIMIT1 (3/4*100)
#define RLIMIT2 100 // punto finale in r del dominio, questo lo devi cambiare se cambi la lunghezza in r del dominio
//void BoundValues (double *v, double x1, double x2, double x3, double t);
/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *
 *********************************************************************** */
{
  /*added by Ema*/
  double T,mu;/*Temperature in K and mean particle weight*/


  #if GEOMETRY != CYLINDRICAL
   #error geometry not valid
  #endif
  #if PHYSICS != MHD
   #error physics not valid
  #endif

  us[iBPHI] = BWALL*x1/100.0;
  us[iBZ] = us[iBR] = 0.0;

  us[iVPHI] = us[iVZ] = us[iVR] = 0.0;
  //us[VX2] = 0.0;
  //us[VX3] = 0.0;
  us[RHO] = 10.0;


  // I set a profile for T, with T constant up to a certain R (which is 3/4 R), and then it decays as cos² to TWALL towards the wall
  // if (x1<=RLIMIT1){
  //   T = T0+TWALL; /*T in Kelvin*/
  // } else if (x1>RLIMIT1){
  //   T = T0* pow(cos(M_PI/2 * (-RLIMIT1+x1)/(RLIMIT2-RLIMIT1)) ,2) + TWALL; /*T in Kelvin*/
  // }
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
 /*ANDARE AVANTI A MODIFICARE QUI!!*/
{
  int  i, j, k;
  double *r_coord, *z_coord;
  double t;
  int  vsign[NVAR]; /*vector containing signs which will be set by Flipsign*/
  double T,mu;/*Temperature in K and mean particle weight, for the usage of macro "KELVIN" see page 45 of the manual*/


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
    // I keep it fixed cause I want to see the heating but not the movement of the gas
    // TOT_LOOP (k,j,i) {
    //   d->Vc[iVR][k][j][i] = 0.0;
    //   d->Vc[iVZ][k][j][i] = 0.0;
    //   d->Vc[iVPHI][k][j][i] = 0.0;
    // }
  }

}
