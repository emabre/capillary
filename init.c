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
#define DENS0 2.5e-6
// #define DENS0 2.5e-7
#define RCAP 0.03
#define DZCAP 0.01 /*the electrodes are wide DZCAP cm*/
#define ZCAP 0.05 /*the capillary is long 2*ZCAP cm and wide 2*RCAP cm*/

/*Auxiliary function to set the temperature*/
void setT(const Data *d, double T, int i, int j, int k);

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

  us[RHO] = (0.1*DENS0/UNIT_DENSITY);
  if (x2<(ZCAP-DZCAP)/UNIT_LENGTH) {
    if (x1<(RCAP/UNIT_LENGTH)) { //in cyl coords x1 is r, x2 is z
      us[iBPHI] = BWALL*x1/(RCAP/UNIT_LENGTH);
      us[RHO] = DENS0/UNIT_DENSITY;
    } else {
      us[iBPHI] = BWALL;
    }
  } else if ( ZCAP-DZCAP<=x2 && x2<=ZCAP ) {
    // the field linearly decreses in z direction (this is provisory, better electrode have to be implemented)
    if (x1<(RCAP/UNIT_LENGTH)) { //in cyl coords x1 is r, x2 is z
      us[iBPHI] = (BWALL*x1/(RCAP/UNIT_LENGTH)) * ( 1 - (x2 - (ZCAP-DZCAP)/UNIT_LENGTH)/(DZCAP/UNIT_LENGTH) );
    } else {
      us[iBPHI] = BWALL * ( 1 - (x2 - (ZCAP-DZCAP)/UNIT_LENGTH)/(DZCAP/UNIT_LENGTH) );
    }
  } else if (x2>ZCAP) {
    // No field outside capillary
    us[iBPHI] = 0.0;
  }
  // us[iBPHI] = 0.0;

  us[iBZ] = us[iBR] = 0.0;
  us[iVPHI] = us[iVZ] = us[iVR] = 0.0;

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
  // double T,mu;/*Temperature in K and mean particle weight, for the usage of macro "KELVIN" see page 45 of the manual*/
  // double mu_all[NX3_TOT][NX2_TOT][NX1_TOT]; /*mean particle weight in the whole domain*/
  double qz,qr,diagonal,sinth,costh;
  double mu;

  /*[Ema] g_time è: "The current integration time."(dalla docuementazione in Doxigen) */
  t = g_time; /*at the moment unused*/

  if (idx_rcap==0 && idx_zcap==0) {
    // print1("inidici prima di algoritmo ricerca bordi interni\n");
    /* I find the indexes of the cells closest to the capillary bounds*/
    idx_rcap = find_idx_closest(grid[0].x_glob, grid[0].gend-grid[0].gbeg+1, RCAP/UNIT_LENGTH);
    idx_zcap = find_idx_closest(grid[1].x_glob, grid[1].gend-grid[1].gbeg+1, ZCAP/UNIT_LENGTH);
    idx_start_electr = find_idx_closest(grid[1].x_glob, grid[1].gend-grid[1].gbeg+1, (ZCAP-DZCAP)/UNIT_LENGTH);

    print1("idx_rcap: %d, idx_zcap: %d, idx_start_electr\n",idx_rcap,idx_zcap, idx_start_electr);

    //print1("grid[0].gbeg:%d, grid[0].gend:%d\n",grid[0].gbeg,grid[0].gend );
    //print1("grid[0].x_glob:%g,grid[0].gend-grid[0].gbeg:%d,rcap%g\n",grid[0].x_glob,grid[0].gend-grid[0].gbeg,RCAP);

    /* Capillary:
                              j=idx_start_electr
                               |     j=idx_zcap (ghost)
        r                      |      |
        ^    |                 |      *
        |    |       wall      |      *
             |                 v      *
i=idx_rcap+1 |****(ghosts)*****o*******|
i=idx_rcap   |                        j=idx_zcap+1 (first outside, not ghost)
             |
i=0          |________________________(axis)
            j=0    -> z
    // I should not change the grid size exacly on the capillary end!
    */
    KTOT_LOOP(k) {
      for (j=0; j<=idx_zcap; j++) {
        for (i=idx_rcap+1; i<NX1_TOT; i++) {
          d->flag[k][j][i] |= FLAG_INTERNAL_BOUNDARY;
        }
      }
    }
  }

  if (side == X1_END){
    if (box->vpos == CENTER) {
    /**********************************************
     side r = rmax
     **********************************************/
     // Setting the Magnetic field
      BOX_LOOP(box,k,j,i){
        d->Vc[iBPHI][k][j][i] = 0.0;
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

      // Setting T
      BOX_LOOP(box,k,j,i){
        setT( d, TWALL, i, j, k);
      }
    } else {
        print1("[Ema]UserDefBoundary: Not setting BCs!!!!\n");
    }
  } else if (side == 0) {
    /**********************************
    ***********************************
    Internal Boundary
    ***********************************
    ***********************************/

    /***********************
    Capillary wall r=cost
    ************************/
    KTOT_LOOP(k) {
      // rho and v
      for (j=0; j<=idx_zcap-1; j++) { // I stop at idx_zcap-1 to exclude the corner cell
        d->Vc[RHO][k][j][idx_rcap+1] = d->Vc[RHO][k][j][idx_rcap];
        d->Vc[iVR][k][j][idx_rcap+1] = -(d->Vc[iVR][k][j][idx_rcap]);
        d->Vc[iVZ][k][j][idx_rcap+1] = d->Vc[iVZ][k][j][idx_rcap];
      }
      // Temperature
      for (j=0; j<=idx_zcap-1; j++) {
        setT( d, TWALL, idx_rcap+1, j, k);
      }
      // magnetic field on capillary wall
      for (j=0; j<idx_start_electr; j++) {
        d->Vc[iBPHI][k][j][idx_rcap+1] = BWALL;
      }
      // magnetic field on electrode (provisory)
      for (j=idx_start_electr; j<=idx_zcap; j++) {
        // Sistemare meglio, usare le posizioni dei punti dove davvero inizia
        //l'elettrodo e le altre cose, anzichè le macro
        // d->Vc[iBPHI][k][j][idx_rcap+1] = BWALL*(1-(grid[1].x_glob[j]-(ZCAP-DZCAP)/UNIT_LENGTH)/(DZCAP/UNIT_LENGTH) );
        d->Vc[iBPHI][k][j][idx_rcap+1] = BWALL*\
            (1-(grid[1].x_glob[j]-grid[1].x_glob[idx_start_electr])/ \
            (grid[1].x_glob[idx_zcap]-grid[1].x_glob[idx_start_electr]));
        // d->Vc[iBPHI][k][j][idx_rcap+1] = BWALL;
      }
    }
    /***********************
    Capillary wall z=cost
    ************************/
    KTOT_LOOP(k) {
      for (i=idx_rcap+2; i<NX1_TOT; i++) { // I start from idx_rcap+2 to exclude the corner cell
        // v and rho
        d->Vc[RHO][k][idx_zcap][i] = d->Vc[RHO][k][idx_zcap+1][i];
        d->Vc[iVR][k][idx_zcap][i] = d->Vc[iVR][k][idx_zcap+1][i];
        d->Vc[iVZ][k][idx_zcap][i] = -(d->Vc[iVZ][k][idx_zcap+1][i]);
        // magnetic field
        d->Vc[iBPHI][k][idx_zcap][i] = 0.0;
        //temperature
        setT( d, TWALL, idx_rcap+1, j, k);
      }
    }
    #if MULTIPLE_GHOSTS != YES
      KTOT_LOOP(k) {
        /*********************
        Corner point (I write 3 possible algorithms to see which is better,
                      decomment the one you like!)
        *********************/
        /****ALGORITHM 1 *****/
        // // I should not change the grid size near the capillary end!
        // diagonal = sqrt( pow(grid[0].dx_glob[idx_rcap+1],2) + pow(grid[1].dx_glob[idx_zcap],2) );
        // // I compute cos(theta) and sin(theta)
        // sinth = grid[0].dx_glob[idx_rcap+1] / diagonal;
        // costh = grid[1].dx_glob[idx_zcap] / diagonal;
        // //qr = 0.5*(rhoA*vrA + rhoB*vrB)
        // qr = d->Vc[RHO][k][idx_zcap][idx_rcap]*d->Vc[iVR][k][idx_zcap][idx_rcap];
        // qr += d->Vc[RHO][k][idx_zcap+1][idx_rcap+1]*d->Vc[iVR][k][idx_zcap+1][idx_rcap+1];
        // qr = 0.5*qr;
        // //qz =  0.5*(rhoA*vzA + rhoB*vzB)
        // qz = d->Vc[RHO][k][idx_zcap][idx_rcap]*d->Vc[iVZ][k][idx_zcap][idx_rcap];
        // qz += d->Vc[RHO][k][idx_zcap+1][idx_rcap+1]*d->Vc[iVZ][k][idx_zcap+1][idx_rcap+1];
        // qz = 0.5*qz;
        // // now I set the actual values
        // d->Vc[RHO][k][idx_zcap][idx_rcap+1] = 0.5*(d->Vc[RHO][k][idx_zcap][idx_rcap]+d->Vc[RHO][k][idx_zcap+1][idx_rcap+1]);
        // d->Vc[iVR][k][idx_zcap][idx_rcap+1] = (qr*sinth+qz*costh)*sinth - 0.5*qr;
        // d->Vc[iVZ][k][idx_zcap][idx_rcap+1] = (qr*sinth+qz*costh)*costh - 0.5*qz;

        /****ALGORITHM 2 *****/
        // /***********************/
        // /*I try to set 0 speed on corner ghost cell (experimental)*/
        // /***********************/
        // d->Vc[RHO][k][idx_zcap][idx_rcap+1] = 0.5*(d->Vc[RHO][k][idx_zcap][idx_rcap]+d->Vc[RHO][k][idx_zcap+1][idx_rcap+1]);
        // d->Vc[iVR][k][idx_zcap][idx_rcap+1] = 0.0;
        // d->Vc[iVZ][k][idx_zcap][idx_rcap+1] = 0.0;
        // /******************/
        // /*Also enforcing zero orthgonal velocity to wall in corner point (experimental!)*/
        // /******************/
        // d->Vc[iVR][k][idx_zcap][idx_rcap] = 0.0;
        // d->Vc[iVZ][k][idx_zcap+1][idx_rcap+1] = 0.0;

        /****ALGORITHM 3 *****/
        /***********************/
        /* Corner ghost cell reflects the cell just below*/
        /***********************/
        // d->Vc[RHO][k][idx_zcap][idx_rcap+1] = d->Vc[RHO][k][idx_zcap][idx_rcap];
        // d->Vc[iVR][k][idx_zcap][idx_rcap+1] = -(d->Vc[iVR][k][idx_zcap][idx_rcap]);
        // d->Vc[iVZ][k][idx_zcap][idx_rcap+1] = d->Vc[iVZ][k][idx_zcap][idx_rcap];
        // /***********************/
        // /* While the cell on the right is enforced to have zero velocity along z */
        // /***********************/
        // d->Vc[iVZ][k][idx_zcap+1][idx_rcap+1] = 0.0;

        /****ALGORITHM 4 *****/
        /***********************/
        /* Ghost cell has 0 speed (and rho is the average of the 2 neighbouring inside
        the domain) and the neighbouring cells reverse their orthogonal speeds*/
        /***********************/
        d->Vc[iVR][k][idx_zcap][idx_rcap] = -(d->Vc[iVR][k][idx_zcap][idx_rcap]);
        d->Vc[iVZ][k][idx_zcap+1][idx_rcap+1] = -(d->Vc[iVZ][k][idx_zcap+1][idx_rcap+1]);
        d->Vc[RHO][k][idx_zcap][idx_rcap+1] = 0.5*(d->Vc[RHO][k][idx_zcap][idx_rcap]\
                                              + d->Vc[RHO][k][idx_zcap+1][idx_rcap+1]);
        d->Vc[iVR][k][idx_zcap][idx_rcap+1] = 0.0;
        d->Vc[iVZ][k][idx_zcap][idx_rcap+1] = 0.0;
      }
    #elif MULTIPLE_GHOSTS == YES

      /****ALGORITHM 5 *****/
      /* Define multiple ghosts on the corner wall cell! */
      /*********************/
      // To the normal corner cell I give the correct reflective boundary for
      // the r direction
      d->Vc[RHO][k][idx_zcap][idx_rcap+1] = d->Vc[RHO][k][idx_zcap][idx_rcap];
      d->Vc[iVR][k][idx_zcap][idx_rcap+1] = -(d->Vc[iVR][k][idx_zcap][idx_rcap]);
      d->Vc[iVZ][k][idx_zcap][idx_rcap+1] = d->Vc[iVZ][k][idx_zcap][idx_rcap];
      //magnetic field
      d->Vc[iBPHI][k][idx_zcap][idx_rcap+1] = 0.0;
      //temperature
      setT( d, TWALL, idx_rcap+1, j, k);

      // I repeat the usual configuration for the correction in r direction
      d_correction[0].Npoints = 1*1*NX3_TOT;
      d_correction[0].i = ARRAY_1D(d_correction[0].Npoints, int);
      d_correction[0].j = ARRAY_1D(d_correction[0].Npoints, int);
      d_correction[0].k = ARRAY_1D(d_correction[0].Npoints, int);
      d_correction[0].Vc = ARRAY_2D( NVAR, d_correction[0].Npoints, double);
      // I assign values to the correction:
      KTOT_LOOP(k) {
        d_correction[0].i[k] = idx_rcap+1;
        d_correction[0].j[k] = idx_zcap;
        d_correction[0].k[k] = k;
        d_correction[0].Vc[iVR][k] = -(d->Vc[iVR][k][idx_zcap][idx_rcap]);
        d_correction[0].Vc[iVZ][k] = d->Vc[iVZ][k][idx_zcap][idx_rcap];
        d_correction[0].Vc[RHO][k] = d->Vc[RHO][k][idx_zcap][idx_rcap];
        #if EOS==IDEAL
            #error double internal ghost not implemented for ideas eos
        #elif EOS==PVTE_LAW
            GetMu(TWALL, d_correction[0].Vc[RHO][k], &mu);
        #endif
        // I cannot correct the Temperature, I must set it same as in the normal
        // *d structure
        d_correction[0].Vc[PRS][k] = d_correction[0].Vc[RHO][k]*TWALL / (KELVIN*mu);
        d_correction[0].Vc[iBZ][k] = 0.0;
        d_correction[0].Vc[iBPHI][k] = 0.0;
        d_correction[0].Vc[iBR][k] = 0.0;
      }

      // In i and j directions I have only one point to fix,
      // but I probably have multiple points in k (NX3_TOT),
      // as I must consider the bcs.
      // If I don't correct all the points in k direction I might have
      // gradients in k direction and flow of matter, momentum, et cetera
      d_correction[1].Npoints = 1*1*NX3_TOT;
      // I initialize the vectors inside d_correction[1]: .i, .j, .k, .Vc
      d_correction[1].i = ARRAY_1D(d_correction[1].Npoints, int);
      d_correction[1].j = ARRAY_1D(d_correction[1].Npoints, int);
      d_correction[1].k = ARRAY_1D(d_correction[1].Npoints, int);
      d_correction[1].Vc = ARRAY_2D( NVAR, d_correction[1].Npoints, double);
      // I assign values to the correction:
      KTOT_LOOP(k) {
        d_correction[1].i[k] = idx_rcap+1;
        d_correction[1].j[k] = idx_zcap;
        d_correction[1].k[k] = k;
        d_correction[1].Vc[iVR][k] = d->Vc[iVR][k][idx_zcap+1][idx_rcap+1];
        d_correction[1].Vc[iVZ][k] = -(d->Vc[iVZ][k][idx_zcap+1][idx_rcap+1]);
        d_correction[1].Vc[RHO][k] = d->Vc[RHO][k][idx_zcap+1][idx_rcap+1];
        #if EOS==IDEAL
            #error double internal ghost not implemented for ideas eos
        #elif EOS==PVTE_LAW
            GetMu(TWALL, d_correction[1].Vc[RHO][k], &mu);
        #endif
        // I cannot correct the Temperature, I must set it same as in the normal
        // *d structure
        d_correction[1].Vc[PRS][k] = d_correction[1].Vc[RHO][k]*TWALL / (KELVIN*mu);
        d_correction[1].Vc[iBZ][k] = 0.0;
        d_correction[1].Vc[iBPHI][k] = 0.0;
        d_correction[1].Vc[iBR][k] = 0.0;
      }
      // No correction for the k direction
      d_correction[2].Npoints = 0;
    #endif

    /*********************
    Set the flag in the whole wall region
    **********************/
    /*** At every step I must set the flag, at the program resets it automatically***/
    KTOT_LOOP(k) {
      for (j=0; j<=idx_zcap; j++) {
        for (i=idx_rcap+1; i<NX1_TOT; i++) {
          d->flag[k][j][i] |= FLAG_INTERNAL_BOUNDARY;
        }
      }
    }
    /*** ***/
  }
}

/*Auxiliary function to set the temperature*/
void setT(const Data *d, double T, int i, int j, int k) {
  double mu;
  /*I don't do a check on the i,j,k indexes, otherwise
  it would mean doing a lot of if cycles, espectially considering
  that this function is called many times */
  #if EOS==IDEAL
    mu = MeanMolecularWeight(d->Vc);
  #elif EOS==PVTE_LAW
    GetMu(T, d->Vc[RHO][k][j][i], &mu);
  #endif
  d->Vc[PRS][k][j][i] = d->Vc[RHO][k][j][i]*T / (KELVIN*mu);
}
