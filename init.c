/* ///////////////////////////////////////////////////////////////////// */
/* fatto da Ema partendo dall'esempio di Field diffusion, 2/8/2017
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"
#include "gamma_transp.h"
#include "capillary_wall.h"
#include "current_table.h"

/*Auxiliary function to set the temperature*/
void setT(const Data *d, double T, int i, int j, int k);

/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *
 *********************************************************************** */
{
  double T,mu;/*Temperature in K and mean particle weight*/
  double curr, Bwall; //Bwall is in code units
  double unit_Mfield, dens0;

  // Just a check that the geometrical settings makes sense:
  if (DZCAP > ZCAP){
    print1("\nElectrode is longer than whole capillary!");
    QUIT_PLUTO(1);
  }

  unit_Mfield = COMPUTE_UNIT_MFIELD(UNIT_VELOCITY, UNIT_DENSITY);
  dens0 = (DENS0)/UNIT_DENSITY;

  curr = current_from_time(0.0);
  // print1("Current from tab: %g", curr);
  // Mag field at the capillary wall, in code units
  Bwall = (BIOTSAV_GAUSS_S_A(curr, RCAP))/unit_Mfield;

  #if GEOMETRY != CYLINDRICAL
   #error geometry not valid
  #endif
  #if PHYSICS != MHD
   #error physics not valid
  #endif

  //Remember: in cyl coords x1 is r, x2 is z

  /* ----------------------------------------------------- 
      Zones not covered in the next lines (except for zone "Everywhere")
    ----------------------------------------------------- */
  us[RHO] = 0.1*dens0;
  /* ----------------------------------------------------- 
      Inside capillary, excluded near-electrode zone
     ----------------------------------------------------- */
  if (x2 < zcap-dzcap && x1 <= rcap) {
    us[iBPHI] = Bwall*x1/rcap;
    us[RHO] = dens0;
  }
  /* ----------------------------------------------------- 
      Inside capillary, in near-electrode zone
     ----------------------------------------------------- */
  if (zcap-dzcap <= x2 && x2 < zcap && x1 < rcap) {
    /* the B field linearly decreses in z direction
    (this is provisory, better electrode have to be implemented) */
    us[iBPHI] = (Bwall*x1/rcap) * ( 1 - (x2 - (zcap-dzcap))/dzcap );
    us[RHO] = dens0;
  }
  /* ------------------------------------------------------
      Above non-electrode wall (internal boundary, outside capillary)
     ------------------------------------------------------ */
  if (x2 < zcap-dzcap && x1>rcap) {
    us[iBPHI] = Bwall;
  }
  /* ------------------------------------------------------
      Above electrode wall (internal boundary, outside capillary)
     ------------------------------------------------------ */
  if ( zcap-dzcap <= x2 && x2 < zcap && x1 >= rcap) {
    us[iBPHI] = Bwall * ( 1 - (x2 - (zcap-dzcap)) / dzcap );
  }
  /* ------------------------------------------------------
      Outside capillary (not in internal boundary)
     ------------------------------------------------------ */
  if (x2 > zcap) {
    // No field outside capillary
    us[iBPHI] = 0.0;
  }
  /* ----------------------------------------------------- 
      Everywhere
     ----------------------------------------------------- */
  us[iBZ] = us[iBR] = 0.0;
  us[iVPHI] = us[iVZ] = us[iVR] = 0.0;
  T = T0;
  #if EOS==IDEAL
      mu = MeanMolecularWeight(us);
  #elif EOS==PVTE_LAW
      GetMu(T, us[RHO], &mu); // GetMu takes T in Kelvin, no need to adim. T
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
  double t_sec, curr, Bwall; //Bwall is in code units, t_sec is in seconds
  int  vsign[NVAR]; /*vector containing signs which will be set by Flipsign*/
  // double T,mu;/*Temperature in K and mean particle weight, for the usage of macro "KELVIN" see page 45 of the manual*/
  // double mu_all[NX3_TOT][NX2_TOT][NX1_TOT]; /*mean particle weight in the whole domain*/
  double qz,qr,diagonal,sinth,costh;
  double mu;
  double unit_Mfield;
  static int first_call=1;
  // static int i_cap_inter_end=-1, j_cap_inter_end=-1,j_elec_start=-1;

  unit_Mfield = COMPUTE_UNIT_MFIELD(UNIT_VELOCITY, UNIT_DENSITY);

  /*[Ema] g_time è: "The current integration time."(dalla docuementazione in Doxigen) */
  t_sec = g_time*(UNIT_LENGTH/UNIT_VELOCITY);

  curr = current_from_time(t_sec);
  // print1("\nCurrent from tab: %g", curr);
  Bwall = BIOTSAV_GAUSS_S_A(curr, RCAP)/unit_Mfield;

  /**********************************
  Find the remarkable indexes (if they had not been found before)
  ***********************************/
  if (capillary_not_set) {
    if (SetRemarkableIdxs(grid)){
      print1("\nError while setting remarkable points!");
      QUIT_PLUTO(1);
    }
  }
  if (first_call) {
    /* Set internal boundary flag on internal boundary points*/
    KTOT_LOOP(k) {
      for (j=0; j<=j_cap_inter_end; j++) {
        for (i=i_cap_inter_end+1; i<NX1_TOT; i++) {
          d->flag[k][j][i] |= FLAG_INTERNAL_BOUNDARY;
        }
      }
    }
    /* Flatten the variables to conveniente values in points
       in internal boundary (except for "ghosts") */
    /* WARNING!! IN CASE OF PRESSURE/TEMPERATURE TABLE INTERPOLATION ERROR, IT MIGHT
       BE ADVISABLE TO CHANGE THE VALUES HERE!!*/
    KTOT_LOOP(k) {
      for (j=0; j<=j_cap_inter_end-1; j++){
        for (i=i_cap_inter_end+2; i<NX1_TOT; i++) {
          d->Vc[RHO][k][j][i] = 1e-2;
          d->Vc[iVZ][k][j][i] = 1e-3;
          d->Vc[iVR][k][j][i] = 1e-3;
          d->Vc[PRS][k][j][i] = 1e-3;
          d->Vc[iBPHI][k][j][i] = 1e-3;
        }
      }
    }
    first_call = 0;
  }

  /**********************************
  ***********************************
  Actual setting of the boundary conditions
  ***********************************
  ***********************************/
  if (side == X1_END){
  /**********************************************
  Side r = rmax
  **********************************************/
    if (box->vpos == CENTER) {
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
    (for RHO, iVR, iVZ, T: I stop at j_cap_inter_end-1 to exclude the corner cell)
    ************************/
    KTOT_LOOP(k) {
      // Fluid "Reflection"
      for (j=0; j<=j_cap_inter_end-1; j++) {
        d->Vc[RHO][k][j][i_cap_inter_end+1] = d->Vc[RHO][k][j][i_cap_inter_end];
        d->Vc[iVR][k][j][i_cap_inter_end+1] = -(d->Vc[iVR][k][j][i_cap_inter_end]);
        d->Vc[iVZ][k][j][i_cap_inter_end+1] = d->Vc[iVZ][k][j][i_cap_inter_end];
      }
      // Temperature
      for (j=0; j<=j_cap_inter_end-1; j++) {
        setT( d, TWALL, i_cap_inter_end+1, j, k);
      }
      // Magnetic field on capillary wall (exclued electrode)
      for (j=0; j<j_elec_start; j++) {
        d->Vc[iBPHI][k][j][i_cap_inter_end+1] = Bwall;
      }
      // Magnetic field on electrode
      for (j=j_elec_start; j<=j_cap_inter_end; j++) {
        /* d->Vc[iBPHI][k][j][i_cap_inter_end+1] = Bwall* \
             (1-(grid[1].x_glob[j]-grid[1].x_glob[j_elec_start])/ \
             (grid[1].x_glob[j_cap_inter_end]-grid[1].x_glob[j_elec_start])); */
        d->Vc[iBPHI][k][j][i_cap_inter_end+1] = Bwall* \
               (1 - (grid[JDIR].x_glob[j]-(zcap_real-dzcap_real))/dzcap );
      }
    }
    /***********************
    Capillary wall z=cost
    (I start from i=i_cap_inter_end+2 to exclude corner cell)
    ************************/
    KTOT_LOOP(k) {
      for (i=i_cap_inter_end+2; i<NX1_TOT; i++) { // I start from i_cap_inter_end+2 to exclude the corner cell
        // Fluid "Reflection"
        d->Vc[RHO][k][j_cap_inter_end][i] = d->Vc[RHO][k][j_cap_inter_end+1][i];
        d->Vc[iVR][k][j_cap_inter_end][i] = d->Vc[iVR][k][j_cap_inter_end+1][i];
        d->Vc[iVZ][k][j_cap_inter_end][i] = -(d->Vc[iVZ][k][j_cap_inter_end+1][i]);
        // Magnetic Field
        d->Vc[iBPHI][k][j_cap_inter_end][i] = 0.0;
        // Temperature
        setT( d, TWALL, i, j_cap_inter_end, k);
      }
    }
    /*********************
    Corner point (I write some possible algorithms to see which is better,
                  decomment the one you like!)
    *********************/
    #if MULTIPLE_GHOSTS != YES
      KTOT_LOOP(k) {
        /****ALGORITHM 1 *****/
        // // I should not change the grid size near the capillary end!
        // diagonal = sqrt( pow(grid[0].dx_glob[i_cap_inter_end+1],2) + pow(grid[1].dx_glob[j_cap_inter_end],2) );
        // // I compute cos(theta) and sin(theta)
        // sinth = grid[0].dx_glob[i_cap_inter_end+1] / diagonal;
        // costh = grid[1].dx_glob[j_cap_inter_end] / diagonal;
        // //qr = 0.5*(rhoA*vrA + rhoB*vrB)
        // qr = d->Vc[RHO][k][j_cap_inter_end][i_cap_inter_end]*d->Vc[iVR][k][j_cap_inter_end][i_cap_inter_end];
        // qr += d->Vc[RHO][k][j_cap_inter_end+1][i_cap_inter_end+1]*d->Vc[iVR][k][j_cap_inter_end+1][i_cap_inter_end+1];
        // qr = 0.5*qr;
        // //qz =  0.5*(rhoA*vzA + rhoB*vzB)
        // qz = d->Vc[RHO][k][j_cap_inter_end][i_cap_inter_end]*d->Vc[iVZ][k][j_cap_inter_end][i_cap_inter_end];
        // qz += d->Vc[RHO][k][j_cap_inter_end+1][i_cap_inter_end+1]*d->Vc[iVZ][k][j_cap_inter_end+1][i_cap_inter_end+1];
        // qz = 0.5*qz;
        // // now I set the actual values
        // d->Vc[RHO][k][j_cap_inter_end][i_cap_inter_end+1] = 0.5*(d->Vc[RHO][k][j_cap_inter_end][i_cap_inter_end]+d->Vc[RHO][k][j_cap_inter_end+1][i_cap_inter_end+1]);
        // d->Vc[iVR][k][j_cap_inter_end][i_cap_inter_end+1] = (qr*sinth+qz*costh)*sinth - 0.5*qr;
        // d->Vc[iVZ][k][j_cap_inter_end][i_cap_inter_end+1] = (qr*sinth+qz*costh)*costh - 0.5*qz;

        /****ALGORITHM 2 *****/
        // /***********************/
        // /*I try to set 0 speed on corner ghost cell (experimental)*/
        // /***********************/
        // d->Vc[RHO][k][j_cap_inter_end][i_cap_inter_end+1] = 0.5*(d->Vc[RHO][k][j_cap_inter_end][i_cap_inter_end]+d->Vc[RHO][k][j_cap_inter_end+1][i_cap_inter_end+1]);
        // d->Vc[iVR][k][j_cap_inter_end][i_cap_inter_end+1] = 0.0;
        // d->Vc[iVZ][k][j_cap_inter_end][i_cap_inter_end+1] = 0.0;
        // /******************/
        // /*Also enforcing zero orthgonal velocity to wall in corner point (experimental!)*/
        // /******************/
        // d->Vc[iVR][k][j_cap_inter_end][i_cap_inter_end] = 0.0;
        // d->Vc[iVZ][k][j_cap_inter_end+1][i_cap_inter_end+1] = 0.0;

        /****ALGORITHM 3 *****/
        /***********************/
        /* Corner ghost cell reflects the cell just below*/
        /***********************/
        // d->Vc[RHO][k][j_cap_inter_end][i_cap_inter_end+1] = d->Vc[RHO][k][j_cap_inter_end][i_cap_inter_end];
        // d->Vc[iVR][k][j_cap_inter_end][i_cap_inter_end+1] = -(d->Vc[iVR][k][j_cap_inter_end][i_cap_inter_end]);
        // d->Vc[iVZ][k][j_cap_inter_end][i_cap_inter_end+1] = d->Vc[iVZ][k][j_cap_inter_end][i_cap_inter_end];
        // /***********************/
        // /* While the cell on the right is enforced to have zero velocity along z */
        // /***********************/
        // d->Vc[iVZ][k][j_cap_inter_end+1][i_cap_inter_end+1] = 0.0;

        /****ALGORITHM 4 *****/
        /***********************/
        /* Ghost cell has 0 speed (and rho is the average of the 2 neighbouring inside
        the domain) and the neighbouring cells reverse their orthogonal speeds*/
        /***********************/
        d->Vc[iVR][k][j_cap_inter_end][i_cap_inter_end] = -(d->Vc[iVR][k][j_cap_inter_end][i_cap_inter_end]);
        d->Vc[iVZ][k][j_cap_inter_end+1][i_cap_inter_end+1] = -(d->Vc[iVZ][k][j_cap_inter_end+1][i_cap_inter_end+1]);
        d->Vc[RHO][k][j_cap_inter_end][i_cap_inter_end+1] = 0.5*(d->Vc[RHO][k][j_cap_inter_end][i_cap_inter_end]\
                                              + d->Vc[RHO][k][j_cap_inter_end+1][i_cap_inter_end+1]);
        d->Vc[iVR][k][j_cap_inter_end][i_cap_inter_end+1] = 0.0;
        d->Vc[iVZ][k][j_cap_inter_end][i_cap_inter_end+1] = 0.0;
      }
    #elif MULTIPLE_GHOSTS == YES

      /****ALGORITHM 5 *****/
      /* Define multiple ghosts on the corner wall cell! */
      /*********************/
      // To the normal corner cell I give the correct reflective boundary for
      // the r direction
      KTOT_LOOP(k) {
        d->Vc[RHO][k][j_cap_inter_end][i_cap_inter_end+1] = d->Vc[RHO][k][j_cap_inter_end][i_cap_inter_end];
        d->Vc[iVR][k][j_cap_inter_end][i_cap_inter_end+1] = -(d->Vc[iVR][k][j_cap_inter_end][i_cap_inter_end]);
        d->Vc[iVZ][k][j_cap_inter_end][i_cap_inter_end+1] = d->Vc[iVZ][k][j_cap_inter_end][i_cap_inter_end];
        //magnetic field
        d->Vc[iBPHI][k][j_cap_inter_end][i_cap_inter_end+1] = 0.0;
        //temperature
        // setT( d, TWALL, i_cap_inter_end+1, j, k);
        setT( d, TWALL, i_cap_inter_end+1, j_cap_inter_end, k); //[Ema?] appena cambiato, vedere se è ok
      }

      // I repeat the usual configuration for the correction in r direction
      d_correction[0].Npoints = 1*1*NX3_TOT;
      d_correction[0].i = ARRAY_1D(d_correction[0].Npoints, int);
      d_correction[0].j = ARRAY_1D(d_correction[0].Npoints, int);
      d_correction[0].k = ARRAY_1D(d_correction[0].Npoints, int);
      d_correction[0].Vc = ARRAY_2D( NVAR, d_correction[0].Npoints, double);
      // I assign values to the correction:
      KTOT_LOOP(k) {
        d_correction[0].i[k] = i_cap_inter_end+1;
        d_correction[0].j[k] = j_cap_inter_end;
        d_correction[0].k[k] = k;
        d_correction[0].Vc[iVR][k] = -(d->Vc[iVR][k][j_cap_inter_end][i_cap_inter_end]);
        d_correction[0].Vc[iVZ][k] = d->Vc[iVZ][k][j_cap_inter_end][i_cap_inter_end];
        d_correction[0].Vc[iVPHI][k] = 0.0;
        d_correction[0].Vc[RHO][k] = d->Vc[RHO][k][j_cap_inter_end][i_cap_inter_end];
        #if EOS==IDEAL
            #error double internal ghost not implemented for ideas eos
        #elif EOS==PVTE_LAW
            GetMu(TWALL, d_correction[0].Vc[RHO][k], &mu); // No need to adim. T (GetMu takes T in K)
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
      d_correction[1].i = ARRAY_1D(d_correction[1].Npoints, int);
      d_correction[1].j = ARRAY_1D(d_correction[1].Npoints, int);
      d_correction[1].k = ARRAY_1D(d_correction[1].Npoints, int);
      d_correction[1].Vc = ARRAY_2D( NVAR, d_correction[1].Npoints, double);
      // I assign values to the correction:
      KTOT_LOOP(k) {
        d_correction[1].i[k] = i_cap_inter_end+1;
        d_correction[1].j[k] = j_cap_inter_end;
        d_correction[1].k[k] = k;
        d_correction[1].Vc[iVR][k] = d->Vc[iVR][k][j_cap_inter_end+1][i_cap_inter_end+1];
        d_correction[1].Vc[iVZ][k] = -(d->Vc[iVZ][k][j_cap_inter_end+1][i_cap_inter_end+1]);
        d_correction[0].Vc[iVPHI][k] = 0.0;
        d_correction[1].Vc[RHO][k] = d->Vc[RHO][k][j_cap_inter_end+1][i_cap_inter_end+1];
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
     Set internal boundary flag on internal boundary points
    **********************/
    /*** At every step I must set the flag, at the program resets it automatically***/
    KTOT_LOOP(k) {
      for (j=0; j<=j_cap_inter_end; j++) {
        for (i=i_cap_inter_end+1; i<NX1_TOT; i++) {
          d->flag[k][j][i] |= FLAG_INTERNAL_BOUNDARY;
        }
      }
    }
    /*** ***/
  }
}

/*-------------------------------------------------------------------------*/
/*Auxiliary function to set the temperature*/
/*-------------------------------------------------------------------------*/
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
  // print1("\nKELVIN:%g", KELVIN);
  // print("\nT:%g",T);
  // print("\nd->Vc[RHO][%d]%d][%i]:%g",d->Vc[RHO][k][j][i],k,j,i);
  // print("\nmu:%g",mu);
  d->Vc[PRS][k][j][i] = d->Vc[RHO][k][j][i]*T / (KELVIN*mu);
}
