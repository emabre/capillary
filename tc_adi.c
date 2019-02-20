#include "pluto.h"
#include "adi.h"
#include "capillary_wall.h"
#include "Thermal_Conduction/tc.h"
#include "pvte_law_heat_capacity.h"

#if THERMAL_CONDUCTION  == ALTERNATING_DIRECTION_IMPLICIT

/****************************************************************************
Function to build the Ip,Im,Jp,Jm, CI, CJ (and also dEdT) for the thermal conduction problem
Note that I must make available for outside dEdT, as I will use it later to
advance the energy in a way that conserves the energy
Note: Harmonic averaging of k is done as suggested in paper P.Sharma,G.W.Hammett,"Preserving Monotonicity in Anisotropic Diffusion"(2007)
*****************************************************************************/
void BuildIJ_TC(const Data *d, Grid *grid, Lines *lines,
                   double **Ip, double **Im, double **Jp,
                   double **Jm, double **CI, double **CJ, double **dEdT) {
  static int first_call=1;
  static double **protoIp, **protoIm, **protoJp, **protoJm, **protoCI, **protoCJ;
  int i,j,k;
  int nv, l;
  double kpar, knor, knor2, phi; // Thermal conductivity
  double v[NVAR];
  double ****Vc;
  double *inv_dri, *inv_dzi;
  double *zL, *zR;
  double *rL, *rR;
  double *ArR, *ArL;
  double *dVr, *dVz;
  double *r, *z, *theta;
  double T;
  int Nlines, lidx, ridx;

  /* -- set a pointer to the primitive vars array --
    I do this because it is done also in other parts of the code
    maybe it makes the program faster or just easier to write/read...*/
  Vc = d->Vc;

  theta = grid[KDIR].x;

  r = grid[IDIR].x;
  z = grid[JDIR].x;
  rL = grid[IDIR].xl;
  rR = grid[IDIR].xr;
  zL = grid[JDIR].xl;
  zR = grid[JDIR].xr;

  // I compose the grid-related part once forever (in static variables) and only update eta
  if (first_call) {
    // Name shorthands
    ArR = grid[IDIR].A;
    ArL = grid[IDIR].A - 1;
    dVr = grid[IDIR].dV;
    dVz = grid[JDIR].dV;
    inv_dzi = grid[JDIR].inv_dxi;
    inv_dri = grid[IDIR].inv_dxi;

    /*Allocatin of memory for the proto variables*/
    protoIp = ARRAY_2D(NX2_TOT, NX1_TOT, double);
    protoIm = ARRAY_2D(NX2_TOT, NX1_TOT, double);
    protoJp = ARRAY_2D(NX2_TOT, NX1_TOT, double);
    protoJm = ARRAY_2D(NX2_TOT, NX1_TOT, double);
    protoCI = ARRAY_2D(NX2_TOT, NX1_TOT, double);
    protoCJ = ARRAY_2D(NX2_TOT, NX1_TOT, double);

    /*[Opt] This is probably useless, it is here just for debugging purposes*/
    TOT_LOOP(k, j, i) {
      protoIp[j][i] = 0.0;
      protoIm[j][i] = 0.0;
      protoJp[j][i] = 0.0;
      protoJm[j][i] = 0.0;
      protoCI[j][i] = 0.0;
      protoCJ[j][i] = 0.0;
      dEdT[j][i] = 0.0;
    }
    KDOM_LOOP(k) {
      LINES_LOOP(lines[IDIR], l, j, i) {
        /* :::: Ip :::: */
        protoIp[j][i] = ArR[i]*inv_dri[i];
        /* [Opt] Here I could use the already computed k to compute
        also Im[1,i+1] (since it needs k at the same interface).
        Doing so I would reduce the calls to TC_kappa by almost a factor 1/2
        (obviously I could do the same for Im/Ip) */
        /* :::: Im :::: */
        protoIm[j][i] = ArL[i]*inv_dri[i-1]; // forse allora Im[j][i] = Ip[j][i-1] ??
        /* :::: Jp :::: */
        protoJp[j][i] = inv_dzi[j];
        /* :::: Jm :::: */
        protoJm[j][i] = inv_dzi[j-1];
        /* :::: CI :::: */
        protoCI[j][i] = dVr[i];
        /* :::: CJ :::: */
        protoCJ[j][i] = dVz[j];
      }
    }
    first_call = 0;
  }


  KDOM_LOOP(k) {

    /* :::: Ip and Im :::: */
    Nlines = lines[IDIR].N;
    for (l = 0; l<Nlines; l++) {
      j = lines->dom_line_idx[l];
      lidx = lines[IDIR].lidx[l];
      ridx = lines[IDIR].ridx[l];

      /* :::: Im on i=lidx :::: */
      for (nv=0; nv<NVAR; nv++)
        v[nv] = Vc[nv][k][j][lidx]; // [Err] no, you should use armonic average of kappa
      TC_kappa( v, r[lidx], z[j], theta[k], &kpar, &knor, &phi);

      for (nv=0; nv<NVAR; nv++)
        v[nv] = Vc[nv][k][j][lidx-1];
      TC_kappa( v, r[lidx-1], z[j], theta[k], &kpar, &knor2, &phi);

      knor = 2/(1/knor + 1/knor2);
      Im[j][lidx] = knor*protoIm[j][lidx];

      /* :::: Ip and Im for internal (non boundary) interfaces :::: */
      for (i=lidx; i<ridx; i++) {
        for (nv=0; nv<NVAR; nv++)
          v[nv] = Vc[nv][k][j][i];  // [Err] no, you should use armonic average of kappa
        TC_kappa( v, r[i], z[j], theta[k], &kpar, &knor, &phi);

        for (nv=0; nv<NVAR; nv++)
          v[nv] = Vc[nv][k][j][i+1];
        TC_kappa( v, r[i+1], z[j], theta[k], &kpar, &knor2, &phi);
        
        knor = 2/(1/knor + 1/knor2);
        Ip[j][i] = knor*protoIp[j][i];
        Im[j][i+1] = knor*protoIm[j][i+1];
      }

      /* :::: Ip on i=ridx :::: */
      for (nv=0; nv<NVAR; nv++)
        v[nv] = Vc[nv][k][j][ridx];
      TC_kappa( v, r[ridx], z[j], theta[k], &kpar, &knor, &phi);  // [Err] no, you should use armonic average of kappa
      
      for (nv=0; nv<NVAR; nv++)
        v[nv] = Vc[nv][k][j][ridx+1];
      TC_kappa( v, r[ridx+1], z[j], theta[k], &kpar, &knor2, &phi);
      
      knor = 2/(1/knor + 1/knor2);
      Ip[j][ridx] = knor*protoIp[j][ridx];
    }

    /* :::: Jp and Jm :::: */
    Nlines = lines[JDIR].N;
    for (l = 0; l<Nlines; l++) {
      i = lines[JDIR].dom_line_idx[l];
      lidx = lines[JDIR].lidx[l];
      ridx = lines[JDIR].ridx[l];

      /* :::: Jm on j=lidx :::: */
      for (nv=0; nv<NVAR; nv++)
        v[nv] = Vc[nv][k][lidx][i];  // [Err] no, you should use armonic average of kappa
      TC_kappa( v, r[i], z[lidx], theta[k], &kpar, &knor, &phi);
      
      for (nv=0; nv<NVAR; nv++)
        v[nv] = Vc[nv][k][lidx-1][i];
      TC_kappa( v, r[i], z[lidx-1], theta[k], &kpar, &knor2, &phi);
      
      knor = 2/(1/knor + 1/knor2);
      Jm[lidx][i] = knor*protoJm[lidx][i];

      /* :::: Jp and Jm for internal (non boundary) interfaces :::: */
      for (j=lidx; j<ridx; j++) {
        for (nv=0; nv<NVAR; nv++)
          v[nv] = Vc[nv][k][j][i];  // [Err] no, you should use armonic average of kappa
        TC_kappa( v, r[i], z[j], theta[k], &kpar, &knor, &phi);
        
        for (nv=0; nv<NVAR; nv++)
          v[nv] = Vc[nv][k][j+1][i];  // [Err] no, you should use armonic average of kappa
        TC_kappa( v, r[i], z[j+1], theta[k], &kpar, &knor2, &phi);
        
        knor = 2/(1/knor + 1/knor2);
        Jp[j][i] = knor*protoJp[j][i];
        Jm[j+1][i] = knor*protoJm[j+1][i]; 
      }

      /* :::: Jp on j=ridx :::: */
      for (nv=0; nv<NVAR; nv++)
        v[nv] = Vc[nv][k][ridx][i];  // [Err] no, you should use armonic average of kappa
      TC_kappa( v, r[i], z[ridx], theta[k], &kpar, &knor, &phi);
      
      for (nv=0; nv<NVAR; nv++)
        v[nv] = Vc[nv][k][ridx+1][i];  // [Err] no, you should use armonic average of kappa
      TC_kappa( v, r[i], z[ridx+1], theta[k], &kpar, &knor2, &phi);
      
      knor = 2/(1/knor + 1/knor2);
      Jp[ridx][i] = knor*protoJp[ridx][i];
    }

    // I separate the computation of CI and CJ just to improve code readability,
    // I could also inglobate them in the previous cycles
    LINES_LOOP(lines[IDIR], l, j, i) {
      for (nv=0; nv<NVAR; nv++)
        v[nv] = Vc[nv][k][j][i];
      #ifdef TEST_ADI
        HeatCapacity_test(v, r[i], z[j], theta[k], &(dEdT[j][i]) );
      #else
        if (GetPV_Temperature(v, &(T) )!=0) {
          #if WARN_ERR_COMP_TEMP
            print1("\nTC_kappa:[Ema]Err.comp.temp");
          #endif
        }
        HeatCapacity(v, T, &(dEdT[j][i]) );
      #endif
      /* :::: CI :::: */
      CI[j][i] = dEdT[j][i]*protoCI[j][i];  
      /* :::: CJ :::: */
      CJ[j][i] = dEdT[j][i]*protoCJ[j][i];
    }

  #ifdef DEBUG_BUILDIJ
    printf("\n[BuildIJ_TC] Im:");
    printmat(Im, NX2_TOT, NX1_TOT);
    printf("\n[BuildIJ_TC] Ip:");
    printmat(Ip, NX2_TOT, NX1_TOT);
    printf("\n[BuildIJ_TC] Jm:");
    printmat(Jm, NX2_TOT, NX1_TOT);
    printf("\n[BuildIJ_TC] Jp:");
    printmat(Jp, NX2_TOT, NX1_TOT);
    printf("\n[BuildIJ_TC] CI:");
    printmat(CI, NX2_TOT, NX1_TOT);
    printf("\n[BuildIJ_TC] CJ:");
    printmat(CJ, NX2_TOT, NX1_TOT);
  #endif
  }
}

/**************************************************************************
 * GetHeatCapacity: Computes the derivative dE/dT (E is the internal energy
 * per unit volume, T is the temperature). This function also normalizes
 * the dEdT.
 * ************************************************************************/
void HeatCapacity_test(double *v, double r, double z, double theta, double *dEdT) {
  double c;
//
  /*[Opt] This should be done with a table or something similar (also I should
    make a more general computation as soon as possible)*/
  /*[Opt] I could also include the normalization in the definition so that I
    do not do the computations twice!*/

  // Specific heat (heat to increase of 1 Kelvin 1 gram of Hydrogen)
    // Low temperature case, no ionization
    c = 3/2*(1/CONST_mp)*CONST_kB;
    // High temperature case, full ionization
    // c = 3/2*(2/CONST_mp)*CONST_kB;

    // Heat capacity (heat to increas of 1 Kelvin 1 cm³ of Hydrogen)
    *dEdT = c*v[RHO]*UNIT_DENSITY;


  // Normalization
  *dEdT = (*dEdT)/(UNIT_DENSITY*CONST_kB/CONST_mp);
}

/****************************************************************************
* Function to build the bcs of lines
* In the current implementation of this function Data *d is not used
* but I leave it there since before or later it might be needed
*****************************************************************************/
void BoundaryADI_TC(Lines lines[2], const Data *d, Grid *grid, double t, int dir) {
  int i,j,l;
  double Twall;
  double Twall_K = g_inputParam[TWALL]; // Wall temperature in Kelvin
  // [Err]
  // double L = 0.02/UNIT_LENGTH;

  // I compute the wall temperature
  Twall = Twall_K / KELVIN;

  if (dir == IDIR) {
    /*-----------------------------------------------*/
    /*----  Set bcs for lines in direction IDIR  ----*/
    /*-----------------------------------------------*/
    for (l=0; l<lines[IDIR].N; l++) {
      j = lines[IDIR].dom_line_idx[l];
      /* :::: Axis ::::*/
      lines[IDIR].lbound[TDIFF][l].kind = NEUMANN_HOM;
      lines[IDIR].lbound[TDIFF][l].values[0] = 0.0;
      // [Opt] (I can avoid making so many if..)
      if (j <= j_cap_inter_end) {
        /* :::: Capillary wall :::: */

        // [Err] Decomment next two lines
        lines[IDIR].rbound[TDIFF][l].kind = DIRICHLET;
        lines[IDIR].rbound[TDIFF][l].values[0] = Twall;

        // [Err] Remove next two lines
        // lines[IDIR].rbound[TDIFF][l].kind = NEUMANN_HOM;
        // lines[IDIR].rbound[TDIFF][l].values[0] = 0.0;
      } else {
        /* :::: Outer domain boundary :::: */
        // lines[IDIR].rbound[TDIFF][l].kind = DIRICHLET;
        // // IS THIS OK?? Before or later change here, and also in its equivalent in init.c
        // lines[IDIR].rbound[TDIFF][l].values[0] = Twall;
        lines[IDIR].rbound[TDIFF][l].kind = NEUMANN_HOM;
        lines[IDIR].rbound[TDIFF][l].values[0] = 0.0;
      }
    }
  } else if (dir == JDIR) {
    /*-----------------------------------------------*/
    /*----  Set bcs for lines in direction JDIR  ----*/
    /*-----------------------------------------------*/
    for (l=0; l<lines[JDIR].N; l++) {
      i = lines[JDIR].dom_line_idx[l];
      if (i <= i_cap_inter_end){
        /* :::: Capillary internal (symmetry plane) ::::*/

        // [Err] Decomment next two lines
        lines[JDIR].lbound[TDIFF][l].kind = NEUMANN_HOM;
        lines[JDIR].lbound[TDIFF][l].values[0] = 0.0;

        // [Err] Remove next two lines
        // lines[JDIR].lbound[TDIFF][l].kind = DIRICHLET;
        // lines[JDIR].lbound[TDIFF][l].values[0] = Twall;
      } else {
        /* :::: Outer capillary wall ::::*/
        lines[JDIR].lbound[TDIFF][l].kind = DIRICHLET;
        lines[JDIR].lbound[TDIFF][l].values[0] = Twall;
      }
      /* :::: Outer domain boundary ::::*/
      lines[JDIR].rbound[TDIFF][l].kind = NEUMANN_HOM;
      lines[JDIR].rbound[TDIFF][l].values[0] = 0.0;
    }
  }
}

#endif