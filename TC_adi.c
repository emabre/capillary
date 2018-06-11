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
*****************************************************************************/
void BuildIJ_TC(const Data *d, Grid *grid, Lines *lines,
                   double **Ip, double **Im, double **Jp,
                   double **Jm, double **CI, double **CJ, double **dEdT) {
  int i,j,k;
  int nv, l;
  double kpar, knor, phi; // Thermal conductivity
  double v[NVAR];
  double ****Vc;
  double *inv_dri, *inv_dzi;
  double *zL, *zR;
  double *rL, *rR;
  double *ArR, *ArL;
  double *dVr, *dVz;
  double *r, *z, *theta;
  double T;

  /* -- set a pointer to the primitive vars array --
    I do this because it is done also in other parts of the code
    maybe it makes the program faster or just easier to write/read...*/
  Vc = d->Vc;

  inv_dzi = grid[JDIR].inv_dxi;
  inv_dri = grid[IDIR].inv_dxi;

  theta = grid[KDIR].x;

  r = grid[IDIR].x;
  z = grid[JDIR].x;
  rL = grid[IDIR].xl;
  rR = grid[IDIR].xr;
  zL = grid[JDIR].xl;
  zR = grid[JDIR].xr;
  ArR = grid[IDIR].A;
  ArL = grid[IDIR].A - 1;
  dVr = grid[IDIR].dV;
  dVz = grid[JDIR].dV;

  /*[Opt] This is probably useless, it is here just for debugging purposes*/
  TOT_LOOP(k, j, i) {
    Ip[j][i] = 0.0;
    Im[j][i] = 0.0;
    Jp[j][i] = 0.0;
    Jm[j][i] = 0.0;
    CI[j][i] = 0.0;
    CJ[j][i] = 0.0;
    dEdT[j][i] = 0.0;
  }
  // lines->lidx
  KDOM_LOOP(k) {
    LINES_LOOP(lines[IDIR], l, j, i) {
      /* :::: Ip :::: */
      for (nv=0; nv<NVAR; nv++)
        v[nv] = 0.5 * (Vc[nv][k][j][i] + Vc[nv][k][j][i+1]);
      TC_kappa( v, rR[i], z[j], theta[k], &kpar, &knor, &phi);
      Ip[j][i] = knor*ArR[i]*inv_dri[i];
      /* [Opt] Here I could use the already computed k to compute
      also Im[1,i+1] (since it needs k at the same interface).
      Doing so I would reduce the calls to TC_kappa by almost a factor 1/2
      (obviously I could do the same for Im/Ip) */

      /* :::: Im :::: */
      for (nv=0; nv<NVAR; nv++)
        v[nv] = 0.5 * (Vc[nv][k][j][i] + Vc[nv][k][j][i-1]);
      TC_kappa( v, rL[i], z[j], theta[k], &kpar, &knor, &phi);
      Im[j][i] = knor*ArL[i]*inv_dri[i-1];

      /* :::: Jp :::: */
      for (nv=0; nv<NVAR; nv++)
        v[nv] = 0.5 * (Vc[nv][k][j][i] + Vc[nv][k][j+1][i]);
      TC_kappa( v, r[i], zR[j], theta[k], &kpar, &knor, &phi);
      Jp[j][i] = knor*inv_dzi[j];

      /* :::: Jm :::: */
      for (nv=0; nv<NVAR; nv++)
        v[nv] = 0.5 * (Vc[nv][k][j][i] + Vc[nv][k][j-1][i]);
      TC_kappa( v, r[i], zL[j], theta[k], &kpar, &knor, &phi);
      Jm[j][i] = knor*inv_dzi[j-1];
      #ifdef TEST_ADI
        HeatCapacity_test(v, r[i], z[j], theta[k], &(dEdT[j][i]) );
      #else
        if (GetPV_Temperature(v, &(T) )!=0) {
          print1("\nTC_kappa:[Ema] Error computing temperature!");
        }
        HeatCapacity(v, T, &(dEdT[j][i]) );
      #endif
      /* :::: CI :::: */
      CI[j][i] = dEdT[j][i]*dVr[i];
      
      /* :::: CJ :::: */
      CJ[j][i] = dEdT[j][i]*dVz[j];
    }
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
void BoundaryTC_ADI(Lines lines[2], const Data *d, Grid *grid, double t) {
  int i,j,l;
  double Twall;
  double Twall_K = g_inputParam[TWALL]; // Wall temperature in Kelvin
  // [Err]
  // double L = 0.02/UNIT_LENGTH;

  // I compute the wall temperature
  Twall = Twall_K / KELVIN;

  // IDIR lines
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
  // JDIR lines
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

#endif