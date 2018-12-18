/* ///////////////////////////////////////////////////////////////////// */
/* fatto da Ema partendo dall'esempio di Field diffusion, 2/8/2017
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"
#include "gamma_transp.h"
#include "capillary_wall.h"
#include "current_table.h"
#include "adi.h"
#include "debug_utilities.h"

#define AS_DIFF 3

/*For defining how the mass density transition from inside to outside of capillary is done*/
#define SMOOTH_COS2 1
#define UNIFORM_FILL 2
#define FROM_FILE 3

/*Auxiliary function to set the temperature*/
void setT(const Data *d, double T, int i, int j, int k);
/*Function to set the initial mass density according to few choices*/
void SetRhoAnalytic(double *rho, double x1, double x2, double x3, int mode);

// [Err] Decomment the fucntion
/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *
 *********************************************************************** */
{ 
  static int first_call = 1;
  double mu;/*Temperature in K and mean particle weight*/
  double curr, Bwall, B_ghostwall; //Bwall, B_ghostwall in code units
  double unit_Mfield;
  double csi = x1/rcap;
  double alpha = g_inputParam[ALPHA_J]; //ratio between delta current density wall-axis and current density on axis
  double T0_K = g_inputParam[T0];
  double vz0 = g_inputParam[VZ0]/UNIT_VELOCITY;

  // Just a check that the geometrical settings makes sense:
  if (DZCAP > ZCAP){
    print1("\nElectrode is longer than whole capillary!");
    QUIT_PLUTO(1);
  }

  unit_Mfield = COMPUTE_UNIT_MFIELD(UNIT_VELOCITY, UNIT_DENSITY);

  curr = current_from_time(0.0);
  // print1("Current from tab: %g", curr);
  // Mag field at the capillary wall, in code units
  Bwall = (BIOTSAV_GAUSS_A_CM(curr, RCAP))/unit_Mfield;

  B_ghostwall = Bwall; // Here I have no access to the grid structures, so I cannot refine B_ghostwall, as I would do in UserDefBoundary

  #if GEOMETRY != CYLINDRICAL
   #error geometry not valid
  #endif
  #if PHYSICS != MHD
   #error physics not valid
  #endif

  //Remember: in cyl coords x1 is r, x2 is z

  /* ******************************************************** */
  /* ******************************************************** */
  /* Setting rho                                              */
  /* ******************************************************** */
  /* ******************************************************** */
  #if (DENS_INITIAL == SMOOTH_COS2 || DENS_INITIAL == UNIFORM_FILL)
    /******************************************/
    /* rho from analytical function           */
    /******************************************/
    SetRhoAnalytic(&(us[RHO]), x1, x2, x3, DENS_INITIAL);

  #elif (DENS_INITIAL == FROM_FILE)
    /******************************************/
    /* rho interpolated from file             */
    /******************************************/
    // for explanatin see page 47 of userguide
    if (first_call) {
      int input_rho[2];

      if (g_inputParam[DENS0] > 0.0) {
        print1("\nIf rho is from file set DENS0 to a negative value!");
        QUIT_PLUTO(1);
      }

      input_rho[0] = RHO;
      input_rho[1] = -1;
      InputDataSet("grid_ic.out", input_rho);
      InputDataRead("rho_ic.flt", "big");
    }
    InputDataInterpolate(us, x1, x2, x3);

  #else
    #error choice for DENS_INITIAL not understood
  #endif

  /* ******************************************************** */
  /* ******************************************************** */
  /* Setting all the quantities but rho                       */
  /* ******************************************************** */
  /* ******************************************************** */
  /* -----------------------------------------------------
      Zones not covered in the next lines (except for zone "Everywhere")
    ----------------------------------------------------- */
  us[iVZ] = 0.0;
  /* -----------------------------------------------------
      Inside capillary, excluded near-electrode zone
     ----------------------------------------------------- */
  if (x2 < zcap-dzcap && x1 <= rcap) {
    us[iBPHI] = B_ghostwall/(1-0.5*alpha) * csi * (1 - alpha*(1 - 0.5*csi*csi));
    us[iVZ] = vz0;
  }
  /* -----------------------------------------------------
      Inside capillary, in near-electrode zone
     ----------------------------------------------------- */
  if (zcap-dzcap <= x2 && x2 < zcap && x1 < rcap) {
    /* the B field linearly decreses in z direction
    (this is provisory, better electrode have to be implemented) */
    us[iBPHI] = (B_ghostwall/(1-0.5*alpha) * csi * (1 - alpha*(1 - 0.5*csi*csi))) * ( 1 - (x2 - (zcap-dzcap))/dzcap );
    us[iVZ] = vz0;
  }
  /* ------------------------------------------------------
      Above non-electrode wall (internal boundary, outside capillary)
     ------------------------------------------------------ */
  if (x2 < zcap-dzcap && x1>rcap) {
    us[iBPHI] = B_ghostwall;
  }
  /* ------------------------------------------------------
      Above electrode wall (internal boundary, outside capillary)
     ------------------------------------------------------ */
  if ( zcap-dzcap <= x2 && x2 < zcap && x1 >= rcap) {
    us[iBPHI] = B_ghostwall * ( 1 - (x2 - (zcap-dzcap)) / dzcap );
  }
  /* ------------------------------------------------------
      Outside capillary, aligned with capillary (not in internal boundary)
     ------------------------------------------------------ */
  if (x2 > zcap && x1 <= rcap) {
    // No field outside capillary
    us[iBPHI] = 0.0;
  }
  /* ------------------------------------------------------
      Outside capillary, above capillary  (not in internal boundary)
     ------------------------------------------------------ */
  if (x2 > zcap && x1 > rcap) {
    // No field outside capillary
    us[iBPHI] = 0.0;
  }
  /* -----------------------------------------------------
      Everywhere
     ----------------------------------------------------- */
  us[iBZ] = us[iBR] = 0.0;
  us[iVPHI] = us[iVR] = 0.0;
  #if EOS==IDEAL
      mu = MeanMolecularWeight(us);
  #elif EOS==PVTE_LAW
      GetMu(T0_K, us[RHO], &mu); // GetMu takes T in Kelvin, no need to adim. T
  #endif
  us[PRS] = us[RHO]*T0_K / (KELVIN*mu); /*for the usage of macro "KELVIN" see page 45 of the manual*/

  first_call = 0;
}

/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/*
 *
 *********************************************************************** */
{ static double first_call = 1;
  static int ncall_an = -1;  // Number of calls to this function

  ncall_an++;

  #if (EN_CONS_CHECK || PRINT_TIME_INFO)
    double t = g_time*(UNIT_LENGTH/UNIT_VELOCITY);
    double dt = g_dt*(UNIT_LENGTH/UNIT_VELOCITY);
  #endif
  #if EN_CONS_CHECK
    double etot=0, Vtot=0;
    double Mtot=0;
    double current = GetCurrADI();
    double en_adv_in_gau, en_tc_in_gau, en_res_in_gau;
    int i, j, k;
    // int nv;
    // double v[NVAR];
    double ****Vc, ****Uc;
    double dV;
    double unit_en = UNIT_DENSITY*UNIT_VELOCITY*UNIT_VELOCITY*UNIT_LENGTH*UNIT_LENGTH*UNIT_LENGTH;
    double *rR, *rL, *dz;
    RBox *box = GetRBox(DOM, CENTER);

    rR = grid[IDIR].xr_glob;
    rL = grid[IDIR].xl_glob;
    dz = grid[JDIR].dx_glob;

    Vc = d->Vc;
    Uc = d->Uc;

    PrimToCons3D(Vc, Uc, box);

    DOM_LOOP (k,j,i) {
      // I do this to exclude points belonging to the wall
      if (i<=i_cap_inter_end || j>j_cap_inter_end) {
        #if GEOMETRY == CYLINDRICAL
        /* Note that I could use instead some element (like dV) of the grid itself,
          I don't do that to make this chunk of code compatible for both the 2015 and 2018 version of PLUTO */
          dV = CONST_PI*(rR[i]*rR[i]- rL[i]*rL[i])*dz[j];
        #else
          #error Only cyl. geom. is implemented for energy conservation computation
        #endif
        etot += dV*Uc[k][j][i][ENG];
        Vtot += dV;
        Mtot += dV*Uc[k][j][i][RHO];
      }
    }

    // I convert values to physical units
    etot *= unit_en;
    Vtot *= UNIT_LENGTH*UNIT_LENGTH*UNIT_LENGTH;
    Mtot *= UNIT_DENSITY*UNIT_LENGTH*UNIT_LENGTH*UNIT_LENGTH;
    en_adv_in_gau = en_adv_in*unit_en;
    en_tc_in_gau = en_tc_in*unit_en;
    en_res_in_gau = en_res_in*unit_en;    

    /* Write to file (remember: prank is the processor rank (0 in serial mode),
      so this chunk of code should work also in parallel mode!).
    */
    if (prank == 0) {
      char fname[512];
      static double tpos = -1.0;
      FILE *fp;

      sprintf (fname, "%s/energy_cons.dat",RuntimeGet()->output_dir);
      if (g_stepNumber == 0) { /* Open for writing only when we’re starting */
        fp = fopen(fname,"w"); /* from beginning */
        fprintf (fp,"# Energy conservation table. Advice: read with R: read.table()\n");
        fprintf (fp,"%6s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n", "", "t", "dt", "volume", "mass",
                 "current", "Etot", "E_adv_in", "E_tc_in", "E_res_in");
      } else {
        /* Append if this is not step 0 */
        if (tpos < 0.0) { /* Obtain time coordinate of to last written row */
          char sline[512];
          if ((fp = fopen(fname,"r")) != NULL) {
            while (fgets(sline, 512, fp)) {} /* read as many lines as you can, to reach the file end*/
            sscanf(sline, "%*d %lf %*e %*e %*e %*e %*e %*e %*e\n",&tpos); /* read tpos (time of the last written row) from sline */
            fclose(fp);
          } else {
            print1("\n[Analysis] I could not open file %s", fname);
          }
        }
        fp = fopen(fname,"a");
      }
      if (g_time > tpos){
      /* Write if current time if > tpos */
      fprintf (fp, "%6d %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e\n", 
               ncall_an, t, dt, Vtot, Mtot, current, etot,
               en_adv_in_gau, en_tc_in_gau, en_res_in_gau);
      }
      fclose(fp);
    }
  #endif
  
  #if PRINT_TIME_INFO
    static time_t tstart;
    time_t tnow;
    double elapsed_time;

    // I compute the time at the first call (approximately at start or restart of pluto)
    if (first_call) {
      time(&tstart);
    }

    time(&tnow);
    elapsed_time = difftime(tnow, tstart);

    if (prank == 0) {
      char fname[512];
      static double tpos = -1.0;
      FILE *fp;

      sprintf (fname, "%s/times.dat",RuntimeGet()->output_dir);
      if (g_stepNumber == 0) { /* Open for writing only when we’re starting */
        fp = fopen(fname,"w"); /* from beginning */
        fprintf (fp,"# Timing table. Advice: read with R: read.table()\n");
        fprintf (fp,"%6s %12s %12s %12s\n", "", "t", "dt", "elapsed_time");
      } else {
        /* Append if this is not step 0 */
        if (tpos < 0.0) { /* Obtain time coordinate of to last written row */
          char sline[512];
          if ((fp = fopen(fname,"r")) != NULL) {
            while (fgets(sline, 512, fp)) {} /* read as many lines as you can, to reach the file end*/
            sscanf(sline, "%*d %lf %*e %*e\n",&tpos); /* read tpos (time of the last written row) from sline */
            fclose(fp);
          } else {
            print1("\n[Analysis] I could not open file %s", fname);
          }
        }
        fp = fopen(fname,"a");
      }
      if (g_time > tpos){
      /* Write if current time if > tpos */
      fprintf (fp, "%6d %12.6e %12.6e %12.6e\n", 
               ncall_an, t, dt, elapsed_time);
      }
      fclose(fp);
    }

  #endif

  first_call = 0;
}

/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid)
/*
 *
 *********************************************************************** */
{
  int  i, j, k;
  #ifdef FLATTEN_B_OUTCAP
    int j_start_flatten;
  #endif
  double t_sec; // t_sec is in seconds
  int  vsign[NVAR]; /*vector containing signs which will be set by Flipsign*/
  // double T,mu;/*Temperature in K and mean particle weight, for the usage of macro "KELVIN" see page 45 of the manual*/
  // double mu_all[NX3_TOT][NX2_TOT][NX1_TOT]; /*mean particle weight in the whole domain*/
  #if MULTIPLE_GHOSTS != YES
    double qz,qr,diagonal,sinth,costh;
  #endif
  static int first_call=1;
  #if (IMPOSE_BWALL==AS_DIFF)
    double t_diff_sec = t_diff*(UNIT_LENGTH/UNIT_VELOCITY); // time at which the diffusion has arrived! (seconds)
  #endif
  #if IMPOSE_TWALL
    double mu;
    double Twall_K = g_inputParam[TWALL]; // Wall temperature in Kelvin
  #endif
  #if IMPOSE_BWALL
    double unit_Mfield;
    double curr, Bwall, B_ghostwall; //Bwall, B_ghostwall in code units,
  #endif

  /*[Ema] g_time è: "The current integration time."(dalla docuementazione in Doxigen) */
  t_sec = g_time*(UNIT_LENGTH/UNIT_VELOCITY);

  #if IMPOSE_BWALL
    unit_Mfield = COMPUTE_UNIT_MFIELD(UNIT_VELOCITY, UNIT_DENSITY);
    // print1("\nCurrent from tab: %g", curr);
  #endif
  #if IMPOSE_BWALL==YES
    curr = current_from_time(t_sec);
    Bwall = BIOTSAV_GAUSS_A_CM(curr, RCAP)/unit_Mfield;
  #elif IMPOSE_BWALL==AS_DIFF
    curr = current_from_time(t_diff_sec);
    Bwall = BIOTSAV_GAUSS_A_CM(curr, RCAP)/unit_Mfield;
  #endif

  #ifdef DEBUG_ACCURATE_BCS
    int nv;
  #endif

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

    SetRBox_capWall(GetNghost());

    // /* Flatten the variables to conveniente values in points
    //    in internal boundary (except for "ghosts") */
    // /* WARNING!! IN CASE OF PRESSURE/TEMPERATURE TABLE INTERPOLATION ERROR, IT MIGHT
    //    BE ADVISABLE TO CHANGE THE VALUES HERE!!*/
    // KTOT_LOOP(k) {
    //   for (j=0; j<=j_cap_inter_end-1; j++){
    //     for (i=i_cap_inter_end+2; i<NX1_TOT; i++) {
    //       d->Vc[RHO][k][j][i] = 1e-2;
    //       d->Vc[iVZ][k][j][i] = 1e-3;
    //       d->Vc[iVR][k][j][i] = 1e-3;
    //       d->Vc[PRS][k][j][i] = 1e-3;
    //       d->Vc[iBPHI][k][j][i] = 1e-3;
    //     }
    //   }
    // }
    first_call = 0;
  }

  /**********************************
  ***********************************
  Actual setting of the boundary conditions
  ***********************************
  ***********************************/
  /* Maybe additional check for the runtime value of the boundary is useless,
     I keep it here as I want to be sure I can change the bc only by editing the
     pluto.ini file, and avoid editing also this file*/
  if (side == X1_END && RuntimeGet()->right_bound[IDIR] == USERDEF){
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
      FlipSign (X1_END, REFLECTIVE, vsign);
      ReflectiveBound (d->Vc[RHO], vsign[RHO], X1_END, CENTER);
      ReflectiveBound (d->Vc[iVZ], vsign[iVZ], X1_END, CENTER);
      ReflectiveBound (d->Vc[iVR], vsign[iVR], X1_END, CENTER);

      BOX_LOOP(box,k,j,i){
      #if IMPOSE_TWALL
        // Setting T
        setT( d, Twall_K, i, j, k);
      #else
        // I reflect pressure, to have no advection of energy through the capillary wall
        ReflectiveBound (d->Vc[PRS], vsign[PRS], X1_END, CENTER);
      #endif
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

  #if ACCURATE_BCS == YES

    #if (IMPOSE_BWALL || IMPOSE_TWALL || !defined(ELECTR_B_NEUM))
      #error IMPOSE_BWALL, IMPOSE_TWALL, not(ELECTR_B_NEUM), are not implemented in the accurate bcs
    #endif

    /***********************
    Capillary wall r=cost (CAP_WALL_INTERNAL)
    ************************/
    FlipSign (X1_END, REFLECTIVE, vsign);
    ReflectiveBoundCap (d->Vc, RHO, vsign[RHO], CAP_WALL_INTERNAL, CENTER);
    ReflectiveBoundCap (d->Vc, iVZ, vsign[iVZ], CAP_WALL_INTERNAL, CENTER);
    ReflectiveBoundCap (d->Vc, iVR, vsign[iVR], CAP_WALL_INTERNAL, CENTER);
    ReflectiveBoundCap (d->Vc, PRS, vsign[PRS], CAP_WALL_INTERNAL, CENTER);
    // [Err] Nutro qualche dubbio su questo.. dovrei riflettere B*r forse
    ReflectiveBoundCap (d->Vc, iBPHI, vsign[iBPHI], CAP_WALL_INTERNAL, CENTER);

    /***********************
    Capillary wall z=cost (CAP_WALL_EXTERNAL)
    ************************/
    FlipSign (X2_BEG, REFLECTIVE, vsign);
    ReflectiveBoundCap (d->Vc, RHO, vsign[RHO], CAP_WALL_EXTERNAL, CENTER);
    ReflectiveBoundCap (d->Vc, iVZ, vsign[iVZ], CAP_WALL_EXTERNAL, CENTER);
    ReflectiveBoundCap (d->Vc, iVR, vsign[iVR], CAP_WALL_EXTERNAL, CENTER);
    ReflectiveBoundCap (d->Vc, PRS, vsign[PRS], CAP_WALL_EXTERNAL, CENTER);
    ReflectiveBoundCap (d->Vc, iBPHI, vsign[iBPHI], CAP_WALL_EXTERNAL, CENTER);

    /*****************************
    Corner of wall
    ******************************/
    // BC  correction for IDIR
    FlipSign (X1_END, REFLECTIVE, vsign);
    ReflectiveBoundCap (d->Vc, RHO, vsign[RHO], CAP_WALL_CORNER_INTERNAL, CENTER);
    ReflectiveBoundCap (d->Vc, iVZ, vsign[iVZ], CAP_WALL_CORNER_INTERNAL, CENTER);
    ReflectiveBoundCap (d->Vc, iVR, vsign[iVR], CAP_WALL_CORNER_INTERNAL, CENTER);
    ReflectiveBoundCap (d->Vc, PRS, vsign[PRS], CAP_WALL_CORNER_INTERNAL, CENTER);
    // [Err] Nutro qualche dubbio su questo.. dovrei riflettere B*r forse
    ReflectiveBoundCap (d->Vc, iBPHI, vsign[iBPHI], CAP_WALL_CORNER_INTERNAL, CENTER);

    // BC correction for JDIR
    FlipSign (X2_BEG, REFLECTIVE, vsign);
    ReflectiveBoundCap (d->Vc, RHO, vsign[RHO], CAP_WALL_CORNER_EXTERNAL, CENTER);
    ReflectiveBoundCap (d->Vc, iVZ, vsign[iVZ], CAP_WALL_CORNER_EXTERNAL, CENTER);
    ReflectiveBoundCap (d->Vc, iVR, vsign[iVR], CAP_WALL_CORNER_EXTERNAL, CENTER);
    ReflectiveBoundCap (d->Vc, PRS, vsign[PRS], CAP_WALL_CORNER_EXTERNAL, CENTER);
    ReflectiveBoundCap (d->Vc, iBPHI, vsign[iBPHI], CAP_WALL_CORNER_EXTERNAL, CENTER);

    // Important! I must set to 0 all the corrections to the variables that I don't want to advance!
    SetNotEvolvedVar(iVPHI);
    SetNotEvolvedVar(iBR);
    SetNotEvolvedVar(iBZ);

    #ifdef DEBUG_ACCURATE_BCS
      for (nv=0, nv<NVAR; nv++;) {
        printf("\n var: %d", nv);
        printmat4d(d->Vc, NX2_TOT, NX1_TOT, nv, 0, -1, -1);
      }
      printcorr(d_correction[IDIR], "d_correction[IDIR]");
      printcorr(d_correction[JDIR], "d_correction[JDIR]");
    #endif

  #else

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
      // Temperature (or pressure)
      for (j=0; j<=j_cap_inter_end-1; j++) {
        #if IMPOSE_TWALL
          setT( d, Twall_K, i_cap_inter_end+1, j, k);
        #else
          // I reflect pressure, to have no advection of energy through the capillary wall
          d->Vc[PRS][k][j][i_cap_inter_end+1] = d->Vc[PRS][k][j][i_cap_inter_end];
        #endif
      }
      // Magnetic field on capillary wall (exclued electrode)
      for (j=0; j<j_elec_start; j++) {
        #if IMPOSE_BWALL
          B_ghostwall = (2*Bwall*rcap_real - d->Vc[iBPHI][k][j][i_cap_inter_end]*grid[IDIR].x_glob[i_cap_inter_end])/grid[IDIR].x_glob[i_cap_inter_end+1];
          d->Vc[iBPHI][k][j][i_cap_inter_end+1] = B_ghostwall;
        #else
          d->Vc[iBPHI][k][j][i_cap_inter_end+1] = d->Vc[iBPHI][k][j][i_cap_inter_end];
          // d->Vc[iBPHI][k][j][i_cap_inter_end+2] = d->Vc[iBPHI][k][j][i_cap_inter_end-1];
        #endif
      }
      // Magnetic field on electrode
      for (j=j_elec_start; j<=j_cap_inter_end; j++) {
        #if (IMPOSE_BWALL && !defined(ELECTR_B_NEUM))
          B_ghostwall = (2*Bwall*rcap_real - d->Vc[iBPHI][k][j][i_cap_inter_end]*grid[IDIR].x_glob[i_cap_inter_end])/grid[IDIR].x_glob[i_cap_inter_end+1];
          /* d->Vc[iBPHI][k][j][i_cap_inter_end+1] = B_ghostwall* \
                (1-(grid[1].x_glob[j]-grid[1].x_glob[j_elec_start])/ \
                (grid[1].x_glob[j_cap_inter_end]-grid[1].x_glob[j_elec_start])); */
          d->Vc[iBPHI][k][j][i_cap_inter_end+1] = B_ghostwall* \
                  (1 - (grid[JDIR].x_glob[j]-(zcap_real-dzcap_real))/dzcap );
        #else
          d->Vc[iBPHI][k][j][i_cap_inter_end+1] = d->Vc[iBPHI][k][j][i_cap_inter_end];
          // d->Vc[iBPHI][k][j][i_cap_inter_end+2] = d->Vc[iBPHI][k][j][i_cap_inter_end-1];
        #endif
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
        #if IMPOSE_BWALL
          d->Vc[iBPHI][k][j_cap_inter_end][i] = 0.0;
        #else
          d->Vc[iBPHI][k][j_cap_inter_end][i] = d->Vc[iBPHI][k][j_cap_inter_end+1][i];
        #endif
        // Temperature (or pressure)
        #if IMPOSE_TWALL
          setT( d, Twall_K, i, j_cap_inter_end, k);
        #else
          // I reflect pressure, to have no advection of energy through the wall
          d->Vc[PRS][k][j_cap_inter_end][i] = d->Vc[PRS][k][j_cap_inter_end+1][i];
        #endif
      }
    }
    /*********************
    Corner point (I write some possible algorithms to see which is better,
                  decomment the one you like!)
    *********************/
    #if MULTIPLE_GHOSTS != YES
      print1("\nBe careful! without multiple ghosts the bc for temperature(pressure) and mag.field is not quite ok, you should implement it more carefully");
      KTOT_LOOP(k) {
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

      /*********************/
      /* Define multiple ghosts on the corner wall cell! */
      /*********************/
      // To the normal corner cell I give the correct reflective boundary for
      // the r direction
      KTOT_LOOP(k) {
        d->Vc[RHO][k][j_cap_inter_end][i_cap_inter_end+1] = d->Vc[RHO][k][j_cap_inter_end][i_cap_inter_end];
        d->Vc[iVR][k][j_cap_inter_end][i_cap_inter_end+1] = -(d->Vc[iVR][k][j_cap_inter_end][i_cap_inter_end]);
        d->Vc[iVZ][k][j_cap_inter_end][i_cap_inter_end+1] = d->Vc[iVZ][k][j_cap_inter_end][i_cap_inter_end];
        // Magnetic field
        #if IMPOSE_BWALL
          d->Vc[iBPHI][k][j_cap_inter_end][i_cap_inter_end+1] = 0.0;
        #else
          d->Vc[iBPHI][k][j_cap_inter_end][i_cap_inter_end+1] = d->Vc[iBPHI][k][j_cap_inter_end][i_cap_inter_end];
        #endif
        // Temperature (or pressure)
        #if IMPOSE_TWALL
          setT( d, Twall_K, i_cap_inter_end+1, j_cap_inter_end, k);
        #else
          // I reflect pressure, to have no advection of energy through the wall
          d->Vc[PRS][k][j_cap_inter_end][i_cap_inter_end+1] = d->Vc[PRS][k][j_cap_inter_end][i_cap_inter_end];
        #endif
      }

      // I allocate memory for the correction structure
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
        #if IMPOSE_BWALL
          d_correction[0].Vc[iBZ][k] = 0.0;
          d_correction[0].Vc[iBPHI][k] = 0.0;
          d_correction[0].Vc[iBR][k] = 0.0;
        #else
          d_correction[0].Vc[iBZ][k] = d->Vc[iBZ][k][j_cap_inter_end][i_cap_inter_end];
          d_correction[0].Vc[iBPHI][k] = d->Vc[iBPHI][k][j_cap_inter_end][i_cap_inter_end];
          d_correction[0].Vc[iBR][k] = d->Vc[iBR][k][j_cap_inter_end][i_cap_inter_end];
        #endif

        #if IMPOSE_TWALL
          #if EOS==IDEAL
              #error double internal ghost not implemented for ideas eos
          #elif EOS==PVTE_LAW
              GetMu(Twall_K, d_correction[0].Vc[RHO][k], &mu); // No need to adim. T (GetMu takes T in K)
          #endif
          // I don't correct the Temperature, I set it same as in the normal
          // *d structure
          d_correction[0].Vc[PRS][k] = d_correction[0].Vc[RHO][k]*Twall_K / (KELVIN*mu);
        #else
          // I reflect pressure, to have no advection of energy through the wall
          d_correction[0].Vc[PRS][k] = d->Vc[PRS][k][j_cap_inter_end][i_cap_inter_end];
        #endif

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
        d_correction[1].Vc[iVPHI][k] = 0.0;
        d_correction[1].Vc[RHO][k] = d->Vc[RHO][k][j_cap_inter_end+1][i_cap_inter_end+1];
        #if IMPOSE_BWALL
          d_correction[1].Vc[iBZ][k] = 0.0;
          d_correction[1].Vc[iBPHI][k] = 0.0;
          d_correction[1].Vc[iBR][k] = 0.0;
        #else
          d_correction[1].Vc[iBZ][k] = d->Vc[iBZ][k][j_cap_inter_end+1][i_cap_inter_end+1];
          d_correction[1].Vc[iBPHI][k] = d->Vc[iBPHI][k][j_cap_inter_end+1][i_cap_inter_end+1];
          d_correction[1].Vc[iBR][k] = d->Vc[iBR][k][j_cap_inter_end+1][i_cap_inter_end+1];
        #endif
        #if IMPOSE_TWALL
          #if EOS==IDEAL
              #error double internal ghost not implemented for ideas eos
          #elif EOS==PVTE_LAW
              GetMu(Twall_K, d_correction[1].Vc[RHO][k], &mu);
          #endif
          // I don't correct the Temperature, I set it same as in the normal
          // *d structure
          d_correction[1].Vc[PRS][k] = d_correction[1].Vc[RHO][k]*Twall_K / (KELVIN*mu);
        #else
          // I reflect pressure, to have no advection of energy through the wall
          d_correction[1].Vc[PRS][k] = d->Vc[PRS][k][j_cap_inter_end+1][i_cap_inter_end+1];
        #endif
      }

      // No correction for the k direction
      d_correction[2].Npoints = 0;
    #endif
  #endif

    #ifdef FLATTEN_B_OUTCAP
      j_start_flatten = (j_cap_inter_end+NX2_TOT)/2;
      KDOM_LOOP(k) {
        for (j=j_start_flatten; j<=JEND; j++)
          for (i=IBEG; i<=IEND; i++)
            d->Vc[iBPHI][k][j][i] = 0.0;
      }
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

/*---------------------------------------------------------------------------------------------*/
/* Auxiliary function to set the mass density in the domain, according to few options          */
/*---------------------------------------------------------------------------------------------*/
void SetRhoAnalytic(double *rho, double x1, double x2, double x3, int mode) {

  if (mode == SMOOTH_COS2) {
    /****************************/
    /* cos² smoothing case      */
    /****************************/
    
    double dens0 = g_inputParam[DENS0]/UNIT_DENSITY;
    // Decay lenths in r and z, for setting density
    double decay_r = 0.1/UNIT_LENGTH;
    double decay_z = 0.5/UNIT_LENGTH;
    #ifdef DEBUG_EMA
      double rho_red_vac = 1;
    #else
      double rho_red_vac = 0.001; // Fraction of rho inside capillary, used to emumate vacuum
    #endif

    /* -----------------------------------------------------
        Zones not covered in the next lines (except for zone "Everywhere")
      ----------------------------------------------------- */
    *rho = rho_red_vac*dens0;
    /* -----------------------------------------------------
        Inside capillary, excluded near-electrode zone
      ----------------------------------------------------- */
    if (x2 < zcap-dzcap && x1 <= rcap) {
      *rho = dens0;
    }
    /* -----------------------------------------------------
        Inside capillary, in near-electrode zone
      ----------------------------------------------------- */
    if (zcap-dzcap <= x2 && x2 < zcap && x1 < rcap) {
      *rho = dens0;
    }
    /* ------------------------------------------------------
        Outside capillary, aligned with capillary (not in internal boundary)
      ------------------------------------------------------ */
    if (x2 > zcap && x1 <= rcap) {
      if (x2 < zcap+decay_z) {
        *rho = (1-rho_red_vac)*dens0;
        *rho *= cos(0.5*CONST_PI*(x2-zcap)/decay_z)*cos(0.5*CONST_PI*(x2-zcap)/decay_z);
        *rho += rho_red_vac*dens0;
      }
    }
    /* ------------------------------------------------------
        Outside capillary, above capillary  (not in internal boundary)
      ------------------------------------------------------ */
    if (x2 > zcap && x1 > rcap) {
      if (x1 < rcap+decay_r && x2 < zcap+decay_z) {
        *rho = (1-rho_red_vac)*dens0;
        *rho *= cos(0.5*CONST_PI*(x2-zcap)/decay_z)*cos(0.5*CONST_PI*(x2-zcap)/decay_z);
        *rho *= cos(0.5*CONST_PI*(x1-rcap)/decay_r)*cos(0.5*CONST_PI*(x1-rcap)/decay_r);
        *rho += rho_red_vac*dens0;
      }
    }

  } else if (mode == UNIFORM_FILL) {
    /******************************************/
    /* uniform filling case (with sharp edge) */
    /******************************************/

    double dens0 = g_inputParam[DENS0]/UNIT_DENSITY;
    #ifdef DEBUG_EMA
      double rho_red_vac = 1;
    #else
      double rho_red_vac = 0.001; // Fraction of rho inside capillary, used to emumate vacuum
    #endif
    /* -----------------------------------------------------
        Everywhere except inside capillary
      ----------------------------------------------------- */
    *rho = rho_red_vac*dens0;
    /* -----------------------------------------------------
        Inside capillary
      ----------------------------------------------------- */
    if (x2 < zcap && x1 <= rcap) {
      *rho = dens0;
    }
  }
}

// Disaster-case copy
// /*---------------------------------------------------------------------------------------------*/
// /* Auxiliary function to set the mass density in the domain, according to few options          */
// /*---------------------------------------------------------------------------------------------*/
// void SetRho(double *rho, double x1, double x2, double x3, int mode) {
//   double dens0 = g_inputParam[DENS0]/UNIT_DENSITY;

//   #ifdef DEBUG_EMA
//     double rho_red_vac = 1;
//   #else
//     double rho_red_vac = 0.001; // Fraction of rho inside capillary, used to emumate vacuum
//   #endif
//   #if DENS_INITIAL == SMOOTH_COS2
//     // Decay lenths in r and z, for setting density
//     double decay_r = 0.1/UNIT_LENGTH;
//     double decay_z = 0.5/UNIT_LENGTH;
//   #endif

//   /* -----------------------------------------------------
//       Zones not covered in the next lines (except for zone "Everywhere")
//     ----------------------------------------------------- */
//   *rho = rho_red_vac*dens0;
//   /* -----------------------------------------------------
//       Inside capillary, excluded near-electrode zone
//      ----------------------------------------------------- */
//   if (x2 < zcap-dzcap && x1 <= rcap) {
//     *rho = dens0;
//   }
//   /* -----------------------------------------------------
//       Inside capillary, in near-electrode zone
//      ----------------------------------------------------- */
//   if (zcap-dzcap <= x2 && x2 < zcap && x1 < rcap) {
//     *rho = dens0;
//   }
//   /* ------------------------------------------------------
//       Outside capillary, aligned with capillary (not in internal boundary)
//      ------------------------------------------------------ */
//   if (x2 > zcap && x1 <= rcap) {
//     #if DENS_INITIAL == SMOOTH_COS2
//       if (x2 < zcap+decay_z) {
//         *rho = (1-rho_red_vac)*dens0;
//         *rho *= cos(0.5*CONST_PI*(x2-zcap)/decay_z)*cos(0.5*CONST_PI*(x2-zcap)/decay_z);
//         *rho += rho_red_vac*dens0;
//       }
//     #elif DENS_INITIAL == UNIFORM_FILL
//       // Do nothing rho is already set outside 
//     #else
//       #error Choice for DENS_INITIAL not understood
//     #endif
//   }
//   /* ------------------------------------------------------
//       Outside capillary, above capillary  (not in internal boundary)
//      ------------------------------------------------------ */
//   if (x2 > zcap && x1 > rcap) {
//     // No field outside capillary
//     #if DENS_INITIAL == SMOOTH_COS2
//       if (x1 < rcap+decay_r && x2 < zcap+decay_z) {
//         *rho = (1-rho_red_vac)*dens0;
//         *rho *= cos(0.5*CONST_PI*(x2-zcap)/decay_z)*cos(0.5*CONST_PI*(x2-zcap)/decay_z);
//         *rho *= cos(0.5*CONST_PI*(x1-rcap)/decay_r)*cos(0.5*CONST_PI*(x1-rcap)/decay_r);
//         *rho += rho_red_vac*dens0;
//       }
//      #elif DENS_INITIAL == UNIFORM_FILL
//       // Do nothing rho is already set outside 
//     #else
//       #error Choice for DENS_INITIAL not understood
//     #endif
//   }
// }