#define  PHYSICS                 MHD
#define  DIMENSIONS              2
#define  COMPONENTS              3
#define  GEOMETRY                CYLINDRICAL
#define  BODY_FORCE              NO
#define  COOLING                 NO
#define  RECONSTRUCTION          LINEAR
#define  TIME_STEPPING           RK2
#define  DIMENSIONAL_SPLITTING   YES
#define  NTRACER                 0
#define  USER_DEF_PARAMETERS     10

/* -- physics dependent declarations -- */

#define  EOS                     PVTE_LAW
#define  ENTROPY_SWITCH          NO
#define  DIVB_CONTROL            NO
#define  BACKGROUND_FIELD        NO
#define  RESISTIVITY             ALTERNATING_DIRECTION_IMPLICIT
#define  THERMAL_CONDUCTION      ALTERNATING_DIRECTION_IMPLICIT
#define  VISCOSITY               NO
#define  ROTATING_FRAME          NO

/* -- user-defined parameters (labels) -- */

#define  ETAX_GAU                0
#define  ETAY_GAU                1
#define  ETAZ_GAU                2
#define  KAPPA_GAU               3
#define  DYN_VISC_GAU            4
#define  TWALL                   5
#define  T0                      6
#define  DENS0                   7
#define  VZ0                     8
#define  ALPHA_J                 9

/* [Beg] user-defined constants (do not change this line) */

#define  UNIT_DENSITY            1.67382e-8
#define  UNIT_LENGTH             1.0e-3
#define  UNIT_VELOCITY           1.0e6
#define  TC_SATURATED_FLUX       NO

/* [End] user-defined constants (do not change this line) */

/* -- supplementary constants (user editable) -- */ 

#define  INITIAL_SMOOTHING         NO
#define  WARNING_MESSAGES          YES
#define  PRINT_TO_FILE             YES
#define  INTERNAL_BOUNDARY         YES
#define  SHOCK_FLATTENING          NO
#define  CHAR_LIMITING             NO
#define  LIMITER                   DEFAULT
#define  ASSIGN_VECTOR_POTENTIAL   NO
#define  UPDATE_VECTOR_POTENTIAL   NO
#define  SHOW_TIME_STEPS           YES
#define  T_CUT_RHOE                800

/* ---------------------------------------------------- */
/* Additional constants */
#define  UNIT_ETA    (4*CONST_PI/(CONST_c*CONST_c)*UNIT_VELOCITY*UNIT_LENGTH)
#define  UNIT_KAPPA  (UNIT_DENSITY*UNIT_VELOCITY*UNIT_LENGTH*CONST_kB/CONST_mp)

/*  Ema's additional macros                            */
#define SPLIT_DIFF_ADV_ADV_DIFF
// #define PROFILE_GPROF_STOPSTEP  40
// #define FREEZE_FLUID

/* ---------------------------------------------------- */
/* ADI scheme settings                                  */
/*
Method for Thermal conduction and Resisitivity (when ADI is chosen), available choices:
  - FRACTIONAL_THETA
  - SPLIT_IMPLICIT
  - STRANG_LIE
  - DOUGLAS_RACHFORD
  - PEACEMAN_RACHFORD_MOD
  - STRANG
*/
#define METHOD_TC                  DOUGLAS_RACHFORD
#define METHOD_RES                 DOUGLAS_RACHFORD
/*
Number of sub-iterations in the whole "adi" scheme (at every sub-iteration the
conservative variables are updated and the kappa/eta re-evaluated)
(total adi-diffusive steps)
*/
#define NSUBS_ADI_TOT              4
/*
Set the period (in number of total adi-diffusive steps) for recomputation of discrete diffusive
*/
#define DIFF_OP_RECOMPUTE_PERIOD   3
/*
Number of sub-iterations for the thermal conduction scheme (the
conservative variables and kappa, are not updated between two iterations)
*/
#define NSUBS_TC                   30
/*
Number of sub-iterations for the magnetic diffusion scheme (the
conservative variables and eta, are not updated between two iterations)
*/
#define NSUBS_RES                  70

/*Theta value for Glowinsky's fractional theta method (a value in ]0,0.5[)*/
// #define FRACTIONAL_THETA_THETA_TC   0.3
// #define FRACTIONAL_THETA_THETA_RES  0.3
/* For a pseudo P-R algoritm: if FRACT==0.5 you have the usual P-R
  otherwise you unbalance the scheme towards the implcit or explicit part
  (keep it in ]0,1[)*/
// #define FRACT_TC                   0.4999999999999
// #define FRACT_RES            0.4999999999999
/*
To set the order of directions in the ADI scheme, allowed values: YES, NO, RANDOM, PERMUTE.
*/
#define FIRST_JDIR_THEN_IDIR       NO
// #define  TEST_ADI
#define JOULE_EFFECT_AND_MAG_ENG   (YES &&  RESISTIVITY==ALTERNATING_DIRECTION_IMPLICIT)
//Keep it YES for now. If YES: power flux is computed inside adi schemes (if NO, outside)
#define POW_INSIDE_ADI             YES
/* If YES: magnetic field related power sources (joule effect, div(B x (B x v))) are not
   computed inside diffusion/advection steps, but outisde, before and after the general
   Strang-splitted step. */
#define MAG_PS_OUTSIDE_SSTEP       NO

/*--------------------------------------------------------------------------*/
/* Other settings                                                           */

/* Macros to impose T (B) on walls also for advection (unphisical!)
   (if NO, conduction and B diffusion can be modeled only via ADI scheme) */
#define IMPOSE_TWALL               NO
#define IMPOSE_BWALL               NO
/* Decide whether the electrode must be set as a hom-Neumann boundary*/
#define ELECTR_B_NEUM
/* To set to 0 the mag field in a region outside capillary*/
// #define FLATTEN_B_OUTCAP
#define RUNTIMESET_CALL            AFTER_SETOUTPUT
#define CONE_LOW_TCKAPPA           (CONST_PI/4)
/* ---------------------------------------------------- */

/*---------------------------------------------------------------------------*/
/* Online analysis (Analysis() function)*/
#define EN_CONS_CHECK              YES
#define PRINT_TIME_INFO            YES
/*---------------------------------------------------------------------------*/

/* ---------------------------------------------------- */
/* Transport parameters, settings for computing them from tables*/
#define ETA_TABLE                  YES
#define KAPPA_TABLE                YES
#define MAKE_ETA_TAB_FILE          YES /* If YES, the ascii table file will be made with python script, */
#define MAKE_KAPPA_TAB_FILE        YES /* instead, if NO it is assumed that the file is already present*/
#define RHO_TAB_MIN                (2.5e-13)  /* You should never go below UNIT_DENSITY*1e-7 */
#define RHO_TAB_MAX                (2.5e-5)  /* You should never go hiher than UNIT_DENSITY*1e7 */
#define N_TAB_RHO                  50
#define T_TAB_MIN                  (0.8*T_CUT_RHOE)
#define T_TAB_MAX                  3.e5
#define N_TAB_T                    120
/* ---------------------------------------------------- */

/* ---------------------------------------------------- */
/*    CAPILLARY GEOMETRY SETTINGS                      */
#define RCAP                       0.05
#define DZCAP                      0.1 /*the electrodes are wide DZCAP cm*/
#define ZCAP                       1.6 /*the capillary is long 2*ZCAP cm and wide 2*RCAP cm*/
/* ---------------------------------------------------- */

/* ---------------------------------------------------- */
/* Parameters for hard/brute-force solutions for certain numerical problems 
- if rho is below RHO_VACUUM than the k is reduced a bit (see tc_kappa.c)
- if T is higher than T_MAX_HARD_RESET than it is put back to T_MAX_HARD_RESET
*/
#define RHO_VACUUM                 2.5e-10  /* You should never go below UNIT_DENSITY*1e-7 */
#define T_MAX_HARD_RESET           (T_TAB_MAX*0.99)  /* It must be below T_TAB_MAX and all other tabulated temperatures */
#define T_LIM_IEN                  1.0e5    /* Limit of temeperature above which EOS is modified (cv not constant anymore)*/
#define BETA_IEN                   2.5   /* When T>T_LIM_IEN the int.enrgy gets multiplied by (T/T_LIM_IEN)^BETA_IEN*/
#define RHO_MIN_HARD_RESET         5.0e-10  /* When RHO<RHO_MIN_HARD_RESET rho is set equal to RHO_MIN_HARD_RESET*/
/* ---------------------------------------------------- */

/* ---------------------------------------------------- */
/* CAPILLARY FILLING SETTINGS (initial condition) */
/* 'DENS_INITIAL' Available: SMOOTH_COS2 (cos² smothing ouside capillary),
                            UNIFORM_FILL (sharp transition with vacuum),
                            FROM_FILE    (from a bin file + grid file, named:  rho_ic.flt, grid_ic.out)*/
#define DENS_INITIAL               FROM_FILE
/* Useful only in case of DENS_INITIAL==FROM_FILE;
   if 1 replot from python the p and rho interpolated,
   if 0 it does nothing (keep 0 as default) */
#define REPLOT_P_RHO               0
/* ---------------------------------------------------- */

/* ---------------------------------------------------- */
#define MAX_LOGSIZE_MIB            4000.0
/* ---------------------------------------------------- */

/*------------------------------------------------------*/
/*     FOR DEBUG    */
// #define DEBUG_EMA
// #define DEBUG_BCS
// #define DEBUG_BUILDIJ
// #define DEBUG_TNEGATIVE
#define WARN_CTP_FAIL  NO
#define WARN_ERR_COMP_TEMP NO
