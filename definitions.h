#define  PHYSICS                 MHD
#define  DIMENSIONS              2
#define  COMPONENTS              3
#define  GEOMETRY                CYLINDRICAL
#define  BODY_FORCE              NO
#define  COOLING                 NO
#define  RECONSTRUCTION          PARABOLIC
#define  TIME_STEPPING           RK3
#define  DIMENSIONAL_SPLITTING   YES
#define  NTRACER                 0
#define  USER_DEF_PARAMETERS     9

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
#define  KAPPA_GAUBOB            3
#define  TWALL                   4
#define  T0                      5
#define  DENS0                   6
#define  VZ0                     7
#define  ALPHA_J                 8

/* [Beg] user-defined constants (do not change this line) */

#define  UNIT_DENSITY            1.67382e-6
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

/* ---------------------------------------------------- */
/* Additional constants */
#define  UNIT_ETA    (4*CONST_PI/(CONST_c*CONST_c)*UNIT_VELOCITY*UNIT_LENGTH)
#define  UNIT_KAPPA  (UNIT_DENSITY*UNIT_VELOCITY*UNIT_LENGTH*CONST_kB/CONST_mp)

/*  Ema's additional macros                            */
// #define PROFILE_GPROF_STOPSTEP      8
// #define FREEZE_FLUID
/* Either give value in ]0,0.5[
   or comment if you prefer to use Peaceman-Rachford scheme */
// #define FRACTIONAL_THETA           0.3
// #define SPLIT_IMPLICIT
/* For a pseudo P-R algoritm: if FRACT==0.5 you have the usual P-R
  otherwise you unbalance the scheme towards the implcit or explicit part
  (keep it in ]0,1[). If you do not define it it will be set to 0.5*/
#define FRACT                      0.5
// to set the order of the ADI scheme, allowed values: YES, NO, RANDOM
#define FIRST_JDIR_THEN_IDIR       NO
// #define  TEST_ADI
#define JOULE_EFFECT_AND_MAG_ENG   (YES &&  RESISTIVITY==ALTERNATING_DIRECTION_IMPLICIT)
//Keep it YES for now. If YES: power flux is computed inside adi schemes (if NO, outside)
#define POW_INSIDE_ADI             YES
/* Macros to impose T (B) on walls also for advection (unphisical!)
   (if NO, conduction and B diffusion can be modeled only via ADI scheme) */
#define IMPOSE_TWALL               YES
#define IMPOSE_BWALL               AS_DIFF
/* Number of subcycles performed by ADI scheme*/
#define NSUBS_ADI                  20
/* Decide whether the electrode must be set as a hom-Neumann boundary*/
#define ELECTR_NEUM
/* To set to 0 the mag field in a region outside capillary*/
// #define FLATTEN_B_OUTCAP
// #define PRESUBS_RES
/* ---------------------------------------------------- */

/* ---------------------------------------------------- */
/*    CAPILLARY GEOMETRY SETTINGS                      */
#define RCAP                       0.05
#define DZCAP                      0.3 /*the electrodes are wide DZCAP cm*/
#define ZCAP                       1.5 /*the capillary is long 2*ZCAP cm and wide 2*RCAP cm*/
/* ---------------------------------------------------- */

/* ---------------------------------------------------- */
/*     ADDITTIONAL OUTPUT                               */
#define WRITE_J                     YES