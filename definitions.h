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

/* ---------------------------------------------------- */
/*  Ema's additional macros                            */
// #define  FREEZE_FLUID
// #define REVERSE_ADI_DIRS
// #define  TEST_ADI
#define JOULE_EFFECT_AND_MAG_ENG YES
/* Macros to impose T (B) on walls also for advection (unphisical!)
   (if NO, conduction and B diffusion can be modeled only via ADI scheme) */
#define IMPOSE_TWALL NO
#define IMPOSE_BWALL NO
/* ---------------------------------------------------- */

/* ---------------------------------------------------- */
/*    CAPILLARY GEOMETRY SETTINGS                      */
#define RCAP 0.05
#define DZCAP 0.5 /*the electrodes are wide DZCAP cm*/
#define ZCAP 1.5 /*the capillary is long 2*ZCAP cm and wide 2*RCAP cm*/
/* ---------------------------------------------------- */
