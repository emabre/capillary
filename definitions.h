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
#define  USER_DEF_PARAMETERS     4

/* -- physics dependent declarations -- */

#define  EOS                     PVTE_LAW
#define  ENTROPY_SWITCH          NO
#define  DIVB_CONTROL            NO
#define  BACKGROUND_FIELD        NO
#define  RESISTIVITY             SUPER_TIME_STEPPING
#define  THERMAL_CONDUCTION      SUPER_TIME_STEPPING
#define  VISCOSITY               NO
#define  ROTATING_FRAME          NO

/* -- user-defined parameters (labels) -- */

#define  ETAX_GAU                0
#define  ETAY_GAU                1
#define  ETAZ_GAU                2
#define  KAPPA_GAUBOB            3  // "..GAUBOB" because it is in gau units with Bob convention(kB not included)

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
#define  VTK_VECTOR_DUMP           YES