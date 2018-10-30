#include "pluto.h"
#include "table_utilities.h"

#define REPRINT_ETA_TAB YES
#define REPRINT_KAPPA_TAB YES
#define ETA_TAB_SCRIPT "transport_tables_scripts/EtaTable_4pluto.py"
#define KAPPA_TAB_SCRIPT "transport_tables_scripts/KappaTable_4pluto.py"
#define ETA_TAB_FILE_NAME "eta.dat"
#define KAPPA_TAB_FILE_NAME "kappa.dat"

static Table2D eta_tab; /*    A 2D table containing pre-computed values of 
                              electr. resistivity stored at equally spaced node 
                              values of Log(T) and Log(rho) .*/
static Table2D kappa_tab; /*    A 2D table containing pre-computed values of 
                              therm. conductivity stored at equally spaced node 
                              values of Log(T) and Log(rho) .*/

/*****************************************************************************/
/* Function to build a table of electrical resistivity using a python script*/
/*****************************************************************************/
void MakeElecResistivityTable() {
  int i,j;
  double rho_min, rho_max, T_min, T_max;
  int N_rho, N_T;
  char table_finame[30] = ETA_TAB_FILE_NAME;
  char command[300];
  // char options[100] = " 800.0 30000.0 6 2.5e-11 2.7e-5 4";
  double **f;
  int logspacing;

  #if MAKE_ETA_TAB_FILE
    sprintf(command, "python3 %s %e %e %d %e %e %d %s", ETA_TAB_SCRIPT,
            (double)(T_TAB_MIN), (double)(T_TAB_MAX), (int) N_TAB_T,
            (double)(RHO_TAB_MIN), (double)(RHO_TAB_MAX), (int)(N_TAB_RHO),
            table_finame);
    system(command);
  #endif

  // Now I read the just made table
  ReadASCIITableSettings(table_finame, &logspacing,
                         &T_min, &T_max, &N_T, 
                         &rho_min, &rho_max, &N_rho);
  if (logspacing!=10) {
    print1("\n> MakeElecResistivityTable(): Error! Only logspacing 10 is supported!");
    QUIT_PLUTO(1);
  }

  print1 ("\n> MakeElecResistivityTable(): Generating table (%d x %d points)",
           N_T, N_rho);
  InitializeTable2D(&eta_tab,
                    T_min, T_max, N_T, 
                    rho_min, rho_max, N_rho);
  
  f = ARRAY_2D(N_rho, N_T, double);
  ReadASCIITableMatrix(table_finame, f, N_T, N_rho);

  for (j = 0; j < eta_tab.ny; j++)
    for (i = 0; i < eta_tab.nx; i++)
      eta_tab.f[j][i] = f[j][i];
  
  eta_tab.interpolation = LINEAR;

  FinalizeTable2D(&eta_tab);

  #if REPRINT_ETA_TAB
    ReprintTable(&eta_tab, table_finame);
  #endif

  FreeArray2D((void *)f);
}

/*************************************************************/
/* Function to get the Electrical res. from table            */
/*************************************************************/
double GetElecResisitivityFromTable(double rho, double T) {
  int    status;
  double eta;

  status = Table2DInterpolate(&eta_tab, T, rho, &eta);
  if (status != 0){
    print ("! GetElecResisitivityFromTable(): table interpolation failure (bound exceeded)\n");
    QUIT_PLUTO(1);
  }
  return eta;
}

/*****************************************************************************/
/* Function to build a table of thermal conductivity using a python script   */
/*****************************************************************************/
void MakeThermConductivityTable() {
  int i,j;
  double rho_min, rho_max, T_min, T_max;
  int N_rho, N_T;
  char table_finame[30] = KAPPA_TAB_FILE_NAME;
  char command[300];
  double **f;
  int logspacing;

  #if MAKE_KAPPA_TAB_FILE
    sprintf(command, "python3 %s %e %e %d %e %e %d %s", KAPPA_TAB_SCRIPT,
            (double)(T_TAB_MIN), (double)(T_TAB_MAX), (int) N_TAB_T,
            (double)(RHO_TAB_MIN), (double)(RHO_TAB_MAX), (int)(N_TAB_RHO),
            table_finame);
    system(command);
  #endif

  // Now I read the just made table
  ReadASCIITableSettings(table_finame, &logspacing,
                         &T_min, &T_max, &N_T, 
                         &rho_min, &rho_max, &N_rho);
  if (logspacing!=10) {
    print1("\n> MakeThermConductivityTable(): Error! Only logspacing 10 is supported!");
    QUIT_PLUTO(1);
  }

  print1 ("\n> MakeThermConductivityTable(): Generating table (%d x %d points)",
           N_T, N_rho);
  InitializeTable2D(&kappa_tab,
                    T_min, T_max, N_T, 
                    rho_min, rho_max, N_rho);
  
  f = ARRAY_2D(N_rho, N_T, double);
  ReadASCIITableMatrix(table_finame, f, N_T, N_rho);

  for (j = 0; j < kappa_tab.ny; j++)
    for (i = 0; i < kappa_tab.nx; i++)
      kappa_tab.f[j][i] = f[j][i];
  
  kappa_tab.interpolation = LINEAR;

  FinalizeTable2D(&kappa_tab);

  #if REPRINT_ETA_TAB
    ReprintTable(&kappa_tab, table_finame);
  #endif

  FreeArray2D((void *)f);
}

/*************************************************************/
/* Function to get the thermal cond. from table              */
/*************************************************************/
double GetThermConductivityFromTable(double rho, double T) {
  int    status;
  double kappa;

  status = Table2DInterpolate(&kappa_tab, T, rho, &kappa);
  if (status != 0){
    print ("! GetThermConductivityFromTable(): table interpolation failure (bound exceeded)\n");
    QUIT_PLUTO(1);
  }
  return kappa;
}