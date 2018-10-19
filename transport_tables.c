#include "pluto.h"
#include "table_utilities.h"

#define REPRINT_ETA_TAB YES
#define ETA_TAB_SCRIPT "transport_tables_scripts/EtaTable_4pluto.py"

static Table2D eta_tab; /**< A 2D table containing pre-computed values of 
                              electr. resistivity stored at equally spaced node 
                              values of Log(T) and Log(rho) .*/

void MakeElecResistivityTable() {
  int i,j;
  double rho_min, rho_max, T_min, T_max;
  int N_rho, N_T;
  char table_finame[30] = "eta.dat";
  char command[200] = "python3 ";
  char options[100] = " 800.0 30000.0 6 2.5e-11 2.7e-5 4";
  double **f;
  int logspacing;

  // I make the table with a python3 script
  strcat(command, ETA_TAB_SCRIPT);
  strcat(command, options);
  strcat(command, " ");
  strcat(command, table_finame);
  print1("%s", command);
  system(command);

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