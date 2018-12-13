#ifndef TRANSPORT_TABLES_H
#define TRANSPORT_TABLES_H

void MakeElecResistivityTable();
int GetElecResisitivityFromTable(double rho, double T, double *eta);
void MakeThermConductivityTable();
int GetThermConductivityFromTable(double rho, double T, double *kappa);

#endif