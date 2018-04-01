#include "pluto.h"
#include "current_table.h"

void WriteCurrForDebug (char *file_basename, Data *d,
                        Grid *grid, int dir){
double unit_J, unit_Mfield;
char *J_name, *file_name;
Data_Arr J;
Grid grid_curr;
FILE *fvtk;
//This variable keeps trace of the number of times this function has been called
static int curr_debug_idx=-1;

curr_debug_idx++;

if (dir == IDIR) {
// Create shifted grid accordingly to the direction we're in

}
if (dir == JDIR) {
// Create shifted grid accordingly to the direction we're in

}

// Assign current density to Data_Arr J;


// Compute the unit of measure
unit_Mfield = COMPUTE_UNIT_MFIELD(UNIT_VELOCITY, UNIT_DENSITY);
unit_J = CONST_c/(4*CONST_PI)*unit_Mfield/UNIT_LENGTH;

// Open file, using curr_debug_idx in the name to be sure you do not overwrite other files
file_name = malloc(strlen(file_basename) + 5 + 4 + 1);
sprintf(file_name,"%s.%d.vtk", file_basename, curr_debug_idx);

fvtk  = OpenBinaryFile(file_name, SZ_Float_Vect, "w");

// Write J to vtk
WriteVTK_Header (fvtk, &grid_curr); NO! rifare anche questa funzione perchè lei usa DOM_LOOP e altr macro, che per me non vanno bene
WriteVTK_Vector (fvtk, J, unit_J, J_name, &grid_curr);

fclose(fvtk);
}