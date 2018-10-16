# Notes on this simulation

## *PVTE_LAW*
PVTE_LAW is being used, so (as userguide says at page 69), it is suggested to use RK time stepping with tvdlf, hll or hllc Riemann solvers. Entropy switch is not supported. And algorithms requiring charactersitic decomposition are not supported (I am not sure what does this mean, Alberto said ti refers to some scheme advancing the solution in time (remember that the characteristics are waves! So it is related to the way one uses to solve Riemann problems))(can I cannot use some implementations of the function "States()"? I don't know which uses characterisic decomposition: the one in limo3_states.c and maybe the one in plm_states.c (when CHAR_LIMITING == YES) (I am not sure whether I can use it with CHAR_LIMITING == NO) )

### Temperature diffusion with *PVTE_LAW*
Also I updated the temperature computation only in the time advancing (and other) routines which work with the rk stepping. So, despite what said in the userguide, you are not simply adviced to use rk time stepping, but *forced to*!
__I should also check thet I improved correctly the temperature computation also in the routines (if needed!) that compute it for super time stepping. Otherwise I shoul use explicit time stepping.__

## *Internal boundary*
When using internal boundary you must not alter the staggered magnetic field values, which means you cannot use (if I understand currectly) STAGGERED_MHD.

## Change in the definition/declaration of the UpdateStage() function and in the "prototypes.h" file.
I removed the "const" specifier in front of the "Data \*d" entry in the UpdateStage() function.
This is useful for modifying the Data structure to apply multiple ghosts on a single internal cell.
To do so, I also had to change the function prototype in prototypes.h.
Thus I edited the file prototypes.h in the original Source folder, since copying it into my simulation directory doesn't really work, the compiler gives the precedence always to the one in the origin src folder!

## To summarize:
**Do not use:**
+ Entropy switch
+ Staggered MHD
+ States() function inside limo3_states.c
+ States() (maybe) function inside plm_states.c with CHAR_LIMITING==YES

**Do use:**
+ RK time stepping

**Better use:**
+ tvdlf, hll or hllc Riemann solvers
+ States() function inside ppm_states.c with CHAR_LIMITING==NO

**Check**
+ Tempertature computation in case of super-time-stepping is correct?

## For profiling the code:
- Insert the lines:
`` CFLAGS += -pg ``
`` LDFLAGS += -pg ``
- in the local_make file, and inside definitions.h define a macro:
`` #define PROFILE_GPROF   pluto_time_step_where_you_want_to_stop ``
- then compile and run pluto as usual;
- then call ``gprof pluto > filename_where_to_store_gprof_output``

## For restarting simulations (you need the bash script cppsim, saved for instance in ~/bin and sourced from ~/bashrc or similar...)
Type in terminal:
``cppsim N ORIG DEST``
where N is the number of dbl data dump from which you want to restart, orig is the folder of the simulation you want to restart, DEST is the name you give to the restarted simulation (the folder will be created ex novo, _don't use the name of an existing folder!_).

