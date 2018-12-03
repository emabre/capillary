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
``cppsim -restart N ORIG DEST``
where N is the number of dbl data dump from which you want to restart, orig is the folder of the simulation you want to restart, DEST is the name you give to the restarted simulation (the folder will be created ex novo, _don't use the name of an existing folder!_). Remember then to restart the simulation with the usual pluto's command (in terminal):
``./pluto -restart N``
and maybe add nohup at the beginning of the line.

## Using a table (e.g. for eta (resistivity) as function of rho, T)
### How pluto makes a table
E.g. for the internal energy:
  - calls: ``InitializeTable2D()`` to initialize the table
  - fills the internal energy values in the table looping: ``rhoe_tab.f[j][i] = InternalEnergyFunc(v,T);``
  - it sets something related to later interpolation with: ``rhoe_tab.interpolation = SPLINE1;``
  - it computes the *cubic spline coefficients*, with some lines of code (Actually, I guess I cannnot copy-paste them, mutatis mutandis, since I would need also an analytical expression of the derivative of the quantity I tabulate. But if I use LINEAR interpolation, then I don't need this part (I hope!))
  - it calls ``FinalizeTable2D()`` which computes the differences between table entries which are adiacent in x or y (I guess it might be useful later either for faster differentiation or for interpolation..)
### How pluto reads the desired values from a table
E.g. for the internal energy:
``status = Table2DInterpolate(&rhoe_tab, T, rho, &rhoe);``
status: 0 if everything works well,
T: temperature (scalar, double),
rho: mass density (scalar, double),
&rhoe: pointer to internal energy (scalar, double*)
&rhoe_tab: pointer to table (it's defined as: ``static Table2D rhoe_tab;``)

### How to make a table using Python
If you go with the terminal inside the folder *transport_tables_scripts* and type (for instance):
``python3 EtaTable_4pluto.py 800.0 30000.0 6 2.5e-11 2.7e-5 4 eta.dat``
You build an ASCII file containing a table (and its required header) suitable to be read from PLUTO, named *eta.tab*, for:
+ $T \in [800, 30000]$, with 6 log-spaced points,
+ $\rho \in [2.5e-11,2.5e-5]\mathrm{g/cm}^3$, with 4 log-spaced points,

## Doing tests on this simulation environment:
You can test that the plasma heating due to ohmic dissipation is correct by checking the temperature change of the plasma.
+ Set the thermal diffusion off (or set a very low thermal conductivity, e.g. 1e2 in gaussian units);
+ Set the advection off (``#define FREEZE_FLUID`` in __definitions.h__);
+ Set the current to a fixed value (e.g. 100A);
+ Set the initial temperature of the gas not close to any non linearity due to ionization (be careful that non linearity in heat capacity due to ionization starts very early, roughly 5000K; also be careful with dissociation, if you implement it).
+ Run the simulation for the time you want
+ Compare the temperature increase with the one you foresee theoretically.
To foresee the temperature increase __for a case where the hydrogen is fully dissociated but no ioniation is present__ use this formula:
$$ \Delta T = 7.3546\cdot 10^9 \frac{\eta I_\mathrm{A}^2 t}{R^4 \rho} $$,
where all the quantities are in cgs gaussian units, except for the current which is in Ampere.
For instance, for the case: $\eta=10^{-14}\mathrm{s}$, $I_\mathrm{A}=100\mathrm{A}$, $t=1\mathrm{ns}$, $R=0.05\mathrm{cm}$, $\rho=2.5\cdot10^{-7}\mathrm{g}/\mathrm{cm}^3$, you get $\Delta T = 470.69\mathrm{K}$.
