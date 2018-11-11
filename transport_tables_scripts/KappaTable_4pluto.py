# -*- coding: utf-8 -*-
import numpy as np
import sys
import importlib
import PlasmaPar_asDevoto as ppdv
import ionization as iz
import constantsGAU_ema as gau
import write_ascii_transport_table as wa

importlib.reload(ppdv)
importlib.reload(gau)
importlib.reload(iz)
importlib.reload(wa)

# <codecell>
if len(sys.argv) != 8: # 8 = 7(provided args) + 1(script name(sys.argv[0])) 
    raise ValueError("Input args must be 7!")

T_min = float(sys.argv[1])
T_max = float(sys.argv[2])
N_T = int(sys.argv[3])
rho_min = float(sys.argv[4])
rho_max = float(sys.argv[5])
N_rho = int(sys.argv[6])
tab_fi = str(sys.argv[7])

# <codecell>
# I make the log spacing
T = np.logspace(np.log10(T_min), np.log10(T_max), N_T)
rho = np.logspace(np.log10(rho_min), np.log10(rho_max), N_rho)

kT = T*gau.kB
ioniz_min = 1e-10

# I compute eta
kappa_Dev = np.zeros((N_T, N_rho))
x = np.zeros((N_T, N_rho))
y = np.zeros((N_T, N_rho))
for rr in range(N_rho):
    for tt in range(N_T):
        # I use this function but I neglect later dissociation (in ppdv.elRes_norm)
        x,y = iz.ionizDissSaha(rho[rr], kT[tt])
        # I put a lower limit on ys
        y = max(ioniz_min, y)
        # Hydrogen thermal conductivity (traslative + reactive, no internal)
        # assuming full dissociation and partial ionization. No B field.
        kappa_Dev[tt,rr] = ppdv.thermCond_tot_norm(y, rho[rr], kT[tt])

# <codecell>
# Print eta table to ASCII file
if (wa.write_ASCIITable_4pluto(tab_fi, 10,
                              T_min, T_max, N_T,
                              rho_min, rho_max, N_rho,
                              kappa_Dev) == 0):
    print('Kappa table successfully created!')
else:
    raise ValueError('Kappa table NOT successfully created!')