# -*- coding: utf-8 -*-
import numpy as np
import sys
import datetime as dt
import importlib
import PlasmaPar_asDevoto as ppdv
import ionization as iz
import constantsGAU_ema as gau

importlib.reload(ppdv)
importlib.reload(gau)
importlib.reload(iz)

def write_ASCIITable_4pluto(tab_fi, logspacing,
                            T_min, T_max, N_T,
                            rho_min, rho_max, N_rho,
                            eta_Dev):
#    print(type(logspacing))

    #    ftab = open(tab_fi, mode='w')
#    
#    ftab.write('# comment line')
#    ftab.write('\n# another comment line')
#    ftab.write('\n{:d}'.format(logspacing))
#    ftab.write('\n{:e}'.format(T_min))
#    ftab.write('\n{:e}'.format(T_max))
#    ftab.write('\n{:d}'.format(N_T))
#    ftab.write('\n{:e}'.format(rho_min))
#    ftab.write('\n{:e}'.format(rho_max))
#    ftab.write('\n{:d}'.format(N_rho))
#    ftab.write('\n{}'.format(eta_Dev))
#    
#    ftab.close()
    
    comm1 = '# ' + str(dt.datetime.now())
    
    comm2 = '\n# numb. meaning: logspacing(integer, only allowed:10),'
    comm2 += 'min(T),max(T),numb. T points,'
    comm2 += 'min(rho),max(rho),numb. rho points,'
    comm2 += 'matrix(T,rho): columns(or fastest running index) span rho'
    
    tab_header = comm1 + comm2
    tab_header += '\n{:d}'.format(logspacing)
    tab_header += '\n{:e}\n{:e}\n{:d}'.format(T_min, T_max, N_T)
    tab_header += '\n{:e}\n{:e}\n{:d}'.format(rho_min, rho_max, N_rho)    
    np.savetxt(tab_fi, eta_Dev, fmt = '%.10e', delimiter = ' ', newline = '\n',
               header = tab_header, comments = '')
    
    return 0

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
eta_Dev = np.zeros((N_T, N_rho))
x = np.zeros((N_T, N_rho))
y = np.zeros((N_T, N_rho))
for rr in range(N_rho):
    for tt in range(N_T):
        # I use this function but I neglect later dissociation (in ppdv.elRes_norm)
        x,y = iz.ionizDissSaha(rho[rr], kT[tt])
        # I put a lower limit on ys
        y = max(ioniz_min, y)
        eta_Dev[tt,rr] = ppdv.elRes_norm(y, rho[rr], kT[tt])

# <codecell>
# Print eta table to ASCII file
if (write_ASCIITable_4pluto(tab_fi, 10,
                              T_min, T_max, N_T,
                              rho_min, rho_max, N_rho,
                              eta_Dev) == 0):
    print('Eta table successfully created!')
else:
    raise ValueError('Eta table successfully created!')