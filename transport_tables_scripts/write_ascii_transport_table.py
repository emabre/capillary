# -*- coding: utf-8 -*-
import datetime as dt
import numpy as np

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
