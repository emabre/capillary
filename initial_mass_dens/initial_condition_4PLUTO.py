import numpy as np
import os
# Note: interpolate (as others) need to be explicitly imported
import scipy.interpolate as interp
import matplotlib.pyplot as plt
import struct
import sys

plt.close("all")

# <codecell>
# Main settings
if len(sys.argv) != 13: # 8 = 7(provided args) + 1(script name(sys.argv[0]))
    raise ValueError("Input args must be 10!")

# Everything in cgs units
UNIT_LENGTH = float(sys.argv[1])
UNIT_DENSITY = float(sys.argv[2])
grid_Nr = int(sys.argv[3])
grid_Nz = int(sys.argv[4])
grid_rbeg = float(sys.argv[5])
grid_rend = float(sys.argv[6])
grid_zbeg = float(sys.argv[7])
grid_zend = float(sys.argv[8])
min_rho_cgs = float(sys.argv[9])
replot_p_rho = bool(int(sys.argv[10]))  # convert to int and then to bool otherwise it will always become True
rho_fi = str(sys.argv[11])
grid_finame = str(sys.argv[12])

# <codecell>
# Other settings
case_dir = 'initial_mass_dens'

T_finame = 'T.raw'
p_finame = 'p.raw'

# Tolerance to declare zero (and thus neglect) the a certain length
pos_zero_toll = 1e-15

# How to fill regions outside from provided data.
# "min" -> put minimum value in domain (or minimum allowed)
# no other possibility implemented up to now
fill_outside = "min"

# <codecell>
# Read data
p_fi = os.path.join(case_dir, p_finame)
T_fi = os.path.join(case_dir, T_finame)

# I read the data
p = np.loadtxt(p_fi, skiprows=2)
T = np.loadtxt(T_fi, skiprows=2)

if np.any(np.abs(T[:,(0,1)]-p[:,(0,1)])>pos_zero_toll):
    raise ValueError("p and T refer to different points, cannot continue!")

if np.all(p[:,2]<pos_zero_toll):
    p = np.delete(p, 2, 1)
else:
    raise ValueError("I cannot remove the third dimension from p, as it seems not zero")
if np.all(T[:,2]<pos_zero_toll):
    T = np.delete(T, 2, 1)
else:
    raise ValueError("I cannot remove the third dimension from T, as it seems not zero")

# <codecell>
# I elaborate the data
mp = 1.6726e-25 # proton mass (g)
kB = 1.3807e-16 # Boltzm. const. (erg/K)

# I convert presure to dyne/cm² to have gcs units
p[:,2] *= 10.0
# I convert x,y to cm to have cgs units
p[:,0] *= 1e2
p[:,1] *= 1e2
T[:,0] *= 1e2
T[:,1] *= 1e2

# I compute rho (mass density)
print("[initial_condition_4PLUTO]I assume this is molecular Hydrogen (no dissociation, no ionization)")
# Mass density in g/cm³
rho = np.zeros(p.shape)
rho[:,(0,1)] = p[:,(0,1)]
rho[:,2] = 2*mp*p[:,2]/(kB*T[:,2])

# <codecell>
# Interpolation
# Compute cell vertexes
z = np.linspace(grid_zbeg, grid_zend, grid_Nz+1)
r = np.linspace(grid_rbeg, grid_rend, grid_Nr+1)
r_mg, z_mg = np.meshgrid(r,z)

# Compute cell centers
z_cc = 0.5*(z[1:]+z[:-1])
r_cc = 0.5*(r[1:]+r[:-1])
r_cc_mg, z_cc_mg = np.meshgrid(r_cc,z_cc)

if fill_outside=="min":
    fill_rho = np.min(rho[:,2])
else:
    raise ValueError("Option for filling regions outside provided domain not understood")
rho_i = interp.griddata(rho[:,(0,1)], rho[:,2],
                        (z_cc_mg.flatten(),r_cc_mg.flatten()),
                        fill_value = fill_rho)
# Where rho_i is lower than the tolerance, I set it equal to the tolerance
rho_i[rho_i<min_rho_cgs] = min_rho_cgs

if replot_p_rho:
    if fill_outside=="min":
        fill_p = np.min(p[:,2])
    else:
        raise ValueError("Option for filling regions outside provided domain not understood")
    p_i = interp.griddata(p[:,(0,1)], p[:,2],
                          (z_cc_mg.flatten(),r_cc_mg.flatten()),
                          fill_value = fill_p)
    plt.figure()
    plt.pcolormesh(z_mg,r_mg, p_i.reshape(z_cc_mg.shape))
    plt.title("p_i")
    plt.xlabel("z / cm")
    plt.ylabel("r / cm")
    plt.colorbar()
    plt.figure()
    plt.pcolormesh(z_mg,r_mg, rho_i.reshape(z_cc_mg.shape))
    plt.title("rho_i")
    plt.xlabel("z / cm")
    plt.ylabel("r / cm")
    plt.colorbar()
    plt.ioff()
    plt.show()


# <codecell>
# Write grid file

# r and z in pluto convention
r_pl = r / UNIT_LENGTH
z_pl = z / UNIT_LENGTH

grid_f = open(grid_finame, 'w')
grid_f.write('# GEOMETRY: CYLINDRICAL')

grid_f.write('\n{:d}'.format(grid_Nr)) # ARE YOU SURE THE ORDER IS CORRECT? (first r then z or the opposite?)
for ii in range(grid_Nr):
    grid_f.write('\n {:d}   {:e}   {:e}'.format(ii+1, r_pl[ii], r_pl[ii+1]))

grid_f.write('\n{:d}'.format(grid_Nz))
for ii in range(grid_Nz):
    grid_f.write('\n {:d}   {:e}   {:e}'.format(ii+1, z_pl[ii], z_pl[ii+1]))

grid_f.write('\n{:d}'.format(1))
grid_f.write('\n {:d}   {:e}   {:e}'.format(1, 0.0, 1.0))

grid_f.close()

# <codecell>
# Write binary file
rho_i_pl = rho_i/UNIT_DENSITY
# '>' means big-endian, '<' is little-endian
rho_bin = struct.pack('>'+'f'*len(rho_i_pl), *rho_i_pl)
with open(rho_fi, 'wb') as rho_f:
    rho_f.write(rho_bin)
