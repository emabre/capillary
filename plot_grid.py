# To plot pluto's grid written in the usual pluto's file grid.out

import numpy as np
import matplotlib.pyplot as plt
import os

# <codecell>
# Options
plt.ion()
plt.close("all")
output_dir = "out"
grid_fi = os.path.join(output_dir, "grid.out")

# <codecell>
# read grid file
gf = open(grid_fi, 'r')
# ii=0
# for line in gf:
#     ii = ii +1 
#     print(line)
#     if ii>20:
#         break
ll = 0
Np = []
grid = []
for line in gf:
    ll += 1
    line_split = line.split()
    if line_split[0][0]=="#":
        # It's just a comment
        next
    elif len(line_split)==1:
        try:
            Np.append(int(line_split[0]))
            current_dim = len(Np)
            grid.append([])
            print("I am reading dimension {}".format(current_dim))
        except ValueError:
            raise ValueError("Unexpected single (non-int)number in line {}".format(ll))
    elif len(line_split) == 3:
        grid[current_dim-1].append([float(line_split[1]), float(line_split[2])])
    else:
        raise ValueError("Unexpected number of columns ad line {}".format(ll))

if len(Np)!=3:
    raise ValueError("Strange grid dimension: {}".format( len(Np)))

x1 = np.array(grid[0])
x2 = np.array(grid[1])
x3 = np.array(grid[2])

## Read unit length
df = open("definitions.h")
for line in df:
    if len(line.split()) == 3:
        if line.split()[0]=="#define" and line.split()[1]=="UNIT_LENGTH":
            UNIT_LENGTH = float(line.split()[2])

# <codecell>
# Some computations
# Ratio between adiacent Dx's
ratio_x1 = (x1[1:,1]-x1[1:,0])/(x1[:-1,1]-x1[:-1,0])
ratio_x2 = (x2[1:,1]-x2[1:,0])/(x2[:-1,1]-x2[:-1,0])
# Meshlines
x1min = x1[0,0]
x1max = x1[-1,1]
x2min = x2[0,0]
x2max = x2[-1,1]
lines_x1 = [x1[0,0]] + [x1[ii,1] for ii in range(Np[0])]
lines_x1 = np.array(lines_x1)
lines_x2 = [x2[0,0]] + [x2[ii,1] for ii in range(Np[1])]
lines_x2 = np.array(lines_x2)

# <codecell>
# Plot
if x3.shape[0]!=1:
    raise ValueError("I cannot plot as the grid is not 2D")

fig1, ax1 = plt.subplots()
fig2, ax2 = plt.subplots()

for ii in range(Np[0]):
    ax1.plot([lines_x1[ii], lines_x1[ii]], [lines_x2[0],lines_x2[-1]], 'k', lw='0.5')
for ii in range(Np[1]):
    ax1.plot([lines_x1[0],lines_x1[-1]], [lines_x2[ii], lines_x2[ii]], 'k', lw='0.5')
ax1.set_xlabel("x1 / UNIT_LENGTH")
ax1.set_ylabel("x2 / UNIT_LENGTH")
# plt.axis('equal')

## Plot in cgs
lines_x1_cgs = UNIT_LENGTH*lines_x1
lines_x2_cgs = UNIT_LENGTH*lines_x2

for ii in range(Np[0]):
    ax2.plot([lines_x1_cgs[ii], lines_x1_cgs[ii]],
               [lines_x2_cgs[0],lines_x2_cgs[-1]],
               'k', lw='0.5')
for ii in range(Np[1]):
    ax2.plot([lines_x1_cgs[0],lines_x1_cgs[-1]],
               [lines_x2_cgs[ii], lines_x2_cgs[ii]],
               'k', lw='0.5')
ax2.set_xlabel("x1 / cm")
ax2.set_ylabel("x2 / cm")
# plt.axis('equal')

# Plot of the ratios of adiacent DeltaX
fig_r, ax_r = plt.subplots()
ax_r.plot(ratio_x1, label=r'$\Delta x1_i/\Delta x1_{i+1}$')
ax_r.plot(ratio_x2, label=r'$\Delta x2_i/\Delta x2_{i+1}$')
ax_r.legend()

plt.show()
