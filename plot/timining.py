#!/usr/bin/env python3
# To plot pluto's simulation timings, written by Ema in out/times.dat

import numpy as np
import matplotlib.pyplot as plt
import os

# <codecell>
output_dir = "out"
times_fi = os.path.join(output_dir, "times.dat")

# Number of calls to Analysis() over which to compute the elapsed time / simulated time
# minimum is 1
an_calls_avg = 3

# <codecell>
if an_calls_avg<1:
    raise ValueError('an_calls_average must be at least 1')
# read time file
times = np.loadtxt(times_fi, comments='#', skiprows=2)
        
# Real elapsed time (hours)
el_time = times[:,3]
el_time_h = el_time/3600.0
# Simulated time (seconds)
t_sim = times[:,1]
# number of calls to Analysis() function in PLUTO
an_calls = times[:,0]
# time-step (seconds)
dt_sim = times[:,2]

if times.shape[0]>an_calls_avg:
    el_time_ov_sim_time_computed = True
    # I compute elapsed time / simulated time
    el_time_ov_sim_time = (el_time[an_calls_avg:]-el_time[:-an_calls_avg])/(t_sim[an_calls_avg:] - t_sim[:-an_calls_avg])
    # Now I express it in hours per microsecond
    el_time_ov_sim_time /= 3600.0/1e-6
else:
    el_time_ov_sim_time_computed = False
    print('quantity "elapsed time / simulated time" not computed, due to insufficient number of entries in times table')

# <codecell>
fig, ax = plt.subplots(nrows=2, ncols=2)

ax[0,0].plot(el_time_h, t_sim, '.-')
ax[0,0].set_ylabel('elapsed time / h')
ax[0,0].set_xlabel('simulated time / s')
ax[0,0].grid()

ax[1,0].plot(an_calls, dt_sim, '.-')
ax[1,0].set_xlabel('Analysis() calls')
ax[1,0].set_ylabel('dt / s')
ax[1,0].grid()

if el_time_ov_sim_time_computed:
    ax[0,1].plot(an_calls[an_calls_avg:], el_time_ov_sim_time, '.-')
    ax[0,1].set_ylabel('(el.time/sim.time) / (h/μs)')
    ax[0,1].set_xlabel('Analysis() calls')
    ax[0,1].grid()

plt.tight_layout()
plt.ioff()
plt.show()