#!/usr/bin/env python
# Copyright (c) 2014--2018, Julien Seguinot <seguinot@vaw.baug.ethz.ch>
# GNU General Public License v3.0+ (https://www.gnu.org/licenses/gpl-3.0.txt)

"""Plot test examples"""

from netCDF4 import Dataset
from matplotlib import pyplot as plt

plt.rc('font', size=8)
funclist = ['cos', 'ramp',
            'lapaz21p', 'odp1012', 'odp1020',
            'domefuji', 'epica',  'vostok',
            'ngrip', 'grip', 'guliya', 'md012444']
exptlist = range(4)

# initialize figure
fig, grid = plt.subplots(1, 2)

# plot various functions
ax = grid[0]
for i, func in enumerate(funclist):
    nc = Dataset(func + '.nc')
    time = nc.variables['time'][:]/1000.
    temp = nc.variables['delta_T'][:]
    nc.close()
    ax.plot(time, temp-5*i, label=func)
    ax.set_xlim((-120, 0))
    ax.set_xlabel('time (ka)')
    ax.set_ylabel('delta_T (K)')
    ax.legend(loc='upper left')

# plot smoothing effect
ax = grid[1]
for i, expt in enumerate(exptlist):
    nc = Dataset('grip-1e%d.nc' % expt)
    time = nc.variables['time'][:]/1000.
    temp = nc.variables['delta_T'][:]
    nc.close()
    ax.plot(time, temp-5*i, label='smooth %d' % (10**expt))
    ax.set_xlim((-120, 0))
    ax.set_xlabel('time (ka)')
    ax.set_ylabel('delta_T (K)')
    ax.legend(loc='upper left')

# add labels and save
plt.show()
