#!/usr/bin/env python2
"""Plot test examples"""

from netCDF4 import Dataset
from matplotlib import pyplot as plt

plt.rc('font', size=8)
funclist = ['cos', 'ramp',
            'lapaz21p', 'odp1012', 'odp1020',
            'domefuji', 'epica',  'vostok',
            'ngrip', 'grip', 'guliya']

# plot
for i, func in enumerate(funclist):
    nc = Dataset(func + '.nc')
    time = nc.variables['time'][:]/1000.
    temp = nc.variables['delta_T'][:]
    nc.close()
    plt.plot(time, temp-5*i, label=func)

# add labels and save
plt.xlim((-120, 0))
plt.xlabel('time (ka)')
plt.ylabel('delta_T (K)')
plt.legend(loc='upper left')
plt.show()
