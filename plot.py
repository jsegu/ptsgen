#!/usr/bin/env python2

from netCDF4 import Dataset
from matplotlib import pyplot as plt

# plot
for func in ['cos', 'ramp', 'grip', 'epica', 'odp1012']:
    nc = Dataset('%s.nc' % func)
    time = nc.variables['time'][:]/1000.
    temp = nc.variables['delta_T'][:]
    nc.close()
    plt.plot(time, temp)

# add labels and save
plt.xlim((-120, 0))
plt.xlabel('time (ka)')
plt.ylabel('delta_T (K)')
plt.show()
