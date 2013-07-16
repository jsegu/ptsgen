#!/usr/bin/env python

from math import pi
from netCDF4 import Dataset as NC
import numpy as np

def generate(func, tmin, tmax, orig, ampl, n=101, output=None, var='T'):
    """Generate NetCDF file"""

    # initialize netCDF file
    nc = NC(output or 'delta_%s.nc' % var, 'w', format='NETCDF3_CLASSIC')
    timedim = nc.createDimension('time', n)
    timevar = nc.createVariable('time', 'f4', 'time')
    timevar.units = 'years'
    tempvar = nc.createVariable('delta_%s' % var, 'f4', 'time')
    tempvar.units = {'T': 'K', 'P': 'm year-1'}[var]

    # write data and close
    t = np.linspace(0, 1, n)
    timevar[:] = tmin + (tmax-tmin)*t
    if func == 'ramp':
      tempvar[:] = orig + ampl*t
    elif func == 'cos':
      tempvar[:] = orig + ampl/2*(1-np.cos(2*pi*t))
    nc.close()

if __name__ == "__main__":
    import argparse

    # Argument parser
    parser = argparse.ArgumentParser(
      description='''Scalar offsets time series generator for PISM''')
    parser.add_argument('func', type=str, help='Function to use',
      choices=['ramp', 'cos'])
    parser.add_argument('tmin', type=float, help='Start time in years')
    parser.add_argument('tmax', type=float, help='End time in years')
    parser.add_argument('orig', type=float, help='Value at origin')
    parser.add_argument('ampl', type=float, help='Amplitude')
    parser.add_argument('-n', '--length', type=int, default=101,
      help='Length of time series (default: %(default)s)')
    parser.add_argument('-o', '--output', help='output file')
    parser.add_argument('-v', '--variable',
      choices=['T', 'P'], default='T',
      help='Variable name (default: %(default)s)')
    args = parser.parse_args()

    # Generate netCDF file
    generate(args.func, args.tmin, args.tmax, args.orig, args.ampl,
      args.length, args.output, args.variable)
