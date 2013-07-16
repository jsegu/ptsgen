#!/usr/bin/env python

from math import pi
from netCDF4 import Dataset as NC
import numpy as np

def generate(output, var, tmin, tmax, n, amp):
    """Generate NetCDF file"""

    # initialize netCDF file
    nc = NC(output, 'w', format='NETCDF3_CLASSIC')
    timedim = nc.createDimension('time', n)
    timevar = nc.createVariable('time', 'f4', 'time')
    timevar.units = 'years'
    tempvar = nc.createVariable('delta_%s' % var, 'f4', 'time')
    tempvar.units = {'T': 'K', 'P': 'm year-1'}[var]

    # write data and close
    timevar[:] = np.linspace(tmin, tmax, n)
    tempvar[:] = (np.cos(2*pi*timevar[:]/(tmax-tmin)) - 1) * amp/2
    nc.close()

if __name__ == "__main__":
    import argparse

    # Argument parser
    parser = argparse.ArgumentParser(
      description='''Scalar offsets time series generator for PISM''')
    parser.add_argument('tmin', type=float, help='Start time in years')
    parser.add_argument('tmax', type=float, help='End time in years')
    parser.add_argument('dt',   type=float, help='Time interval in years')
    parser.add_argument('amp',  type=float, help='Amplitude in variable unit')
    parser.add_argument('-o', '--output', help='output file')
    parser.add_argument('-v', '--variable',
      choices=['T', 'P'], default='T',
      help='Variable name (default: %(default)s)')
    args = parser.parse_args()

    # Parse arguments
    tmin = args.tmin
    tmax = args.tmax
    dt   = args.dt
    amp  = args.amp
    var  = args.variable

    # Output file name
    output = args.output or 'cos%gk-cool%02g.nc' % ((tmax-tmin)/1000, amp)

    # Generate netCDF file
    generate(output, var, tmin, tmax, dt, amp)
