#!/usr/bin/env python

import os
import urllib
from math import pi
from netCDF4 import Dataset as NC
import numpy as np


def generate(func, tmin, tmax, orig, ampl, n=101, output=None, var='T'):
    """Generate NetCDF file"""

    # initialize netCDF file
    nc = NC(output or 'delta_%s.nc' % var, 'w', format='NETCDF3_CLASSIC')
    nc.createDimension('time', 0)
    timevar = nc.createVariable('time', 'f4', 'time')
    timevar.units = 'years'
    tempvar = nc.createVariable('delta_%s' % var, 'f4', 'time')
    tempvar.units = {'T': 'K', 'P': 'm year-1'}[var]

    # in case of a regular function
    if func in ('ramp', 'cos'):
        t = np.linspace(0, 1, n)
        timevar[:] = tmin + (tmax-tmin)*t
        if func == 'ramp':
            tempvar[:] = orig + ampl*t
        elif func == 'cos':
            tempvar[:] = orig + ampl/2*(1-np.cos(2*pi*t))

    # in case of a proxy record
    else:
        if func == 'epica':
            if not os.path.isfile('epica.txt'):
                url = ('ftp://ftp.ncdc.noaa.gov/pub/data/paleo/'
                       'icecore/antarctica/epica_domec/edc3deuttemp2007.txt')
                urllib.urlretrieve(url, 'epica.txt')
            time, temp = np.genfromtxt('epica.txt', dtype='f4',
                                       delimiter=(4, 13, 17, 13, 13),
                                       skip_header=104, skip_footer=1,
                                       usecols = (2, 4), unpack=True)
        elif func == 'grip':
            if not os.path.isfile('grip.txt'):
                url = ('ftp://ftp.ncdc.noaa.gov/pub/data/paleo/icecore/'
                       'greenland/summit/grip/isotopes/gripd18o.txt')
                urllib.urlretrieve(url, 'grip.txt')
            time, d18o = np.genfromtxt('grip.txt', skip_header=37,
                                       usecols=(2, 1), unpack=True)
            temp = -11.88*(d18o-d18o[0]) - 0.1925*(d18o**2-d18o[0]**2)
        time = -time[::-1]
        temp = temp[::-1]
        avep = (-30e3 < time) * (time < -20e3)
        tmin = temp[avep].mean()
        temp = temp/tmin
        timevar[:] = time
        tempvar[:] = orig + ampl*temp

    # close file
    nc.close()


if __name__ == "__main__":
    import argparse

    # Argument parser
    parser = argparse.ArgumentParser(
        description='''Scalar offsets time series generator for PISM''')
    parser.add_argument('func', type=str, help='Function to use',
                        choices=['ramp', 'cos', 'epica', 'grip'])
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
