#!/usr/bin/env python
"""PISM time series generator"""

import os
import urllib
from math import pi
from netCDF4 import Dataset as NC
import numpy as np

data_sources = {
    'epica':    'ftp://ftp.ncdc.noaa.gov/pub/data/paleo/'
                'icecore/antarctica/epica_domec/edc3deuttemp2007.txt',
    'grip':     'ftp://ftp.ncdc.noaa.gov/pub/data/paleo/icecore/'
                'greenland/summit/grip/isotopes/gripd18o.txt',
    'odp1012':  'ftp://ftp.ncdc.noaa.gov/pub/data/paleo/'
                'contributions_by_author/herbert2001/odp1012.txt'}


def retrieve(rec):
    """Retrieve record data from the web unless local copy exists."""
    filename = rec + '.txt'
    if not os.path.isfile(filename):
        urllib.urlretrieve(data_sources[rec], filename)
    return filename


def extract(rec):
    """Extract temperature anomaly data from local file"""
    if rec == 'epica':
        time, temp = np.genfromtxt('epica.txt', delimiter=(4, 13, 17, 13, 13),
                                   skip_header=104, skip_footer=1,
                                   usecols = (2, 4), unpack=True)
    elif rec == 'grip':
        time, temp = np.genfromtxt('grip.txt', skip_header=37,
                                   usecols=(2, 1), unpack=True)
        temp = -11.88*(temp-temp[0]) - 0.1925*(temp**2-temp[0]**2)
    elif rec == 'odp1012':
        time, temp = np.genfromtxt('odp1012.txt', delimiter='\t',
                                   skip_header=1, usecols=(6, 8),
                                   missing_values=-999, usemask=True,
                                   unpack=True)
        temp.mask += time.mask
        time = time.compressed()*1000
        temp = temp.compressed()-temp[0]
    return time, temp


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
        retrieve(func)
        time, temp = extract(func)
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
                        choices=['ramp', 'cos', 'epica', 'grip', 'odp1012'])
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
