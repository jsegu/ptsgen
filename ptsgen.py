#!/usr/bin/env python
"""PISM time series generator"""

import os
import urllib
from math import pi
from netCDF4 import Dataset as NC
import numpy as np

_noaa = 'ftp://ftp.ncdc.noaa.gov/pub/data/paleo/'
data_sources = {
    'lapaz21p': _noaa + 'contributions_by_author/herbert2001/lapaz21p.txt',
    'odp1012':  _noaa + 'contributions_by_author/herbert2001/odp1012.txt',
    'odp1020':  _noaa + 'contributions_by_author/herbert2001/odp1020.txt',
    'domefuji': _noaa + 'icecore/antarctica/domefuji/df2012isotope-temperature.txt',
    'epica':    _noaa + 'icecore/antarctica/epica_domec/edc3deuttemp2007.txt',
    'vostok':   _noaa + 'icecore/antarctica/vostok/deutnat.txt',
    'ngrip':    _noaa + 'icecore/greenland/summit/ngrip/isotopes/ngrip-d18o-50yr.txt',
    'grip':     _noaa + 'icecore/greenland/summit/grip/isotopes/gripd18o.txt',
    'guliya':   _noaa + 'icecore/trop/guliya/guliya1997.txt'}


def retrieve(rec):
    """Retrieve record data from the web unless local copy exists."""
    filename = rec + '.txt'
    if not os.path.isfile(filename):
        urllib.urlretrieve(data_sources[rec], filename)
    return filename


def extract(rec):
    """Extract temperature anomaly data from local file"""
    txtkw = {
        'domefuji': {'skip_header': 1795, 'usecols': (0, 4)},
        'epica':    {'delimiter': (4, 13, 17, 13, 13),
                     'skip_header': 104, 'skip_footer': 1, 'usecols': (2, 4)},
        'vostok':   {'skip_header': 111, 'usecols': (1, 3)},
        'ngrip':    {'skip_header': 80, 'usecols': (0, 1),
                     'converters': {0: lambda s: float(s.replace(',',''))}},
        'grip':     {'skip_header': 37, 'usecols': (2, 1)},
        'guliya':   {'skip_header': 445, 'skip_footer': 33, 'usecols': (0, 1)},
        'lapaz21p': {'delimiter': '\t', 'skip_header': 1, 'usecols': (2, 5),
                     'missing_values': -999, 'usemask': True},
        'odp1012':  {'delimiter': '\t', 'skip_header': 1, 'usecols': (6, 8),
                     'missing_values': -999, 'usemask': True},
        'odp1020':  {'delimiter': '\t', 'skip_header': 1, 'usecols': (4, 7),
                     'missing_values': -999, 'usemask': True}}[rec]
    time, temp = np.genfromtxt(rec + '.txt', unpack=True, **txtkw)
    if rec == 'ngrip':
        temp = temp[::2]
        time = time[::2]
    if rec == 'domefuji':
        time *= 1000
    elif rec in ('grip', 'ngrip'):
        temp = -11.88*(temp-temp[0]) - 0.1925*(temp**2-temp[0]**2)
    elif rec == 'guliya':
        time *= 1000
        temp -= temp[0]
    elif rec in ('lapaz21p', 'odp1012', 'odp1020'):
        time.mask += temp.mask
        temp.mask += time.mask
        time = time.compressed()*1000
        temp = temp.compressed()-temp[0]
    return time, temp


def generate(func, tmin, tmax, orig, ampl, n=101, output=None,
             scale_interval=(-32e3, -22e3), var='delta_T', unit='K',
             smoothing=None):
    """Generate NetCDF file"""

    # initialize netCDF file
    nc = NC(output or 'delta_%s.nc' % var, 'w', format='NETCDF3_CLASSIC')
    nc.createDimension('time', 0)
    timevar = nc.createVariable('time', 'f4', 'time')
    timevar.units = 'years'
    tempvar = nc.createVariable(var, 'f4', 'time')
    tempvar.units = unit

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
        t1, t2 = scale_interval
        temp /= temp[(t1 < time) * (time < t2)].mean()
        if time[-1] < tmax:
            time = np.append(time, tmax)
            temp = np.append(temp, temp[-1])
        if time[0] > tmin:
            time = np.insert(time, 0, tmin)
            temp = np.insert(temp, 0, temp[0])
        timevar[:] = time
        tempvar[:] = orig + ampl*temp

    # optional smoothing
    if smoothing:
        window = np.ones(smoothing)/smoothing
        tempvar[:] = np.convolve(tempvar[:], window, mode='same')

    # close file
    nc.close()


if __name__ == "__main__":
    import argparse

    # Argument parser
    parser = argparse.ArgumentParser(
        description='''Scalar time series generator for PISM''')
    parser.add_argument('func', type=str, help='Function to use',
                        choices=['ramp', 'cos',
                                 'lapaz21p', 'odp1012', 'odp1020',
                                 'domefuji', 'epica', 'vostok',
                                 'ngrip', 'grip', 'guliya'])
    parser.add_argument('tmin', type=float, help='Start time in years')
    parser.add_argument('tmax', type=float, help='End time in years')
    parser.add_argument('orig', type=float, help='Value at origin')
    parser.add_argument('ampl', type=float, help='Amplitude')
    parser.add_argument('-n', '--length', type=int, default=101,
                        help='Length of time series (default: %(default)s)')
    parser.add_argument('-o', '--output', help='output file')
    parser.add_argument('-i', '--scale-interval', type=float, nargs=2,
                        default=(-32e3, -22e3), metavar=('T1', 'T2'),
                        help='Record scaling time interval (default: %(default)s)')
    parser.add_argument('-v', '--variable', default='delta_T',
                        help='Variable name (default: %(default)s)')
    parser.add_argument('-u', '--unit', default='K',
                        help='Variable unit (default: %(default)s)')
    parser.add_argument('-s', '--smoothing', type=int, default=None,
                        help='Optional smoothing window lenght (default: %(default)s)')
    args = parser.parse_args()

    # Generate netCDF file
    generate(args.func, args.tmin, args.tmax, args.orig, args.ampl,
             args.length, args.output, args.scale_interval, args.variable,
             args.unit, args.smoothing)
