#!/usr/bin/env python
# Copyright (c) 2018, Julien Seguinot <seguinot@vaw.baug.ethz.ch>
# GNU General Public License v3.0+ (https://www.gnu.org/licenses/gpl-3.0.txt)

"""
Prepare scalar modifiers for the Parallel Ice Sheet Model.
"""

import sys
from time import strftime
from math import pi
from netCDF4 import Dataset as NC
import numpy as np


_NOAA = 'https://www1.ncdc.noaa.gov/pub/data/paleo/'
_NOAA_ANTARCTICA = _NOAA + 'icecore/antarctica/'
_NOAA_GREENLAND = _NOAA + 'icecore/greenland/'
RECORDS = dict(
    lapaz21p=dict(
        fname=_NOAA+'contributions_by_author/herbert2001/lapaz21p.txt',
        delimiter='\t', missing_values=-999, skip_header=1, usecols=(2, 5),
        usemask=True),
    odp1012=dict(
        fname=_NOAA+'contributions_by_author/herbert2001/odp1012.txt',
        delimiter='\t', missing_values=-999, skip_header=1, usecols=(6, 8),
        usemask=True),
    odp1020=dict(
        fname=_NOAA+'contributions_by_author/herbert2001/odp1020.txt',
        delimiter='\t', missing_values=-999, skip_header=1, usecols=(4, 7),
        usemask=True),
    domefuji=dict(
        fname=_NOAA_ANTARCTICA+'domefuji/df2012isotope-temperature.txt',
        skip_header=1795, usecols=(0, 4)),
    epica=dict(
        fname=_NOAA_ANTARCTICA+'epica_domec/edc3deuttemp2007.txt',
        delimiter=(4, 13, 17, 13, 13), skip_header=104, skip_footer=1,
        usecols=(2, 4)),
    vostok=dict(
        fname=_NOAA_ANTARCTICA+'vostok/deutnat.txt',
        skip_header=111, usecols=(1, 3)),
    ngrip=dict(
        fname=_NOAA_GREENLAND+'summit/ngrip/isotopes/ngrip-d18o-50yr.txt',
        converters={0: lambda s: float(s.replace(',', ''))},
        skip_header=80, usecols=(0, 1)),
    grip=dict(
        fname=_NOAA_GREENLAND+'summit/grip/isotopes/gripd18o.txt',
        skip_header=37, usecols=(2, 1)),
    guliya=dict(
        fname=_NOAA+'icecore/trop/guliya/guliya1997.txt',
        skip_header=445, skip_footer=33, usecols=(0, 1)),
    md012444=dict(
        fname='https://doi.pangaea.de/10.1594/PANGAEA.771891?format=textfile',
        delimiter='\t', skip_header=15, usecols=(1, 2)))


def extract(rec):
    """Extract anomaly data from local file"""
    time, data = np.genfromtxt(encoding='latin1', unpack=True, **RECORDS[rec])
    if rec == 'ngrip':
        data = data[::2]
        time = time[::2]
    if rec == 'domefuji':
        time *= 1000
    elif rec in ('grip', 'ngrip'):
        data = -11.88*(data-data[0]) - 0.1925*(data**2-data[0]**2)
    elif rec in ('guliya', 'md012444'):
        time *= 1000
        data -= data[0]
    elif rec in ('lapaz21p', 'odp1012', 'odp1020'):
        time.mask += data.mask
        data.mask += time.mask
        time = time.compressed()*1000
        data = data.compressed()-data[0]
    if rec == 'md012444':
        argsort = time.argsort()
        time = time[argsort]
        data = data[argsort]
    return time, data


def generate(func, tmin, tmax, orig, ampl, n=101, output=None,
             scale_interval=(-32e3, -22e3), var='delta_T', unit='K',
             smoothing=None):
    """Generate NetCDF file"""

    # initialize netCDF file
    nc = NC(output or 'delta_%s.nc' % var, 'w', format='NETCDF3_CLASSIC')
    nc.createDimension('time', 0)
    timevar = nc.createVariable('time', 'f4', 'time')
    timevar.units = 'years'
    datavar = nc.createVariable(var, 'f4', 'time')
    datavar.units = unit

    # in case of a regular function
    if func in ('ramp', 'cos'):
        t = np.linspace(0, 1, n)
        time = tmin + (tmax-tmin)*t
        if func == 'ramp':
            data = t
        elif func == 'cos':
            data = (1-np.cos(2*pi*t))/2

    # in case of a proxy record
    else:
        time, data = extract(func)
        time = -time[::-1]
        data = data[::-1]
        t1, t2 = scale_interval
        data /= data[(t1 < time) * (time < t2)].mean()
        if time[-1] < tmax:
            time = np.append(time, tmax)
            data = np.append(data, data[-1])
        if time[0] > tmin:
            time = np.insert(time, 0, tmin)
            data = np.insert(data, 0, data[0])

    # optional smoothing
    if smoothing:
        window = np.ones(smoothing)/smoothing
        data = np.convolve(data, window, mode='same')

    # linear scaling
    data = orig + ampl*data

    # assign values to variables
    timevar[:] = time
    datavar[:] = data

    # add history line
    nc.history = strftime('%Y-%m-%d %H:%M:%S %Z: ') + ' '.join(sys.argv)

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
                                 'ngrip', 'grip', 'guliya',
                                 'md012444'])
    parser.add_argument('tmin', type=float, help='Start time in years')
    parser.add_argument('tmax', type=float, help='End time in years')
    parser.add_argument('orig', type=float, help='Value at origin')
    parser.add_argument('ampl', type=float, help='Amplitude')
    parser.add_argument('-n', '--length', type=int, default=101,
                        help='Length of time series (default: %(default)s)')
    parser.add_argument('-o', '--output', help='output file')
    parser.add_argument('-i', '--scale-interval', type=float, nargs=2,
                        default=(-32e3, -22e3), metavar=('T1', 'T2'),
                        help='Record scaling interval (default: %(default)s)')
    parser.add_argument('-v', '--variable', default='delta_T',
                        help='Variable name (default: %(default)s)')
    parser.add_argument('-u', '--unit', default='K',
                        help='Variable unit (default: %(default)s)')
    parser.add_argument('-s', '--smoothing', type=int, default=None,
                        help='Smoothing window lenght (default: %(default)s)')
    args = parser.parse_args()

    # Generate netCDF file
    generate(args.func, args.tmin, args.tmax, args.orig, args.ampl,
             args.length, args.output, args.scale_interval, args.variable,
             args.unit, args.smoothing)
