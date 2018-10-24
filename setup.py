#!/usr/bin/env python
# Copyright (c) 2018, Julien Seguinot <seguinot@vaw.baug.ethz.ch>
# GNU General Public License v3.0+ (https://www.gnu.org/licenses/gpl-3.0.txt)

"""Build script for pism-palseries."""

import setuptools

with open('README.rst', 'r') as f:
    README = f.read()

setuptools.setup(
    name='pism-palseries',
    version='0.1.0',
    author='Julien Seguinot',
    author_email='seguinot@vaw.baug.ethz.ch',
    description='Prepare scalar modifiers for the Parallel Ice Sheet Model',
    long_description=README,
    long_description_content_type='text/x-rst',
    url='http://github.com/juseg/pism-palseries',
    license='gpl-3.0',
    install_requires=['netCDF4', 'numpy'],
    py_modules=['pism_palseries'],
    scripts=['pism_palseries.py'],
)
