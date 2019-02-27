#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2019, by the GROMACS development team, led by
# Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
# and including many others, as listed in the AUTHORS file in the
# top-level source directory and at http://www.gromacs.org.
#
# GROMACS is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation; either version 2.1
# of the License, or (at your option) any later version.
#
# GROMACS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with GROMACS; if not, see
# http://www.gnu.org/licenses, or write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
#
# If you want to redistribute modifications to GROMACS, please
# consider that scientific software is very special. Version
# control is crucial - bugs must be traceable. We will be happy to
# consider code for inclusion in the official distribution, but
# derived work must not be called official GROMACS. Details are found
# in the README & COPYING files - if they are missing, get the
# official version at http://www.gromacs.org.
#
# To help us fund GROMACS development, we humbly ask that you cite
# the research papers on the package. Check out http://www.gromacs.org.

# Python setuptools script to build and install the gmxapi Python interface
# from a GROMACS installation directory.

# Usage note: things go smoothly when we stick to the setup.py convention of
# having a package source directory with the same name as the package at the
# same level as the setup.py script and only expect `pip install .` in the
# setup.py directory. If we play with the layout more, it is hard to keep all
# of the `pip` and `setup.py` cases working as expected. This is annoying
# because running the Python interpreter immediately from the same directory
# can find the uninstalled source instead of the installed package. We can
# ease this pain by building an sdist in the enclosing CMake build scope
# and encouraging users to `pip install the_sdist.archive`. Otherwise, we
# just have to document that we only support full build-install of the Python
# package from the directory containing setup.py, which may clutter that
# directory with some artifacts.

from setuptools import setup

setup(
    name='gmxapi',

    # TODO: replace with CMake variables from GMXAPI version.
    version='0.1.0.dev0',
    python_requires='>=3.4, <4',
    setup_requires=['setuptools>=28'],

    packages=['gmxapi'],

    author='M. Eric Irrgang',
    author_email='info@gmxapi.org',
    description='gmxapi Python interface for GROMACS',
    license='LGPL',
    url='http://gmxapi.org/',

    # The installed package will contain compiled C++ extensions that cannot be loaded
    # directly from a zip file.
    zip_safe=False
)
