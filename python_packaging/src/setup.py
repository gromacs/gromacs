#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2019,2020, by the GROMACS development team, led by
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

import os

# Allow setup.py to be run when scikit-build is not installed, such as to
# produce source distribution archives with `python setup.py sdist`
try:
    from skbuild import setup
except ImportError:
    from distutils.core import setup

usage = """
The `gmxapi` package requires an existing GROMACS installation, version 2020 or higher.
To specify the GROMACS installation to use, provide a GMXTOOLCHAINDIR
environment variable when running setup.py or `pip`.

Example:
    GMXTOOLCHAINDIR=/path/to/gromacs/share/cmake/gromacs pip install gmxapi

If you have multiple builds of GROMACS, distinguished by a suffix `$SUFFIX`, the
tool chain directory will use that suffix.

Example:
    GMXTOOLCHAINDIR=/path/to/gromacs/share/cmake/gromacs$SUFFIX pip install gmxapi

In the example, `gmxapi` is downloaded automatically from pypi.org. You can
replace `gmxapi` with a local directory or archive file to build from a source
distribution.

setup.py will use the location of GMXTOOLCHAINDIR to locate the
gmxapi library configured during GROMACS installation. Alternatively, if
gmxapi_DIR is provided, or if GMXRC has been "sourced", the toolchain file
location may be deduced. Note, though, that if multiple GROMACS installations
exist in the same location (with different suffixes) only the first one will be
used when guessing a toolchain, because setup.py does not know which corresponds
to the gmxapi support library.

If specifying GMXTOOLCHAINDIR and gmxapi_DIR, the tool chain directory must be 
located within a subdirectory of gmxapi_DIR.

Refer to project web site for complete documentation.

"""


class GmxapiInstallError(Exception):
    """Error processing setup.py for gmxapi Python package."""


gmx_toolchain_dir = os.getenv('GMXTOOLCHAINDIR')
gmxapi_DIR = os.getenv('gmxapi_DIR')
if gmxapi_DIR is None:
    # Infer from GMXRC exports, if available.
    gmxapi_DIR = os.getenv('GROMACS_DIR')

def _find_first_gromacs_suffix(directory):
    dir_contents = os.listdir(directory)
    for entry in dir_contents:
        if entry.startswith('gromacs'):
            return entry.strip('gromacs')

if gmx_toolchain_dir is None:
    # Try to guess from standard GMXRC environment variables.
    if gmxapi_DIR is not None:
        if os.path.exists(gmxapi_DIR) and os.path.isdir(gmxapi_DIR):
            share_cmake = os.path.join(gmxapi_DIR, 'share', 'cmake')
            suffix = _find_first_gromacs_suffix(share_cmake)
            if suffix is not None:
                gmx_toolchain_dir = os.path.join(share_cmake, 'gromacs' + suffix)

if gmx_toolchain_dir is None:
    print(usage)
    raise GmxapiInstallError('Could not configure for GROMACS installation. Provide GMXTOOLCHAINDIR.')

suffix = os.path.basename(gmx_toolchain_dir).strip('gromacs')
gmx_toolchain = os.path.abspath(os.path.join(gmx_toolchain_dir, 'gromacs-toolchain' + suffix + '.cmake'))

if gmxapi_DIR is None:
    # Example: given /usr/local/gromacs/share/cmake/gromacs/gromacs-toolchain.cmake,
    # we would want /usr/local/gromacs.
    # Note that we could point more directly to the gmxapi-config.cmake but,
    # so far, we have relied on CMake automatically looking into
    # <package>_DIR/share/cmake/<package>/ for such a file.
    # We would need a slightly different behavior for packages that link against
    # libgromacs directly, as sample_restraint currently does.
    gmxapi_DIR = os.path.join(os.path.dirname(gmx_toolchain), '..', '..', '..')

gmxapi_DIR = os.path.abspath(gmxapi_DIR)

if not os.path.exists(gmxapi_DIR) or not os.path.isdir(gmxapi_DIR):
    print(usage)
    raise GmxapiInstallError('Please set a valid gmxapi_DIR.')

if gmxapi_DIR != os.path.commonpath([gmxapi_DIR, gmx_toolchain]):
    raise GmxapiInstallError('GROMACS toolchain file {} is not in gmxapi_DIR {}'.format(
        gmx_toolchain,
        gmxapi_DIR
    ))

cmake_platform_hints = '-DCMAKE_TOOLCHAIN_FILE={}'.format(gmx_toolchain)
# Note that <package>_ROOT is not standard until CMake 3.12
# Reference https://cmake.org/cmake/help/latest/policy/CMP0074.html#policy:CMP0074
cmake_gmxapi_hint = '-Dgmxapi_ROOT={}'.format(gmxapi_DIR)
cmake_args = [cmake_platform_hints, cmake_gmxapi_hint]

setup(
    name='gmxapi',

    # TODO: single-source version information (currently repeated in gmxapi/version.py)
    version='0.2.0b1',
    python_requires='>=3.6, <3.9',
    install_requires=['networkx>=2.0',
                      'numpy>=1'],

    packages=['gmxapi', 'gmxapi.simulation'],

    cmake_args=cmake_args,

    author='M. Eric Irrgang',
    author_email='info@gmxapi.org',
    description='gmxapi Python interface for GROMACS',
    license='LGPL',
    url='http://gmxapi.org/',

    # The installed package will contain compiled C++ extensions that cannot be loaded
    # directly from a zip file.
    zip_safe=False
)
