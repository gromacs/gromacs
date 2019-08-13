#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2015,2016,2017,2018,2019, by the GROMACS development team, led by
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

# These options need to match the node label in Jenkins and the
# capabilities in releng/agents.py for the agent where the analysis is
# intended to run.
build_options = ['clang-8', 'clang-static-analyzer-8']

# Policy global variables
use_stdlib_through_env_vars = False

def do_build(context):
    cmake_opts = {
            'CMAKE_BUILD_TYPE': 'Debug',
            # Build tests as part of the all target.
            'GMX_DEVELOPER_BUILD': 'ON',
            'GMX_GPU': 'OFF',
            'GMX_OPENMP': 'OFF',
            'GMX_SIMD': 'None',
            'GMX_USE_RDTSCP': 'OFF',
            'GMX_FFT_LIBRARY': 'fftpack',
            'GMX_CLANG_ANALYZER' : 'ON'
        }
    context.env.append_to_env_var('CXXFLAGS', '-stdlib=libc++')

    context.run_cmake(cmake_opts)
    context.build_target(target=None)
    context.process_clang_analyzer_results()
