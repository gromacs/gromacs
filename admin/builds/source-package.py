#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2015,2016, by the GROMACS development team, led by
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

import os.path

build_out_of_source = True

def do_build(context):
    cmake_opts = {
            'GMX_BUILD_HELP': 'ON',
            'CMAKE_BUILD_TYPE': 'Debug',
            'GMX_SIMD': 'None',
            'GMX_THREAD_MPI': 'OFF',
            'GMX_OPENMP': 'OFF',
            'GMX_GPU': 'OFF'
        }
    if context.params.get('RELEASE', Parameter.bool):
        cmake_opts['GMX_BUILD_TARBALL'] = 'ON'

    context.run_cmake(cmake_opts)
    context.build_target(target='gmx')
    context.build_target(target='man')
    context.build_target(target='completion')
    context.build_target(target='install-guide')

    context.build_target(target='package_source')

    cpack_config_path = os.path.join(context.workspace.build_dir, 'CPackSourceConfig.cmake')
    cpack_config = context.read_cmake_variable_file(cpack_config_path)
    package_name = cpack_config['CPACK_PACKAGE_FILE_NAME'] + '.tar.gz'
    version = cpack_config['CPACK_PACKAGE_VERSION']
    context.write_package_info(Project.GROMACS, package_name, version)
