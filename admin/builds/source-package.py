#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2015, by the GROMACS development team, led by
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

import os
import re

build_out_of_source = True

def do_build(context):
    cmake_opts = {
            'GMX_BUILD_HELP': 'ON',
            # TODO: Consider encapsulating this within releng.
            # The environment variable is intended to be set as a build
            # parameter (or otherwise in the Jenkins job configuration).
            'GMX_BUILD_TARBALL': os.getenv('RELEASE', None),
            'CMAKE_BUILD_TYPE': 'Debug',
            'GMX_SIMD': 'None',
            'GMX_THREAD_MPI': 'OFF',
            'GMX_OPENMP': 'OFF',
            'GMX_GPU': 'OFF'
        }

    context.run_cmake(cmake_opts)
    context.build_target(target='gmx')
    context.build_target(target='man')
    context.build_target(target='completion')
    context.build_target(target='install-guide')

    context.build_target(target='package_source')

    cpack_config = read_source_package_config(context.workspace.build_dir)
    package_name = cpack_config['CPACK_PACKAGE_FILE_NAME'] + '.tar.gz'
    package_info = {
            'SOURCE_PACKAGE_FILE_NAME': package_name,
            'SOURCE_PACKAGE_VERSION': cpack_config['CPACK_PACKAGE_VERSION'],
            'SOURCE_MD5SUM': context.compute_md5(package_name)
        }
    log_path = context.workspace.get_path_for_logfile('package-info.log')
    context.write_property_file(log_path, package_info)

def read_source_package_config(build_dir):
    values = dict()
    set_re = r'SET\((\w+)\s*"(.*)"\)\s*'
    with open(os.path.join(build_dir, 'CPackSourceConfig.cmake'), 'r') as fp:
        for line in fp:
            match = re.match(set_re, line)
            if match:
                values[match.group(1)] = match.group(2)
    return values
