#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2016, by the GROMACS development team, led by
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

build_options = ['gcc-4.8', 'cmake-3.4.3']
extra_projects = [Project.REGRESSIONTESTS]

def run_build(context, cmake_opts):
    context.chdir(context.workspace.build_dir)
    context.run_cmake(cmake_opts)
    context.build_target(target=None)

    context.chdir(context.workspace.get_project_dir(Project.REGRESSIONTESTS))
    cmd = ['perl', 'gmxtest.pl', 'all']
    if cmake_opts['GMX_DOUBLE'] == 'ON':
        cmd += ['-double']
    context.run_cmd(cmd, failure_message='Regression tests failed to execute')

def do_build(context):
    cmake_opts=dict()
    cmake_opts['CMAKE_BUILD_TYPE'] = 'Reference'
    context.env.set_env_var('GMX_NO_TERM', '1')
    context.env.prepend_path_env(os.path.join(context.workspace.build_dir, 'bin'))

    cmake_opts['GMX_DOUBLE'] = 'ON'
    run_build(context, cmake_opts)

    cmake_opts['GMX_DOUBLE'] = 'OFF'
    run_build(context, cmake_opts)

    context.workspace.upload_revision(project=Project.REGRESSIONTESTS, file_glob='*reference*')
