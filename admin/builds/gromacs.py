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

extra_projects = [Project.REGRESSIONTESTS]

def do_build(context):
    cmake_opts=dict()
    cmake_opts['GMX_DEFAULT_SUFFIX'] = 'OFF'
    cmake_opts['CMAKE_BUILD_TYPE'] = 'Debug'

    if context.params.build_type is None:
        pass
    elif context.params.build_type == BuildType.REFERENCE:
        cmake_opts['CMAKE_BUILD_TYPE'] = 'Reference'
    elif context.params.build_type == BuildType.OPTIMIZED:
        cmake_opts['CMAKE_BUILD_TYPE'] = 'RelWithAssert'
    elif context.params.build_type == BuildType.PERFORMANCE:
        cmake_opts['CMAKE_BUILD_TYPE'] = 'Release'
    elif context.params.build_type == BuildType.ASAN:
        cmake_opts['CMAKE_BUILD_TYPE'] = 'ASAN'
    elif context.params.build_type == BuildType.TSAN:
        cmake_opts['CMAKE_BUILD_TYPE'] = 'TSAN'

    if context.params.phi:
        cmake_opts['CMAKE_TOOLCHAIN_FILE'] = 'Platform/XeonPhi'

    if context.params.double:
        cmake_opts['GMX_DOUBLE'] = 'ON'

    if context.params.simd is None:
        cmake_opts['GMX_SIMD'] = 'None'
    else:
        cmake_opts['GMX_SIMD'] = context.params.simd
    if context.params.gpu:
        cmake_opts['GMX_GPU'] = 'ON'
        cmake_opts.update(context.get_cuda_cmake_options())
    else:
        cmake_opts['GMX_GPU'] = 'OFF'
    if context.params.thread_mpi is False:
        cmake_opts['GMX_THREAD_MPI'] = 'OFF'
    if context.params.mpi:
        cmake_opts['GMX_MPI'] = 'ON'
    if context.params.openmp is False:
        cmake_opts['GMX_OPENMP'] = 'OFF'

    if context.params.fft_library is None:
        pass
    elif context.params.fft_library == FftLibrary.MKL:
        cmake_opts['GMX_FFT_LIBRARY'] = 'mkl'
    elif context.params.fft_library == FftLibrary.FFTPACK:
        cmake_opts['GMX_FFT_LIBRARY'] = 'fftpack'
    if context.params.external_linalg:
        cmake_opts['GMX_EXTERNAL_BLAS'] = 'ON'
        cmake_opts['GMX_EXTERNAL_LAPACK'] = 'ON'

    if context.params.x11:
        cmake_opts['GMX_X11'] = 'ON'

    if context.params.mdrun_only:
        cmake_opts['GMX_BUILD_MDRUN_ONLY'] = 'ON'

    context.env.add_env_var('GMX_NO_TERM', '1')

    context.run_cmake(cmake_opts)
    context.build_target(target=None, keep_going=True)
    context.build_target(target='tests', keep_going=True)

    if not context.params.memcheck:
        context.run_ctest(args=['--output-on-failure'])
    else:
        context.run_ctest(args=['-LE', 'GTest', '--output-on-failure'])
        context.run_ctest(args=['-L', 'GTest', '--output-on-failure'], memcheck=True)

    if not context.params.mdrun_only:
        context.env.prepend_path_env(os.path.join(context.workspace.build_dir, 'bin'))
        os.chdir(context.workspace.get_project_dir(Project.REGRESSIONTESTS))

        if not context.params.mpi and context.params.thread_mpi is not False:
            use_tmpi = True

        cmd = 'perl gmxtest.pl -mpirun mpirun -xml -nosuffix all'
        if context.params.build_type == BuildType.ASAN:
            cmd+=' -parse asan_symbolize.py'

        # setting this stuff below is just a temporary solution,
        # it should all be passed as a proper the runconf from outside
        # The whole mechanism should be rethought in #1587.
        if context.params.phi:
            cmd += ' -ntomp 28'
        elif context.params.openmp:
            # OpenMP should always work when compiled in! Currently not set if
            # not explicitly set
            cmd += ' -ntomp 2'

        if context.params.gpu:
            if context.params.mpi or use_tmpi:
                gpu_id = '01' # for (T)MPI use the two GT 640-s
            else:
                gpu_id = '0' # use GPU #0 by default
            cmd += ' -gpu_id ' + gpu_id

        if any([x in ('-np', '-nt') for x in context.params.extra_gmxtest_args]):
            # Use values passed
            pass
        elif context.params.mpi:
            cmd += ' -np 2'
        elif use_tmpi:
            cmd += ' -nt 2'
        if context.params.double:
            cmd += ' -double'
        cmd += ' ' + ' '.join(context.params.extra_gmxtest_args)
        context.run_cmd_with_env(cmd, shell=True, failure_message='Regression tests failed to execute')
