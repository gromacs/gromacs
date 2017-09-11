#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2015,2016,2017, by the GROMACS development team, led by
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

# These are accessible later in the script, just like other
# declared options, via e.g. context.opts.release.
extra_options = {
    'mdrun-only': Option.simple,
    'static': Option.simple,
    'reference': Option.simple,
    'release': Option.simple,
    'release-with-assert': Option.simple,
    'release-with-debug-info': Option.simple,
    'asan': Option.simple,
    'mkl': Option.simple,
    'fftpack': Option.simple,
    'double': Option.simple,
    'thread-mpi': Option.bool,
    'gpu': Option.bool,
    'opencl': Option.bool,
    'openmp': Option.bool,
    'nranks': Option.string,
    'npme': Option.string,
    'gpu_id': Option.string
}

extra_projects = [Project.REGRESSIONTESTS]

def do_build(context):
    cmake_opts = dict()
    cmake_opts['GMX_COMPILER_WARNINGS'] = 'ON'
    cmake_opts['GMX_DEFAULT_SUFFIX'] = 'OFF'
    cmake_opts['CMAKE_BUILD_TYPE'] = 'Debug'
    cmake_opts['GMX_USE_RDTSCP'] = 'DETECT'

    if context.opts.reference:
        cmake_opts['CMAKE_BUILD_TYPE'] = 'Reference'
    elif context.opts['release']:
        cmake_opts['CMAKE_BUILD_TYPE'] = 'Release'
    elif context.opts['release-with-assert']:
        cmake_opts['CMAKE_BUILD_TYPE'] = 'RelWithAssert'
    elif context.opts['release-with-debug-info']:
        cmake_opts['CMAKE_BUILD_TYPE'] = 'RelWithDebInfo'
    elif context.opts.asan:
        cmake_opts['CMAKE_BUILD_TYPE'] = 'ASAN'
    elif context.opts.tsan:
        cmake_opts['CMAKE_BUILD_TYPE'] = 'TSAN'

    if context.opts.static:
        cmake_opts['BUILD_SHARED_LIBS'] = 'OFF'

    if context.opts.phi:
        cmake_opts['CMAKE_TOOLCHAIN_FILE'] = 'Platform/XeonPhi'

    if context.opts.double:
        cmake_opts['GMX_DOUBLE'] = 'ON'

    if context.opts.simd is None:
        cmake_opts['GMX_SIMD'] = 'None'
    else:
        cmake_opts['GMX_SIMD'] = context.opts.simd
    if context.opts.gpu or context.opts.opencl:
        cmake_opts['GMX_GPU'] = 'ON'
        if context.opts.opencl:
            context.env.set_env_var('CUDA_PATH', context.env.cuda_root)
            context.env.set_env_var('AMDAPPSDKROOT', context.env.amdappsdk_root)
            cmake_opts['GMX_USE_OPENCL'] = 'ON'
        else:
            cmake_opts['CUDA_TOOLKIT_ROOT_DIR'] = context.env.cuda_root
            cmake_opts['CUDA_HOST_COMPILER'] = context.env.cuda_host_compiler
    else:
        cmake_opts['GMX_GPU'] = 'OFF'
    if context.opts.thread_mpi is False:
        cmake_opts['GMX_THREAD_MPI'] = 'OFF'
    if context.opts.mpi:
        cmake_opts['GMX_MPI'] = 'ON'
    if context.opts.openmp is False:
        cmake_opts['GMX_OPENMP'] = 'OFF'

    if context.opts.mkl:
        cmake_opts['GMX_FFT_LIBRARY'] = 'mkl'
    elif context.opts.fftpack:
        cmake_opts['GMX_FFT_LIBRARY'] = 'fftpack'
    if context.opts.mkl or context.opts.atlas:
        cmake_opts['GMX_EXTERNAL_BLAS'] = 'ON'
        cmake_opts['GMX_EXTERNAL_LAPACK'] = 'ON'

    if context.opts.x11:
        cmake_opts['GMX_X11'] = 'ON'

    # At least hwloc on Jenkins produces a massive amount of reports about
    # memory leaks, which cannot be reasonably suppressed because ASAN cannot
    # produce a reasonable stack trace for them.
    if context.opts.asan:
        cmake_opts['GMX_HWLOC'] = 'OFF'

    regressiontests_path = context.workspace.get_project_dir(Project.REGRESSIONTESTS)

    if context.job_type == JobType.RELEASE:
        cmake_opts['REGRESSIONTEST_PATH'] = regressiontests_path
    else:
        if context.opts.mdrun_only:
            cmake_opts['GMX_BUILD_MDRUN_ONLY'] = 'ON'

    context.env.set_env_var('GMX_NO_TERM', '1')

    context.run_cmake(cmake_opts)
    context.build_target(target=None, keep_going=True)

    # TODO: Consider if it would be better to split this into a separate build
    # script, since it is somewhat different, even though it benefits from some
    # of the same build options.
    if context.job_type == JobType.RELEASE:
        context.build_target(target='check', keep_going=True)
        context.build_target(target='install')
        if context.opts.mdrun_only:
            context.workspace.clean_build_dir()
            cmake_opts['REGRESSIONTEST_PATH'] = None
            cmake_opts['GMX_BUILD_MDRUN_ONLY'] = 'ON'
            context.run_cmake(cmake_opts)
            context.build_target(target=None, keep_going=True)
            context.build_target(target='check', keep_going=True)
            context.build_target(target='install')
        gmxrc_cmd = '. ' + os.path.join(context.workspace.install_dir, 'bin', 'GMXRC')
        context.env.run_env_script(gmxrc_cmd)
        cmd = [os.path.join(regressiontests_path, 'gmxtest.pl'), '-nosuffix', 'all']
        if context.opts.mpi:
            cmd += ['-np', '1']
        if context.opts.double:
            cmd += ['-double']
        if context.opts.mdrun_only:
            cmd += ['-mdrun', 'mdrun']
        context.run_cmd(cmd, failure_message='Regression tests failed to execute')
        # TODO: Add testing for building the template.
        # TODO: Generalize the machinery here such that it can easily be used
        # also for non-release builds.
    else:
        context.build_target(target='tests', keep_going=True)

        context.run_ctest(args=['--output-on-failure'], memcheck=context.opts.asan)

        context.build_target(target='install')
        # TODO: Consider what could be tested about the installed binaries.

        if not context.opts.mdrun_only:
            context.env.prepend_path_env(os.path.join(context.workspace.build_dir, 'bin'))
            context.chdir(regressiontests_path)

            use_tmpi = not context.opts.mpi and context.opts.thread_mpi is not False

            cmd = 'perl gmxtest.pl -mpirun mpirun -xml -nosuffix all'

            # setting this stuff below is just a temporary solution,
            # it should all be passed as a proper the runconf from outside
            # The whole mechanism should be rethought in #1587.
            if context.opts.phi:
                cmd += ' -ntomp 28'
            elif context.opts.openmp:
                # OpenMP should always work when compiled in! Currently not set if
                # not explicitly set
                cmd += ' -ntomp 2'

            if context.opts.gpu_id:
                cmd += ' -gpu_id ' + context.opts.gpu_id

            if context.opts.nranks:
                nranks = context.opts.nranks
            else:
                nranks = '2'

            if context.opts.npme:
                cmd += ' -npme ' + context.opts.npme

            if context.opts.mpi:
                cmd += ' -np ' + nranks
            elif use_tmpi:
                cmd += ' -nt ' + nranks
            if context.opts.double:
                cmd += ' -double'
            if context.opts.asan:
                context.env.set_env_var('ASAN_OPTIONS', 'detect_leaks=0')
            context.run_cmd(cmd, shell=True, failure_message='Regression tests failed to execute')
