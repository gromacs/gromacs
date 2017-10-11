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

import os
import re

build_out_of_source = True

extra_options = {
    'source-md5': Option.string
}

def do_build(context):
    cmake_opts = {
            'GMX_BUILD_HELP': 'ON',
            'GMX_BUILD_MANUAL': 'ON',
            'SOURCE_MD5SUM': context.opts.source_md5,
            'CMAKE_BUILD_TYPE': 'Debug',
            'GMX_GPU': 'OFF',
            'GMX_OPENMP': 'OFF',
            'GMX_SIMD': 'None',
            'GMX_USE_RDTSCP': 'OFF'
        }
    release = (context.job_type == JobType.RELEASE)
    if release:
        cmake_opts['GMX_BUILD_TARBALL'] = 'ON'
    elif context.job_type == JobType.GERRIT:
        cmake_opts['GMX_COMPACT_DOXYGEN'] = 'ON'
    cmake_opts.update(context.get_doc_cmake_options(
        doxygen_version='1.8.5', sphinx_version='1.4.1'))
    context.run_cmake(cmake_opts);
    context.build_target(target='gmx', parallel=True,
            continue_on_failure=True)

    context.build_target(target='manual', parallel=False,
            target_descr='PDF manual', continue_on_failure=True)
    logfile = os.path.join(context.workspace.build_dir, 'docs/manual/gromacs.log')
    if os.path.isfile(logfile):
        with open(logfile, 'r') as f:
            manual_log = f.read()
        if re.search(r'LaTeX Warning: Reference .* on page .* undefined', manual_log):
            context.mark_unstable('undefined references in PDF manual')
    context.publish_logs([logfile])

    # check-source is not necessary for a release build, and building these
    # separately causes many of the Doxygen targets to get built twice if run
    # from a tarball.
    if not release:
        context.build_target(target='doxygen-all', parallel=False,
                target_descr='Doxygen documentation', continue_on_failure=True)
        context.build_target(target='check-source', parallel=False,
                failure_string='check-source failed to run', continue_on_failure=True)
        logs = []
        for target in ('check-source', 'doxygen-xml', 'doxygen-user',
                'doxygen-lib', 'doxygen-full'):
            logfile = os.path.join(context.workspace.build_dir,
                    'docs/doxygen/{0}.log'.format(target))
            if os.path.isfile(logfile) and os.stat(logfile).st_size > 0:
                context.mark_unstable('{0} produced warnings'.format(target))
            logs.append(logfile)
        context.publish_logs(logs, category='doxygen')
        if context.failed:
            return

    context.build_target(target='sphinx-input', parallel=False,
            failure_string='Generating Sphinx input failed',
            continue_on_failure=True)
    context.build_target(target='sphinx-programs', parallel=False,
            failure_string='Running gmx help -export rst failed',
            continue_on_failure=True)
    if context.failed:
        return

    sphinx_targets = [ ('webpage-sphinx', 'html', 'HTML') ]
    if not release:
        sphinx_targets.extend((
                ('man', 'man', 'man page'),
                ('install-guide', 'install', 'install-guide')
            ))
    logs = []
    for target, log, descr in sphinx_targets:
        context.build_target(target=target, parallel=False,
                failure_string='Sphinx: {0} generation failed'.format(descr),
                continue_on_failure=True)
        logfile = os.path.join(context.workspace.build_dir,
                'docs/sphinx-{0}.log'.format(log))
        if os.path.isfile(logfile) and os.stat(logfile).st_size > 0:
            context.mark_unstable('Sphinx: {0} generation produced warnings'.format(descr))
        logs.append(logfile)
    context.publish_logs(logs, category='sphinx')
    if context.failed:
        return

    context.build_target(target='webpage', parallel=False)
    if context.failed:
        return

    ignore_urls = ['html-full', 'html-user', 'html-lib', '.tar.gz', '_sources']
    cmd = ['linkchecker', 'docs/html/index.html', '-f',
            context.workspace.get_project_dir(Project.GROMACS) + '/docs/linkcheckerrc',
            '-Fxml']
    for url in ignore_urls:
        cmd.extend(['--ignore-url', url])

    context.run_cmd(cmd, ignore_failure=False)

    logfile = os.path.join(context.workspace.build_dir,
        'docs/linkchecker-out.xml')
    if os.path.isfile(logfile):
        with open(logfile, 'r') as f:
            manual_log = f.read()
            if re.search(r'URLError:', manual_log):
                context.mark_unstable('non resolvable URL in webpage')
            if re.search(r'warnings', manual_log):
                context.mark_unstable('empty pages in web documentation')
        context.publish_logs([logfile], category='webpage')

    if context.failed:
        return

    if release:
        version_info = context.read_cmake_variable_file('VersionInfo.cmake')
        version = version_info['GMX_VERSION_STRING']
        package_name = 'website-' + version
        context.make_archive(package_name, root_dir='docs/html', prefix=package_name)
