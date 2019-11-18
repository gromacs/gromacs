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

def do_build(context):
    context.env.set_env_var('CLANG_FORMAT', context.env.get_clang_format_command('7'))
    clangformat_log = context.workspace.get_path_for_logfile('clang-format.log', category='clang-format')
    copyright_log  = context.workspace.get_path_for_logfile('copyright.log', category='copyright')
    cmd = ['admin/clang-format.sh', 'check', '--rev=HEAD^', '--warnings=' + clangformat_log]
    ret = context.run_cmd(cmd, use_return_code=True)
    if ret == 1:
        with open(clangformat_log, 'r') as f:
            warnings = f.readlines()
        if len(warnings) <= 5:
            details = [x.rstrip() for x in warnings]
        else:
            format_count = 0
            for w in warnings:
                if 'needs formatting' in w:
                    format_count += 1
            details = []
            if format_count > 0:
                details.append('formatting issues in {0} files'.format(format_count))
        context.mark_unstable(reason='clang-format.sh found issues', details=details)
    elif ret != 0:
        raise BuildError('clang-format.sh failed to run')

    cmd = ['admin/copyright.sh', 'check', '--rev=HEAD^', '--warnings=' + copyright_log]
    ret = context.run_cmd(cmd, use_return_code=True)
    if ret == 1:
        with open(copyright_log, 'r') as f:
            warnings = f.readlines()
        if len(warnings) <= 5:
            details = [x.rstrip() for x in warnings]
        else:
            cpyear_count = 0
            cpheader_count = 0
            for w in warnings:
                if 'copyright year' in w:
                    cpyear_count += 1
                if 'copyright header' in w:
                    cpheader_count += 1
            details = []
            if cpyear_count > 0:
                details.append('copyright year missing in {0} files'.format(cpyear_count))
            if cpheader_count > 0:
                details.append('copyright header issues in {0} files'.format(cpheader_count))
        context.mark_unstable(reason='copyright.sh found issues', details=details)
    elif ret != 0:
        raise BuildError('copyright.sh failed to run')
