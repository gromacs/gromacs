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

import json

extra_options = {
    'md5sum': Option.string
}

def do_build(context):
    info_path = 'cmake/gmxVersionInfo.cmake'
    cmd = [context.env.cmake_command, '-P', info_path]
    info_json = context.run_cmd(cmd, use_output=True)
    values = json.loads(info_json)
    old_md5sum = values['regressiontest-md5sum']
    new_md5sum = context.opts.md5sum
    if new_md5sum != old_md5sum:
        context.replace_in_file(info_path, r'set\(REGRESSIONTEST_MD5SUM "(\w*)"',
                lambda x: do_replacement(x, new_md5sum))
        context.workspace.upload_revision(project=Project.GROMACS, file_glob=info_path)

def do_replacement(match, new_md5sum):
    result = match.group(0)
    start = match.start(1) - match.start(0)
    end = match.end(1) - match.end(0)
    return result[:start] + new_md5sum + result[end:]
