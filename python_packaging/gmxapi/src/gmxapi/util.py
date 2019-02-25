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
"""Utility functions supporting the Gromacs Python interface.
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

__all__ = ['which']

import os

from gmxapi import exceptions

def which(command):
    """
    Get the full path of an executable that can be resolved by the shell.

    :param command: executable in the user's PATH
    :return: Absolute path of executable.

    Ref: https://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
    """
    try:
        command_path = os.fsencode(command)
    except:
        raise exceptions.ValueError("Argument must be representable on the command line.")
    if os.path.exists(command_path):
        command_path = os.path.abspath(command_path)
        if os.access(command_path, os.X_OK):
            return command_path
    else:
        # Try to find the executable on the default PATH
        try:
            # available in Python 3
            from shutil import which
        except:
            # Python 2 compatibility, from
            # https://stackoverflow.com/a/377028/5351807
            def which(program):
                import os
                is_exe = lambda fpath: os.path.isfile(fpath) and os.access(fpath, os.X_OK)

                fpath, fname = os.path.split(program)
                if fpath:
                    if is_exe(program):
                        return program
                else:
                    for path in os.environ["PATH"].split(os.pathsep):
                        exe_file = os.path.join(path, program)
                        if is_exe(exe_file):
                            return exe_file
                return None

        return which(command)

def _test():
    import doctest
    doctest.testmod()

if __name__ == "__main__":
    _test()
