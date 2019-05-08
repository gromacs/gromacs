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

"""Test gmxapi functionality described in roadmap.rst."""

import os
import shutil
import stat
import tempfile

import gmxapi as gmx
import pytest
from gmxapi.version import has_feature


@pytest.mark.skipif(not has_feature('fr3'),
                    reason="Feature level not met.")
def test_fr3():
    """FR3: Output proxy can be used as input."""
    with tempfile.TemporaryDirectory() as directory:
        file1 = os.path.join(directory, 'input')
        file2 = os.path.join(directory, 'output')

        # Make a shell script that acts like the type of tool we are wrapping.
        scriptname = os.path.join(directory, 'clicommand.sh')
        with open(scriptname, 'w') as fh:
            fh.write('\n'.join(['#!' + shutil.which('bash'),
                                '# Concatenate an input file and a string argument to an output file.',
                                '# Mock a utility with the tested syntax.',
                                '#     clicommand.sh "some words" -i inputfile -o outputfile',
                                'echo $1 | cat $3 - > $5\n']))
        os.chmod(scriptname, stat.S_IRWXU)

        line1 = 'first line'
        filewriter1 = gmx.commandline_operation(scriptname,
                                                arguments=[line1],
                                                input_files={'-i': os.devnull},
                                                output_files={'-o': file1})

        line2 = 'second line'
        filewriter2 = gmx.commandline_operation(scriptname,
                                                arguments=[line2],
                                                input_files={'-i': filewriter1.output.file['-o']},
                                                output_files={'-o': file2})

        filewriter2.run()
        # Check that the files have the expected lines
        with open(file1, 'r') as fh:
            lines = [text.rstrip() for text in fh]
        assert len(lines) == 1
        assert lines[0] == line1
        with open(file2, 'r') as fh:
            lines = [text.rstrip() for text in fh]
        assert len(lines) == 2
        assert lines[0] == line1
        assert lines[1] == line2
