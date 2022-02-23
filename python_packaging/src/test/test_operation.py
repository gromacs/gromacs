#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2019- The GROMACS Authors
# and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
# Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
# https://www.gnu.org/licenses, or write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
#
# If you want to redistribute modifications to GROMACS, please
# consider that scientific software is very special. Version
# control is crucial - bugs must be traceable. We will be happy to
# consider code for inclusion in the official distribution, but
# derived work must not be called official GROMACS. Details are found
# in the README & COPYING files - if they are missing, get the
# official version at https://www.gromacs.org.
#
# To help us fund GROMACS development, we humbly ask that you cite
# the research papers on the package. Check out https://www.gromacs.org.

"""Tests the interfaces defined in operation.py and behavior of simple operations.

See also test_commandline.py for integration tests.
"""

import os
import shutil
import stat
import tempfile
from typing import NamedTuple

import gmxapi as gmx
from gmxapi import commandline_operation
from gmxapi.operation import ResultDescription


def test_comparison():
    int_array_description = ResultDescription(
        dtype=int,
        width=3
    )
    same_description = ResultDescription(
        dtype=int,
        width=3
    )
    different_type_description = ResultDescription(
        dtype=float,
        width=3
    )
    different_width_description = ResultDescription(
        dtype=int,
        width=1
    )
    assert int_array_description == same_description
    assert int_array_description != different_width_description
    assert int_array_description != different_type_description

    class CompatibleRepresentation(NamedTuple):
        dtype: type
        width: int

    equivalent_description = CompatibleRepresentation(dtype=int, width=3)
    assert int_array_description == equivalent_description
    assert equivalent_description == int_array_description
    assert equivalent_description != different_type_description
    assert different_type_description != equivalent_description
    assert equivalent_description != different_width_description
    assert different_width_description != equivalent_description


def test_scalar():
    operation = gmx.make_constant(42)
    assert isinstance(operation.dtype, type)
    assert operation.dtype == int
    assert operation.result() == 42


def test_list():
    list_a = [1, 2, 3]

    # TODO: test input validation
    list_result = gmx.concatenate_lists(sublists=[list_a])
    assert list_result.dtype == gmx.datamodel.NDArray
    # Note: this is specifically for the built-in tuple type.
    # Equality comparison may work differently for different sequence types.
    assert tuple(list_result.result()) == tuple(list_a)
    assert len(list_result.result()) == len(list_a)

    list_result = gmx.concatenate_lists([list_a, list_a])
    assert len(list_result.result()) == len(list_a) * 2
    assert tuple(list_result.result()) == tuple(list_a + list_a)

    list_b = gmx.ndarray([42])

    list_result = gmx.concatenate_lists(sublists=[list_b])
    assert list_result.result()[0] == 42

    list_result = gmx.join_arrays(front=list_a, back=list_b)
    assert len(list_result.result()) == len(list_a) + 1
    assert tuple(list_result.result()) == tuple(list(list_a) + [42])


def test_data_dependence(cleandir):
    """Test dependent sequence of operations.

    Confirm that data dependencies correctly establish resolvable execution dependencies.

    In a sequence of two operations, write a two-line file one line at a time.
    Use the output of one operation as the input of another.
    """
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
        filewriter1 = commandline_operation(scriptname,
                                            arguments=[line1],
                                            input_files={'-i': os.devnull},
                                            output_files={'-o': file1})

        line2 = 'second line'
        filewriter2 = commandline_operation(scriptname,
                                            arguments=[line2],
                                            input_files={'-i': filewriter1.output.file['-o']},
                                            output_files={'-o': file2})

        filewriter2.run()
        # Check that the files have the expected lines
        with open(filewriter1.output.file['-o'].result(), 'r') as fh:
            lines = [text.rstrip() for text in fh]
        assert len(lines) == 1
        assert lines[0] == line1
        with open(filewriter2.output.file['-o'].result(), 'r') as fh:
            lines = [text.rstrip() for text in fh]
        assert len(lines) == 2
        assert lines[0] == line1
        assert lines[1] == line2

        try:
            # If we're running with MPI, check that the file is accessible
            # on all ranks, but make sure that we don't delete the temporary
            # directory before all ranks have looked for it.
            from mpi4py import MPI
            MPI.COMM_WORLD.barrier()
        except ImportError:
            pass
