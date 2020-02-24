#!/usr/bin/env python
#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2019,2020, by the GROMACS development team, led by
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

"""Tests for dynamically defined operations.

Test the command-line wrapping functionality in gmxapi.commandline. These
operation factories are written with user-facing tools and exercise a lot of
the high-level machinery of the package, effectively serving as integration
tests of the operation-building utilities in the modules depended on by
commandline.py.
"""

import shutil
import unittest

from gmxapi import commandline


class SimpleCliTestCase(unittest.TestCase):
    """Test creation and execution of the basic cli() command line wrapper."""

    def test_true(self):
        """Test a command known to produce a return code of 0."""
        command = shutil.which('true')
        operation = commandline.cli(command=[command], shell=False)

        # Note: 'stdout' and 'stderr' not mapped.
        # Note: getitem not implemented.
        # assert not 'stdout' in operation.output
        # assert not 'stderr' in operation.output
        assert not hasattr(operation.output, 'stdout')
        assert not hasattr(operation.output, 'stderr')

        # Check for the attributes that we _do_ expect.
        assert hasattr(operation.output, 'erroroutput')
        assert hasattr(operation.output, 'returncode')

        operation.run()
        # assert operation.output.returncode.result() == 0
        assert operation.output.returncode.result() == 0

    def test_false_explicit(self):
        """Test a command known to produce a return code of 1."""
        command = shutil.which('false')
        operation = commandline.cli(command=[command], shell=False)
        # Explicitly run the operation.
        operation.run()
        assert operation.output.returncode.result() == 1

    def test_false_implicit(self):
        command = shutil.which('false')
        operation = commandline.cli(command=[command], shell=False)
        # Allow the operation to be executed implicitly to satisfy data constraint.
        assert operation.output.returncode.result() == 1

    def test_command_with_arguments(self):
        """Test that cli() can wrap a command with arguments."""
        # TODO: (FR5+) do we want to pipeline or checkpoint stdout somehow?
        operation = commandline.cli(command=[shutil.which('echo'), 'hi', 'there'], shell=False)
        assert operation.output.returncode.result() == 0

    def test_command_with_stdin(self):
        """Test that cli() can handle string input."""
        stdin = 'hi\nthere\n'
        subcommand = '{wc} -l | {grep} -q 2'.format(wc=shutil.which('wc'), grep=shutil.which('grep'))

        operation = commandline.cli(command=['/bin/sh', '-c', subcommand], shell=False, stdin=stdin)
        assert operation.output.returncode.result() == 0
        operation = commandline.commandline_operation('/bin/sh', ['-c', subcommand], stdin=stdin)
        assert operation.output.returncode.result() == 0

        subcommand = '{wc} -l | {grep} -q 1'.format(wc=shutil.which('wc'), grep=shutil.which('grep'))

        operation = commandline.cli(command=['/bin/sh', '-c', subcommand], shell=False, stdin=stdin)
        assert operation.output.returncode.result() != 0
        operation = commandline.commandline_operation('/bin/sh', ['-c', subcommand], stdin=stdin)
        assert operation.output.returncode.result() != 0


class CommandLineOperationSimpleTestCase(unittest.TestCase):
    """Test the command line wrapper operation factory."""

    def test_true(self):
        operation = commandline.commandline_operation(executable='true')
        # Note: 'stdout' and 'stderr' not mapped.
        # Note: getitem not implemented.
        # assert not 'stdout' in operation.output
        # assert not 'stderr' in operation.output
        assert not hasattr(operation.output, 'stdout')
        assert not hasattr(operation.output, 'stderr')
        assert hasattr(operation.output, 'file')
        assert hasattr(operation.output, 'erroroutput')
        assert hasattr(operation.output, 'returncode')
        assert operation.output.returncode.result() == 0

    def test_false(self):
        operation = commandline.commandline_operation(executable='false')
        assert operation.output.returncode.result() == 1

    def test_echo(self):
        # TODO: (FR5+) do we want to pipeline or checkpoint stdout somehow?
        operation = commandline.commandline_operation(executable='echo',
                                                      arguments=['hi there'])
        assert operation.output.returncode.result() == 0


if __name__ == '__main__':
    unittest.main()
