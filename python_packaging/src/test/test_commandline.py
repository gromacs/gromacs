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

"""Tests for dynamically defined operations.

Test the command-line wrapping functionality in gmxapi.commandline. These
operation factories are written with user-facing tools and exercise a lot of
the high-level machinery of the package, effectively serving as integration
tests of the operation-building utilities in the modules depended on by
commandline.py.
"""
import logging
import os
import shutil
import stat

import pytest

import gmxapi as gmx
import gmxapi.operation
from gmxapi import commandline
from gmxapi.testsupport import scoped_chdir

try:
    from mpi4py import MPI

    comm = MPI.COMM_WORLD
    rank_number = comm.Get_rank()
    comm_size = comm.Get_size()
except ImportError:
    comm = None
    rank_number = 0
    comm_size = 1
    rank_tag = ''
    MPI = None
else:
    rank_tag = 'rank{}:'.format(rank_number)


# if rank_number == 1:
#     import pydevd_pycharm
#     pydevd_pycharm.settrace('localhost', port=12345, stdoutToServer=True, stderrToServer=True)


def test_true_base(cleandir):
    """Test a command known to produce a return code of 0."""
    work_dir = cleandir
    if comm_size > 1:
        work_dir = comm.bcast(work_dir, root=0)
    if rank_number != 0:
        assert work_dir != cleandir
    with scoped_chdir(work_dir):
        current_dir = os.getcwd()
        logging.warning(f'{current_dir}')
        command = shutil.which('true')
        kwargs = dict(command=[command], shell=False)
        if gmx.version.has_feature('cli_env_kwarg'):
            kwargs['env'] = dict(os.environ)
        operation = commandline.cli(**kwargs)

        # Note: getitem not implemented.
        # assert 'stdout' in operation.output
        # assert 'stderr' in operation.output
        assert hasattr(operation.output, 'stdout')
        assert hasattr(operation.output, 'stderr')
        assert not hasattr(operation.output, 'erroroutput')
        assert hasattr(operation.output, 'returncode')

        operation.run()
        # assert operation.output.returncode.result() == 0
        assert operation.output.returncode.result() == 0


def test_false_explicit(cleandir):
    """Test a command known to produce a return code of 1."""
    command = shutil.which('false')
    kwargs = dict(command=[command], shell=False)
    if gmx.version.has_feature('cli_env_kwarg'):
        kwargs['env'] = dict(os.environ)
    operation = commandline.cli(**kwargs)    # Explicitly run the operation.
    operation.run()
    assert operation.output.returncode.result() == 1


def test_false_implicit(cleandir):
    command = shutil.which('false')
    kwargs = dict(command=[command], shell=False)
    if gmx.version.has_feature('cli_env_kwarg'):
        kwargs['env'] = dict(os.environ)
    operation = commandline.cli(**kwargs)
    # Allow the operation to be executed implicitly to satisfy data constraint.
    assert operation.output.returncode.result() == 1


def test_command_with_arguments(cleandir):
    """Test that cli() can wrap a command with arguments."""
    kwargs = dict(command=[shutil.which('echo'), 'hi', 'there'], shell=False)
    if gmx.version.has_feature('cli_env_kwarg'):
        kwargs['env'] = dict(os.environ)
    operation = commandline.cli(**kwargs)
    assert operation.output.returncode.result() == 0


def test_command_with_stdin(cleandir):
    """Test that cli() can handle string input."""
    stdin = 'hi\nthere\n'
    subcommand = '{wc} -l | {grep} -q 2'.format(wc=shutil.which('wc'), grep=shutil.which('grep'))

    kwargs = dict(command=['/bin/sh', '-c', subcommand], shell=False, stdin=stdin)
    if gmx.version.has_feature('cli_env_kwarg'):
        kwargs['env'] = dict(os.environ)
    operation = commandline.cli(**kwargs)
    assert operation.output.returncode.result() == 0
    operation = commandline.commandline_operation('/bin/sh', ['-c', subcommand], stdin=stdin)
    assert operation.output.returncode.result() == 0

    subcommand = '{wc} -l | {grep} -q 1'.format(wc=shutil.which('wc'), grep=shutil.which('grep'))

    kwargs = dict(command=['/bin/sh', '-c', subcommand], shell=False, stdin=stdin)
    if gmx.version.has_feature('cli_env_kwarg'):
        kwargs['env'] = dict(os.environ)
    operation = commandline.cli(**kwargs)
    assert operation.output.returncode.result() != 0
    operation = commandline.commandline_operation('/bin/sh', ['-c', subcommand], stdin=stdin)
    assert operation.output.returncode.result() != 0


def test_true(cleandir):
    operation = commandline.commandline_operation(executable='true')
    # Note: getitem not implemented.
    # assert 'stdout' in operation.output
    # assert 'stderr' in operation.output
    assert not hasattr(operation.output, 'erroroutput')
    assert hasattr(operation.output, 'file')
    assert hasattr(operation.output, 'stdout')
    assert hasattr(operation.output, 'stderr')
    assert hasattr(operation.output, 'returncode')
    assert operation.output.returncode.result() == 0


def test_false(cleandir):
    operation = commandline.commandline_operation(executable='false')
    assert operation.output.returncode.result() == 1


def test_echo(cleandir):
    operation = commandline.commandline_operation(executable='echo',
                                                  arguments=['hi there'])
    assert operation.output.stdout.result() == 'hi there\n'
    assert operation.output.returncode.result() == 0


def test_file_dependency_chain(cleandir):
    """Test the command line wrapper input/output file handling.

    Operation output can be used as operation input.
    """
    file1 = os.path.join(cleandir, 'input')
    file2 = os.path.join(cleandir, 'output')

    # Make a shell script that acts like the type of tool we are wrapping.
    scriptname = os.path.join(cleandir, 'clicommand.sh')
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
    assert isinstance(filewriter1.output, gmxapi.operation.DataProxyBase)
    assert filewriter1.output.ensemble_width == 1

    line2 = 'second line'
    filewriter2 = gmx.commandline_operation(scriptname,
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

    # Make sure that the temporary directory is not removed before all ranks have done
    # the file checks.
    if comm_size > 1:
        comm.barrier()


def test_file_dependency_ensemble_chain(cleandir):
    """Test the command line wrapper input/output file handling.

    Operation output can be used as operation input.
    """
    file1 = (os.path.join(cleandir, 'input1'), os.path.join(cleandir, 'input2'))
    file2 = 'output'

    # Make a shell script that acts like the type of tool we are wrapping.
    scriptname = os.path.join(cleandir, 'clicommand.sh')
    with open(scriptname, 'w') as fh:
        fh.write('\n'.join(['#!' + shutil.which('bash'),
                            '# Concatenate an input file and a string argument to an output file.',
                            '# Mock a utility with the tested syntax.',
                            '#     clicommand.sh "some words" -i inputfile -o outputfile',
                            'echo $1 | cat $3 - > $5\n']))
    os.chmod(scriptname, stat.S_IRWXU)

    line1 = ['first line A', 'first line B']
    filewriter1 = gmx.commandline_operation(scriptname,
                                            arguments=[[line] for line in line1],
                                            input_files={'-i': os.devnull},
                                            output_files=[{'-o': name} for name in file1])
    assert isinstance(filewriter1.output, gmxapi.operation.DataProxyBase)
    assert filewriter1.output.ensemble_width == 2

    line2 = 'second line'
    filewriter2 = gmx.commandline_operation(scriptname,
                                            arguments=[line2],
                                            input_files={'-i': filewriter1.output.file['-o']},
                                            output_files={'-o': file2})
    assert filewriter2.output.ensemble_width == 2
    filewriter2.run()

    # Check that the files have the expected lines
    for member in range(2):
        with open(filewriter1.output.file['-o'].result()[member], 'r') as fh:
            lines = [text.rstrip() for text in fh]
        assert len(lines) == 1
        assert lines[0] == line1[member]

    outputs = filewriter2.output.file['-o'].result()
    assert len(outputs) == 2
    assert outputs[0] != outputs[1]
    for member in range(2):
        path = outputs[member]
        assert os.path.exists(path)
        with open(path, 'r') as fh:
            lines = [text.rstrip() for text in fh]
            assert len(lines) == 2
            assert lines[0] == line1[member]
            assert lines[1] == line2
    # Make sure that the temporary directory is not removed before all ranks have done
    # the file checks.
    if comm_size > 1:
        comm.barrier()


def test_failure(cleandir):
    """The operation should not deliver file output if the subprocess fails."""
    file1 = os.path.join(cleandir, 'input')
    file2 = os.path.join(cleandir, 'output')

    # Make a shell script that acts like the type of tool we are wrapping.
    scriptname = os.path.join(cleandir, 'clicommand.sh')
    with open(scriptname, 'w') as fh:
        fh.write('\n'.join(['#!' + shutil.which('bash'),
                            '# Concatenate an input file and a string argument to an output file.',
                            '# Mock a utility with the tested syntax.',
                            '#     clicommand.sh "some words" -i inputfile -o outputfile',
                            'exit 1\n']))
    os.chmod(scriptname, stat.S_IRWXU)

    filewriter1 = gmx.commandline_operation(scriptname,
                                            input_files={'-i': os.devnull},
                                            output_files={'-o': file1})

    filewriter2 = gmx.commandline_operation(scriptname,
                                            input_files={'-i': filewriter1.output.file['-o']},
                                            output_files={'-o': file2})

    with pytest.warns(UserWarning, match=r'Trouble getting .*'):
        # We expect this to generate errors, warnings, and log messages.
        # We assert the exception and warning, and suppress the logged error
        # to avoid misleading output from the test suite.

        class _Filter(logging.Filter):
            def filter(self, record: logging.LogRecord):
                if record.name == 'gmxapi.operation' and 'Trouble getting' in record.msg:
                    return 0
                else:
                    return 1

        _filter = _Filter()
        logging.getLogger('gmxapi.operation').addFilter(_filter)

        with pytest.raises(KeyError):
            filewriter2.run()

        logging.getLogger('gmxapi.operation').removeFilter(_filter)
