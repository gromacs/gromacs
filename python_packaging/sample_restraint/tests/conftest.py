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

"""Configuration and fixtures for pytest."""

import json
import logging
import os
import shutil
import tempfile
import warnings
from contextlib import contextmanager

import pytest


def pytest_addoption(parser):
    """Add a command-line user option for the pytest invocation."""
    parser.addoption(
        '--rm',
        action='store',
        default='always',
        choices=['always', 'never', 'success'],
        help='Remove temporary directories "always", "never", or on "success".'
    )


@pytest.fixture(scope='session')
def remove_tempdir(request):
    """pytest fixture to get access to the --rm CLI option."""
    return request.config.getoption('--rm')


@contextmanager
def scoped_chdir(dir):
    oldpath = os.getcwd()
    os.chdir(dir)
    try:
        yield dir
        # If the `with` block using scoped_chdir produces an exception, it will
        # be raised at this point in this function. We want the exception to
        # propagate out of the `with` block, but first we want to restore the
        # original working directory, so we skip `except` but provide a `finally`.
    finally:
        os.chdir(oldpath)


@contextmanager
def _cleandir(remove_tempdir):
    """Context manager for a clean temporary working directory.

    Arguments:
        remove_tempdir (str): whether to remove temporary directory "always",
                              "never", or on "success"

    The context manager will issue a warning for each temporary directory that
    is not removed.
    """

    newpath = tempfile.mkdtemp()

    def remove():
        shutil.rmtree(newpath)

    def warn():
        warnings.warn('Temporary directory not removed: {}'.format(newpath))

    if remove_tempdir == 'always':
        callback = remove
    else:
        callback = warn
    try:
        with scoped_chdir(newpath):
            yield newpath
        # If we get to this line, the `with` block using _cleandir did not throw.
        # Clean up the temporary directory unless the user specified `--rm never`.
        # I.e. If the user specified `--rm success`, then we need to toggle from `warn` to `remove`.
        if remove_tempdir != 'never':
            callback = remove
    finally:
        callback()


@pytest.fixture
def cleandir(remove_tempdir):
    """Provide a clean temporary working directory for a test.

    Example usage:

        import os
        import pytest

        @pytest.mark.usefixtures("cleandir")
        def test_cwd_starts_empty():
            assert os.listdir(os.getcwd()) == []
            with open("myfile", "w") as f:
                f.write("hello")

        def test_cwd_also_starts_empty(cleandir):
            assert os.listdir(os.getcwd()) == []
            assert os.path.abspath(os.getcwd()) == os.path.abspath(cleandir)
            with open("myfile", "w") as f:
                f.write("hello")

        @pytest.mark.usefixtures("cleandir")
        class TestDirectoryInit(object):
            def test_cwd_starts_empty(self):
                assert os.listdir(os.getcwd()) == []
                with open("myfile", "w") as f:
                    f.write("hello")

            def test_cwd_also_starts_empty(self):
                assert os.listdir(os.getcwd()) == []
                with open("myfile", "w") as f:
                    f.write("hello")

    Ref: https://docs.pytest.org/en/latest/fixture.html#using-fixtures-from-classes-modules-or-projects
    """
    with _cleandir(remove_tempdir) as newdir:
        yield newdir


@pytest.fixture(scope='session')
def gmxcli():
    # TODO: (#2896) Find a more canonical way to identify the GROMACS commandline wrapper binary.
    #  We should be able to get the GMXRC contents and related hints from a gmxapi
    #  package resource or from module attributes of a ``gromacs`` stub package.
    allowed_command_names = ['gmx', 'gmx_mpi']
    command = None
    for command_name in allowed_command_names:
        if command is not None:
            break
        command = shutil.which(command_name)
        if command is None:
            gmxbindir = os.getenv('GMXBIN')
            if gmxbindir is None:
                gromacsdir = os.getenv('GROMACS_DIR')
                if gromacsdir is not None and gromacsdir != '':
                    gmxbindir = os.path.join(gromacsdir, 'bin')
            if gmxbindir is None:
                gmxapidir = os.getenv('gmxapi_DIR')
                if gmxapidir is not None and gmxapidir != '':
                    gmxbindir = os.path.join(gmxapidir, 'bin')
            if gmxbindir is not None:
                gmxbindir = os.path.abspath(gmxbindir)
                command = shutil.which(command_name, path=gmxbindir)
    if command is None:
        message = "Tests need 'gmx' command line tool, but could not find it on the path."
        raise RuntimeError(message)
    try:
        assert os.access(command, os.X_OK)
    except Exception as E:
        raise RuntimeError('"{}" is not an executable gmx wrapper program'.format(command)) from E
    yield command


@pytest.fixture(scope='class')
def spc_water_box(gmxcli, remove_tempdir):
    """Provide a TPR input file for a simple simulation.

    Prepare the MD input in a freshly created working directory.
    """
    import gmxapi as gmx

    # TODO: (#2896) Fetch MD input from package / library data.
    # Example:
    #     import pkg_resources
    #     # Note: importing pkg_resources means setuptools is required for running this test.
    #     # Get or build TPR file from data bundled via setup(package_data=...)
    #     # Ref https://setuptools.readthedocs.io/en/latest/setuptools.html#including-data-files
    #     from gmx.data import tprfilename

    with _cleandir(remove_tempdir) as tempdir:

        testdir = os.path.dirname(__file__)
        with open(os.path.join(testdir, 'testdata.json'), 'r') as fh:
            testdata = json.load(fh)

        # TODO: (#2756) Don't rely on so many automagical behaviors (as described in comments below)

        structurefile = os.path.join(tempdir, 'structure.gro')
        # We let `gmx solvate` use the default solvent. Otherwise, we would do
        #     gro_input = testdata['solvent_structure']
        #     with open(structurefile, 'w') as fh:
        #         fh.write('\n'.join(gro_input))
        #         fh.write('\n')

        topfile = os.path.join(tempdir, 'topology.top')
        top_input = testdata['solvent_topology']
        # `gmx solvate` will append a line to the provided file with the molecule count,
        # so we strip the last line from the input topology.
        with open(topfile, 'w') as fh:
            fh.write('\n'.join(top_input[:-1]))
            fh.write('\n')

        assert os.path.exists(topfile)
        solvate = gmx.commandline_operation(gmxcli,
                                            arguments=['solvate', '-box', '5', '5', '5'],
                                            # We use the default solvent instead of specifying one.
                                            # input_files={'-cs': structurefile},
                                            output_files={'-p': topfile,
                                                          '-o': structurefile,
                                                          }
                                            )
        assert os.path.exists(topfile)

        if solvate.output.returncode.result() != 0:
            logging.debug(solvate.output.erroroutput.result())
            raise RuntimeError('solvate failed in spc_water_box testing fixture.')

        # Choose an exactly representable dt of 2^-9 ps (approximately 0.002)
        dt = 2.**-9.
        mdp_input = [('integrator', 'md'),
                     ('dt', dt),
                     ('cutoff-scheme', 'Verlet'),
                     ('nsteps', 2),
                     ('nstxout', 1),
                     ('nstvout', 1),
                     ('nstfout', 1),
                     ('tcoupl', 'v-rescale'),
                     ('tc-grps', 'System'),
                     ('tau-t', 1),
                     ('ref-t', 298)]
        mdp_input = '\n'.join([' = '.join([str(item) for item in kvpair]) for kvpair in mdp_input])
        mdpfile = os.path.join(tempdir, 'md.mdp')
        with open(mdpfile, 'w') as fh:
            fh.write(mdp_input)
            fh.write('\n')
        tprfile = os.path.join(tempdir, 'topol.tpr')
        # We don't use mdout_mdp, but if we don't specify it to grompp,
        # it will be created in the current working directory.
        mdout_mdp = os.path.join(tempdir, 'mdout.mdp')

        grompp = gmx.commandline_operation(gmxcli, 'grompp',
                                           input_files={
                                               '-f': mdpfile,
                                               '-p': solvate.output.file['-p'],
                                               '-c': solvate.output.file['-o'],
                                               '-po': mdout_mdp,
                                           },
                                           output_files={'-o': tprfile})
        tprfilename = grompp.output.file['-o'].result()
        if grompp.output.returncode.result() != 0:
            logging.debug(grompp.output.erroroutput.result())
            raise RuntimeError('grompp failed in spc_water_box testing fixture.')

        # TODO: more inspection of grompp errors...
        assert os.path.exists(tprfilename)
        yield tprfilename
