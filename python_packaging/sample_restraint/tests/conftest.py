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

"""Configuration and fixtures for pytest."""

import json
import logging
import os

import pytest

pytest_plugins = ("gmxapi.testsupport",)


@pytest.fixture(scope="class")
def spc_water_box(gmxcli, tmp_path_factory):
    """Provide a TPR input file for a simple simulation.

    Prepare the MD input in a freshly created working directory.
    """
    import gmxapi.runtime
    from gmxapi.testsupport import scoped_chdir

    # TODO: (#2896) Fetch MD input from package / library data.
    # Example:
    #     import pkg_resources
    #     # Note: importing pkg_resources means setuptools is required for running this test.
    #     # Get or build TPR file from data bundled via setup(package_data=...)
    #     # Ref https://setuptools.readthedocs.io/en/latest/setuptools.html#including-data-files
    #     from gmx.data import tprfilename

    with scoped_chdir(tmp_path_factory.mktemp("spc_water_box")) as tempdir:

        testdir = os.path.dirname(__file__)
        with open(os.path.join(testdir, "testdata.json"), "r") as fh:
            testdata = json.load(fh)

        # TODO: (#2756) Don't rely on so many automagical behaviors (as described in comments below)

        structurefile = os.path.join(tempdir, "structure.gro")
        # We let `gmx solvate` use the default solvent. Otherwise, we would do
        #     gro_input = testdata['solvent_structure']
        #     with open(structurefile, 'w') as fh:
        #         fh.write('\n'.join(gro_input))
        #         fh.write('\n')

        topfile = os.path.join(tempdir, "topology.top")
        top_input = testdata["solvent_topology"]
        # `gmx solvate` will append a line to the provided file with the molecule count,
        # so we strip the last line from the input topology.
        with open(topfile, "w") as fh:
            fh.write("\n".join(top_input[:-1]))
            fh.write("\n")

        assert os.path.exists(topfile)
        assert gmxapi.version.api_is_at_least(0, 3, 1)
        solvate = gmxapi.commandline_operation(
            gmxcli,
            arguments=["solvate", "-box", "5", "5", "5"],
            # We use the default solvent instead of specifying one.
            # input_files={'-cs': structurefile},
            output_files={
                "-p": topfile,
                "-o": structurefile,
            },
        )
        assert os.path.exists(topfile)

        if solvate.output.returncode.result() != 0:
            logging.debug(solvate.output.stderr.result())
            raise RuntimeError("solvate failed in spc_water_box testing fixture.")

        # Choose an exactly representable dt of 2^-9 ps (approximately 0.002)
        dt = 2.0**-9.0
        mdp_input = [
            ("integrator", "md"),
            ("dt", dt),
            ("cutoff-scheme", "Verlet"),
            ("nsteps", 2),
            ("nstxout", 1),
            ("nstvout", 1),
            ("nstfout", 1),
            ("tcoupl", "v-rescale"),
            ("tc-grps", "System"),
            ("tau-t", 1),
            ("ref-t", 298),
        ]
        mdp_input = "\n".join(
            [" = ".join([str(item) for item in kvpair]) for kvpair in mdp_input]
        )
        mdpfile = os.path.join(tempdir, "md.mdp")
        with open(mdpfile, "w") as fh:
            fh.write(mdp_input)
            fh.write("\n")
        tprfile = os.path.join(tempdir, "topol.tpr")
        # We don't use mdout_mdp, but if we don't specify it to grompp,
        # it will be created in the current working directory.
        mdout_mdp = os.path.join(tempdir, "mdout.mdp")

        grompp = gmxapi.commandline_operation(
            gmxcli,
            "grompp",
            input_files={
                "-f": mdpfile,
                "-p": solvate.output.file["-p"],
                "-c": solvate.output.file["-o"],
                "-po": mdout_mdp,
            },
            output_files={"-o": tprfile},
        )
        tprfilename = grompp.output.file["-o"].result()
        if grompp.output.returncode.result() != 0:
            logging.debug(grompp.output.stderr.result())
            raise RuntimeError("grompp failed in spc_water_box testing fixture.")

        # TODO: more inspection of grompp errors...
        assert os.path.exists(tprfilename)
        yield tprfilename
