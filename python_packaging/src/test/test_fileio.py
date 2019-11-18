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

"""Test gmx.fileio submodule"""
import os
import tempfile

import gmxapi
import pytest
from gmxapi.simulation.fileio import TprFile
from gmxapi.simulation.fileio import read_tpr
from gmxapi.exceptions import UsageError


@pytest.mark.usefixtures('cleandir')
def test_tprfile_read_old(spc_water_box):
    tpr_filename = spc_water_box
    with pytest.raises(UsageError):
        TprFile(tpr_filename, 'x')
    with pytest.raises(UsageError):
        TprFile()
    tprfile = TprFile(tpr_filename, 'r')
    with tprfile as fh:
        cpp_object = fh._tprFileHandle
        assert cpp_object is not None
        params = cpp_object.params().extract()
        assert "nsteps" in params
        assert "foo" not in params


@pytest.mark.usefixtures('cleandir')
def test_read_tpr(spc_water_box):
    tpr_filename = spc_water_box
    tpr_source = gmxapi.read_tpr(tpr_filename)
    assert 'parameters' in tpr_source.output
    assert hasattr(tpr_source.output, 'parameters')
    parameters = tpr_source.output.parameters.result()
    assert 'nsteps' in parameters
    assert 'foo' not in parameters
    assert parameters['nsteps'] == 2
    assert tpr_source.output.parameters['nsteps'].result() == 2


@pytest.mark.usefixtures('cleandir')
def test_core_tprcopy_alt(spc_water_box):
    """Test gmx.core.copy_tprfile() for update of end_time.

    Set a new end time that is 5000 steps later than the original. Read dt
    from file to avoid floating point round-off errors.

    Transitively test gmx.fileio.read_tpr()
    """
    tpr_filename = spc_water_box
    additional_steps = 5000
    sim_input = read_tpr(tpr_filename)
    params = sim_input.parameters.extract()
    dt = params['dt']
    nsteps = params['nsteps']
    init_step = params['init-step']
    initial_endtime = (init_step + nsteps) * dt
    new_endtime = initial_endtime + additional_steps*dt
    _, temp_filename = tempfile.mkstemp(suffix='.tpr')
    gmxapi._gmxapi.rewrite_tprfile(source=tpr_filename, destination=temp_filename, end_time=new_endtime)
    tprfile = TprFile(temp_filename, 'r')
    with tprfile as fh:
        params = read_tpr(fh).parameters.extract()
        dt = params['dt']
        nsteps = params['nsteps']
        init_step = params['init-step']
        assert (init_step + nsteps) * dt == new_endtime

    os.unlink(temp_filename)


@pytest.mark.usefixtures('cleandir')
def test_write_tpr_file(spc_water_box):
    """Test gmx.fileio.write_tpr_file() using gmx.core API.
    """
    tpr_filename = spc_water_box
    additional_steps = 5000
    sim_input = read_tpr(tpr_filename)
    params = sim_input.parameters.extract()
    nsteps = params['nsteps']
    init_step = params['init-step']
    # Choose a new nsteps to check integer parameter setting.
    new_nsteps = init_step + additional_steps
    # Choose a new dt to check floating point parameter setting
    new_dt = params['dt'] * 2.

    sim_input.parameters.set('nsteps', new_nsteps)
    sim_input.parameters.set('dt', new_dt)

    _, temp_filename = tempfile.mkstemp(suffix='.tpr')
    gmxapi.simulation.fileio.write_tpr_file(temp_filename, input=sim_input)
    tprfile = TprFile(temp_filename, 'r')
    with tprfile as fh:
        params = read_tpr(fh).parameters.extract()
        # Note that we have chosen an exactly representable dt for spc_water_box.
        # Otherwise, we would have to use pytest.approx with a suitable tolerance.
        assert params['dt'] == new_dt
        assert params['nsteps'] != nsteps
        assert params['nsteps'] == new_nsteps

    os.unlink(temp_filename)
