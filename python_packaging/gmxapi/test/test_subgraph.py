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
import os

import pytest

import gmxapi as gmx

try:
    from mpi4py import MPI

    comm = MPI.COMM_WORLD
    rank_number = comm.Get_rank()
    comm_size = comm.Get_size()
except ImportError:
    comm = None
    rank_number = 0
    comm_size = 1
    rank_tag = ""
    MPI = None
else:
    rank_tag = "rank{}:".format(rank_number)

mpi_support = pytest.mark.skipif(
    comm_size > 1
    and gmx.utility.config()["gmx_mpi_type"] == "library"
    and not gmx.version.has_feature("mpi_comm_integration"),
    reason="Multi-rank MPI contexts require gmxapi 0.4.",
)


@gmx.function_wrapper(output={"data": float})
def add_float(a: float, b: float) -> float:
    return a + b


@gmx.function_wrapper(output={"data": bool})
def less_than(lhs: float, rhs: float, output=None):
    output.data = lhs < rhs


def test_subgraph_function():
    subgraph = gmx.subgraph(variables={"float_with_default": 1.0, "bool_data": True})
    with subgraph:
        # Define the update for float_with_default to come from an add_float operation.
        subgraph.float_with_default = add_float(
            subgraph.float_with_default, 1.0
        ).output.data
        subgraph.bool_data = less_than(
            lhs=subgraph.float_with_default, rhs=6.0
        ).output.data
    operation_instance = subgraph()
    operation_instance.run()
    assert operation_instance.values["float_with_default"] == 2.0

    loop = gmx.while_loop(operation=subgraph, condition=subgraph.bool_data)
    handle = loop()
    assert handle.output.float_with_default.result() == 6


_previous_cpt = [""] * comm_size


@gmx.function_wrapper(output={"data": bool})
def _is_new(current_cpt: str, output):
    # Is the simulation new or is it continuing from a checkpoint?
    global _previous_cpt
    # TODO: Don't use private details of the publishing data proxy.
    # Either use another subgraph variable for previous_cpt, distribute function_wrapper tasks by MPI rank,
    # or develop the interface for session resources provided to tasks.
    executing_member = output._client_identifier
    previous_cpt = _previous_cpt[executing_member]
    _previous_cpt[executing_member] = current_cpt
    if current_cpt != previous_cpt and os.path.exists(previous_cpt):
        output.data = False
    else:
        output.data = True


@mpi_support
def test_subgraph_simulation_extension(spc_water_box, mdrun_kwargs):
    tpr_list = gmx.read_tpr([spc_water_box] * comm_size)
    input_list = gmx.modify_input(tpr_list, parameters={"nsteps": 10**6})
    subgraph = gmx.subgraph(
        variables={
            "new": True,
            "checkpoint": "",
            # 'previous': '',
        }
    )
    with subgraph:
        md = gmx.mdrun(
            input_list,
            runtime_args={
                "-cpi": subgraph.checkpoint,
                "-maxh": "0.001",
                "-noappend": None,
            },
        )

        subgraph.new = _is_new(md.output.checkpoint).output.data

        subgraph.checkpoint = md.output.checkpoint

    trajectory_continuation_loop = gmx.while_loop(
        operation=subgraph, condition=subgraph.new
    )()

    trajectory_continuation_loop.run()
    _cpt_output = trajectory_continuation_loop.output.checkpoint
    final_checkpoint = _cpt_output.result()
    loop_condition = trajectory_continuation_loop.output.new.result()
    if comm_size > 1:
        final_checkpoint = final_checkpoint[rank_number]
        loop_condition = loop_condition[rank_number]
    assert isinstance(final_checkpoint, str)
    assert isinstance(loop_condition, bool)
    assert os.path.exists(final_checkpoint)


def test_local_tools_and_assumptions():
    const = gmx.make_constant(1.0)
    assert add_float(const, const).output.data.result() == 2
    assert gmx.logical_not(less_than(const, const).output.data).result()
    # Note: It may not be safe to assume that keyword argument order (lhs, rhs) is preserved.
    assert less_than(const, add_float(const, const).output.data).output.data.result()
