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

"""Test the utilities in the gmxapi.runtime module."""

import logging

import gmxapi.utility
from gmxapi.runtime import (
    scoped_resources,
    BaseContext,
    ResourceRequirements,
    ResourceAssignment,
)

try:
    from mpi4py import MPI

    rank_number = MPI.COMM_WORLD.Get_rank()
    comm_size = MPI.COMM_WORLD.Get_size()
except ImportError:
    rank_number = 0
    comm_size = 1
    rank_tag = ""
    MPI = None
else:
    rank_tag = "rank{}:".format(rank_number)

logger = logging.getLogger()


def test_no_requested_comm():
    """Indicate multiple subtasks, but don't request any comm resources."""
    base = BaseContext.instance()
    assert base.communicator().Get_rank() == rank_number
    assert base.communicator().Get_size() == comm_size
    assignment = None
    # Indicate multiple subtasks, but don't request any comm resources.
    for num_tasks in range(1, comm_size + 1):
        requirements = ResourceRequirements(
            communication=tuple({} for _ in range(num_tasks))
        )
        with scoped_resources(base, requirements=requirements) as assignment:
            assert assignment._base_comm_size == base.communicator().Get_size()
            assert assignment.base_rank() == base.communicator().Get_rank()
            assert set(assignment._group_ids) == set(range(num_tasks))
            group_id = assignment._group_ids[rank_number]
            assert group_id >= 0
            assert group_id < num_tasks
            roots = assignment._roots
            assert len(set(roots) - {None}) == num_tasks
            assert set(roots) - {None} == set(range(num_tasks))
            if rank_number == 0:
                assert group_id == 0
            if rank_number == comm_size - 1:
                assert group_id == num_tasks - 1
            if assignment.is_subtask_root():
                _count = 1
            else:
                _count = 0
            num_roots = base.communicator().allreduce(_count)
            assert num_roots == num_tasks
            assert not assignment.subtask_id()
            assert not assignment.subtask_rank()
            assert set(assignment._group_ids) == set(range(num_tasks))
    assert isinstance(assignment, ResourceAssignment)
    message = " ".join(
        [
            f"*{rank}*" if rank == assignment.base_rank() else f"{rank}"
            for rank in range(comm_size)
        ]
    )
    logger.info(message)


def test_no_ensemble_comm():
    """Indicate multiple subtasks, but don't request an ensemble comm."""
    base = BaseContext.instance()
    for num_tasks in range(1, comm_size + 1):
        requirements = ResourceRequirements(
            communication=tuple({"subtask_comm": True} for _ in range(num_tasks))
        )
        with scoped_resources(base, requirements=requirements) as assignment:
            if assignment.subtask_rank() == 0:
                assert assignment.is_subtask_root()
            subtasks = base.communicator().allgather(assignment.subtask_id())
            subtasks = set(subtasks) - {None}
            assert subtasks == set(range(num_tasks))
            subtask_communicator = assignment.communicator()
            if assignment.subtask_rank() is not None:
                assert subtask_communicator.Get_rank() == assignment.subtask_rank()
                _count = 1
            else:
                _count = 0
            total = base.communicator().allreduce(_count)
            if gmxapi.utility.config()["gmx_mpi_type"] == "library":
                assert total == comm_size
            else:
                if subtask_communicator:
                    assert subtask_communicator.Get_size() == 1
                assert total == num_tasks
            if assignment.subtask_id():
                total = subtask_communicator.allreduce(_count)
                if gmxapi.utility.config()["gmx_mpi_type"] == "library":
                    assert total >= 1
                    assert total <= comm_size - (num_tasks - 1)
                else:
                    assert total == 1
                    assert assignment.subtask_rank() == 0


def test_no_subtask_comm():
    """Indicate multiple subtasks, but don't request subtask comms."""
    base = BaseContext.instance()
    for num_tasks in range(1, comm_size + 1):
        requirements = ResourceRequirements(
            communication=tuple({"ensemble_comm": 1} for _ in range(num_tasks))
        )
        with scoped_resources(base, requirements=requirements) as assignment:
            assert not assignment.subtask_id()
            assert not assignment.subtask_rank()
            assert assignment._resource_id is not None
            labels = base.communicator().allgather(assignment._resource_id)
            assert len(set(labels)) == 1
            if assignment.is_subtask_root():
                assert assignment.ensemble_communicator()
            if assignment.ensemble_communicator():
                assert assignment.ensemble_communicator().Get_size() == num_tasks
                assert assignment.is_subtask_root()
                assert (
                    assignment._group_ids[rank_number]
                    == assignment.ensemble_communicator().Get_rank()
                )
                _count = 1
                assert assignment.ensemble_communicator().allreduce(_count) == num_tasks
            else:
                _count = 0
            total = base.communicator().allreduce(_count)
            assert total == num_tasks


def test_assign_ensemble():
    """Indicate an ensemble of subtasks with all resources."""
    base = BaseContext.instance()
    # for num_tasks in range(3, 4):
    for num_tasks in range(1, comm_size + 1):
        requirements = ResourceRequirements(
            communication=tuple(
                {"ensemble_comm": 2, "subtask_comm": True} for _ in range(num_tasks)
            )
        )
        with scoped_resources(base, requirements=requirements) as assignment:
            subtask_rank = assignment.subtask_rank()
            if subtask_rank is not None:
                subtask_comm = assignment.communicator()
                assert subtask_comm.Get_rank() == subtask_rank
                assert subtask_comm.Get_size() > subtask_rank
            if subtask_rank == 0:
                assert bool(assignment.ensemble_communicator())
            if assignment.ensemble_communicator():
                assert subtask_rank == 0
                assert (
                    assignment.subtask_id()
                    == assignment.ensemble_communicator().Get_rank()
                )
