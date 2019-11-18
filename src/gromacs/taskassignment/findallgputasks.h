/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017,2018,2019, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
/*! \internal
 * \file
 * \brief Declares routine for collecting all GPU tasks found on ranks of a node.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_taskassignment
 */
#ifndef GMX_TASKASSIGNMENT_FINDALLGPUTASKS_H
#define GMX_TASKASSIGNMENT_FINDALLGPUTASKS_H

#include <vector>

namespace gmx
{

enum class GpuTask;
enum class TaskTarget;
class PhysicalNodeCommunicator;
template<typename T>
class ArrayRef;
//! Container of compute tasks suitable to run on a GPU e.g. on each rank of a node.
using GpuTasksOnRanks = std::vector<std::vector<GpuTask>>;

/*! \brief Returns container of all tasks on this rank
 * that are eligible for GPU execution.
 *
 * \param[in]  haveGpusOnThisPhysicalNode Whether there are any GPUs on this physical node.
 * \param[in]  nonbondedTarget            The user's choice for mdrun -nb for where to assign
 *                                        short-ranged nonbonded interaction tasks.
 * \param[in]  pmeTarget                  The user's choice for mdrun -pme for where to assign
 *                                        long-ranged PME nonbonded interaction tasks.
 * \param[in]  bondedTarget               The user's choice for mdrun -bonded for where to assign tasks.
 * \param[in]  updateTarget               The user's choice for mdrun -update for where to assign tasks.
 * \param[in]  useGpuForNonbonded         Whether GPUs will be used for nonbonded interactions.
 * \param[in]  useGpuForPme               Whether GPUs will be used for PME interactions.
 * \param[in]  rankHasPpTask              Whether this rank has a PP task
 * \param[in]  rankHasPmeTask             Whether this rank has a PME task
 */
std::vector<GpuTask> findGpuTasksOnThisRank(bool       haveGpusOnThisPhysicalNode,
                                            TaskTarget nonbondedTarget,
                                            TaskTarget pmeTarget,
                                            TaskTarget bondedTarget,
                                            TaskTarget updateTarget,
                                            bool       useGpuForNonbonded,
                                            bool       useGpuForPme,
                                            bool       rankHasPpTask,
                                            bool       rankHasPmeTask);

/*! \brief Returns container of all tasks on all ranks of this node
 * that are eligible for GPU execution.
 *
 * Perform all necessary communication for preparing for task
 * assignment. Separating this aspect makes it possible to unit test
 * the logic of task assignment. */
GpuTasksOnRanks findAllGpuTasksOnThisNode(ArrayRef<const GpuTask>         gpuTasksOnThisRank,
                                          const PhysicalNodeCommunicator& physicalNodeComm);

} // namespace gmx

#endif
