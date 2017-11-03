/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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
/*! \defgroup module_taskassignment Assigning simulation tasks to hardware (taskassignment)
 * \ingroup group_mdrun
 * \brief Provides code that manages assignment of simulation tasks to hardware.
 */
/*! \libinternal
 * \file
 * \brief Declares high-level functionality for managing assigning
 * tasks on ranks of a node to hardware on that node, and the factory
 * function to build the correct flavours of gmx::INodeTaskAssigner
 * required to implement the user's requirements.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_taskassignment
 * \inlibraryapi
 */
#ifndef GMX_TASKASSIGNMENT_TASKASSIGNMENT_H
#define GMX_TASKASSIGNMENT_TASKASSIGNMENT_H

#include <vector>

struct gmx_hw_info_t;
struct t_commrec;

namespace gmx
{

class MDLogger;

/*! \brief Types of compute tasks that can be run on a GPU.
 *
 * These names refer to existing practice in GROMACS, which is not
 * strictly accurate. */
enum class GpuTask : int
{
    //! Short-ranged interactions.
    Nonbonded,
    //! Long-ranged interactions.
    Pme
};

/*! \libinternal
 * \brief Specifies the GPU deviceID_ available for task_ to use. */
struct GpuTaskMapping
{
    //! The type of this GPU task.
    GpuTask task_;
    //! Device ID on this node to which this GPU task is mapped.
    int     deviceId_;
};

//! Container of GPU tasks on a rank, specifying the task type and GPU device ID, e.g. potentially ready for consumption by the modules on that rank.
using GpuTaskAssignment = std::vector <GpuTaskMapping>;
//! Container of compute tasks suitable to run on a GPU e.g. on each rank of a node.
using GpuTasksOnRanks = std::vector< std::vector<GpuTask> >;
//! Container of RankGpuTaskAssignments e.g. for all ranks on a node.
using GpuTaskAssignments = std::vector<GpuTaskAssignment>;

/*! \brief Coordinate the final stages of task assignment and
 * reporting, and return the assignment for this rank.
 *
 * Communicates between ranks on a node to coordinate task assignment
 * between them onto available hardware, e.g. accelerators.
 *
 * Releases the taskAssigner once its work is complete.
 *
 * \param[in]  gpuIdsToUse                The compatible GPUs that the user permitted us to use.
 * \param[in]  userGpuTaskAssignment      The user-specified assignment of GPU tasks to device IDs.
 * \param[in]  hardwareInfo               The detected hardware
 * \param[in]  mdlog                      Logging object to write to.
 * \param[in]  cr                         Communication object.
 * \param[in]  gpuTasksOnThisRank         Information about what GPU tasks
 *                                        exist on this rank.
 *
 * \returns  A GPU task assignment for this rank.
 *
 * \throws   std::bad_alloc          If out of memory.
 *           InconsistentInputError  If user and/or detected inputs are inconsistent.
 */
GpuTaskAssignments::value_type
runTaskAssignment(const std::vector<int>     &gpuIdsToUse,
                  const std::vector<int>     &userGpuTaskAssignment,
                  const gmx_hw_info_t        &hardwareInfo,
                  const MDLogger             &mdlog,
                  const t_commrec            *cr,
                  const std::vector<GpuTask> &gpuTasksOnThisRank);

//! Function for whether the task of \c mapping has value \c TaskType.
template<GpuTask TaskType>
bool hasTaskType(const GpuTaskMapping &mapping)
{
    return mapping.task_ == TaskType;
}

} // namespace

#endif
