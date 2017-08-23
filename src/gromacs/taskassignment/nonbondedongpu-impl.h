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
#ifndef GMX_TASKASSIGNMENT_NONBONDEDONGPU_IMPL_H
#define GMX_TASKASSIGNMENT_NONBONDEDONGPU_IMPL_H
/*! \internal \file
 * \brief
 * Declares helper functionality for task assignment on GPUs
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_taskassignment
 */

#include <string>
#include <vector>

#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/taskassignment/inodetaskassigner.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

//! Helper template function, works with both GpuTasksAssignment and GpuTasksOnRanks objects.
template <typename T>
size_t countGpuTasksOnThisNode(const std::vector<std::vector<T> > &gpuTasksOnRanksOfThisNode);

// Declare that these templates will be instantiated in a single
// translation unit.
extern template
size_t countGpuTasksOnThisNode<GpuTask>(const std::vector<std::vector<GpuTask> > &);
extern template
size_t countGpuTasksOnThisNode<GpuTaskMapping>(const std::vector<std::vector<GpuTaskMapping> > &);

/*! \brief Build data structure of types of GPU tasks on a rank,
 * together with the mapped GPU device IDs, for all GPU tasks on all
 * the ranks of this node.
 *
 * \param[in]   flatGpuDeviceIds           A flat vector of GPU device IDs whose
 *                                         order maps GPU tasks to GPU device IDs.
 * \param[in]   gpuTasksOnRanksOfThisNode  For each rank on this node, the set of tasks
 *                                         that are eligible to run on GPUs.
 */
GpuTaskAssignments
buildTaskAssignment(const GpuTasksOnRanks  &gpuTasksOnRanksOfThisNode,
                    const std::vector<int> &flatGpuDeviceIds);

/*! \brief Implement simple round-robin assignment of all GPU tasks on
 * this node to ids of GPUs on this node, e.g. for mdrun -nb auto.
 *
 * Later, if we identify desirable assignments that e.g. reflect the
 * relative strengths of GPUs or cost of GPU tasks then this code
 * can become more complex.
 *
 * If there are limitations on whether ranks can share GPUs,
 * they could be implemented here.
 *
 * \throws InconsistentInputError when there are more GPU tasks than
 * compatible GPUs, and the former is not a multiple of the latter. */
std::vector<int>
makeFlatGpuTaskAssignment(size_t                  numGpuTasksOnThisNode,
                          const std::vector<int> &compatibleGpus);

} // namespace

#endif
