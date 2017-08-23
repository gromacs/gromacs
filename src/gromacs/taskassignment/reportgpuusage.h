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
/*! \internal
 * \file
 * \brief Declares routine for reporting GPU usage.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_taskassignment
 */
#ifndef GMX_TASKASSIGNMENT_REPORTGPUUSAGE_H
#define GMX_TASKASSIGNMENT_REPORTGPUUSAGE_H

#include <cstdlib>

#include "gromacs/taskassignment/taskassignment.h"

namespace gmx
{

class MDLogger;

/*! \brief Log a report on how GPUs are being used on
 * the ranks of the physical node of rank 0 of the simulation.
 *
 * \todo It could be useful to report also whether any nodes differed,
 * and in what way.
 *
 * \param[in]  mdlog                               Logging object.
 * \param[in]  userSetGpuIds                       Whether the user selected the GPU ids
 * \param[in]  gpuTaskAssignmentOnRanksOfThisNode  The selected GPU IDs.
 * \param[in]  numGpuTasksOnThisNode               The number of GPU tasks on this node.
 * \param[in]  numPpRanks                          Number of PP ranks on this node
 * \param[in]  bPrintHostName                      Print the hostname in the usage information
 *
 * \throws     std::bad_alloc if out of memory */
void
reportGpuUsage(const MDLogger                &mdlog,
               bool                           userSetGpuIds,
               const GpuTaskAssignments      &gpuTaskAssignmentOnRanksOfThisNode,
               size_t                         numGpuTasksOnThisNode,
               size_t                         numPpRanks,
               bool                           bPrintHostName);


} // namespace

#endif
