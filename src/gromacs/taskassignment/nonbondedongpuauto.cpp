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
/*! \internal \file
 * \brief
 * Defines gmx::INodeTaskAssigner interface, its concrete classes, and factory function.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_taskassignment
 */
#include "gmxpre.h"

#include "nonbondedongpuauto.h"

#include <algorithm>
#include <string>
#include <vector>

#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/taskassignment/inodetaskassigner.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/sysinfo.h"

#include "nonbondedongpu-impl.h"

struct gmx_gpu_info_t;

namespace gmx
{

NonbondedOnGpuAuto::NonbondedOnGpuAuto(bool                 forceNonbondedOnPhysicalGpu,
                                       const gmx_hw_info_t &hardwareInfo)
    : forceNonbondedOnPhysicalGpu_(forceNonbondedOnPhysicalGpu), hardwareInfo_(hardwareInfo),
      nonbondedOnGpu_(false) {}

void
NonbondedOnGpuAuto::decideWhetherToUseGpus(const bool usingVerletScheme,
                                           const bool gpuAccelerationIsUseful)
{
    bool nonbondedOnGpu = usingVerletScheme && gpuAccelerationIsUseful;
    if (forceNonbondedOnPhysicalGpu_ && !nonbondedOnGpu)
    {
        GMX_THROW(InconsistentInputError
                      ("GPU acceleration required, but that is not supported with the given input "
                      "settings. Either do not request GPUs, or change your inputs to something supported."));
    }
    // Implements strong exception guarantee.
    nonbondedOnGpu_ = nonbondedOnGpu;
}

bool NonbondedOnGpuAuto::areNonbondedOnGpu() const
{
    return nonbondedOnGpu_;
}

GpuTaskAssignments
NonbondedOnGpuAuto::assignGpuTasksToDeviceIds(const GpuTasksOnRanks &gpuTasksOnRanks) const
{
    auto numGpuTasksOnThisNode = countGpuTasksOnThisNode(gpuTasksOnRanks);
    if (numGpuTasksOnThisNode == 0)
    {
        // No GPU tasks found on any rank on this node, so
        // return an empty assignment for each rank.
        return GpuTaskAssignments(gpuTasksOnRanks.size());
    }

    if (hardwareInfo_.compatibleGpus.empty())
    {
        if (forceNonbondedOnPhysicalGpu_)
        {
            // TODO The annotation of node name could perhaps be the
            // responsibility of the code that catches the exception,
            // e.g. a logger.
            char host[STRLEN];
            gmx_gethostname(host, STRLEN);

            auto message = formatString("Work on GPUs was requested on a rank on node %s, but no "
                                        "compatible GPUs were detected on that node. All nodes with such ranks need "
                                        "to have GPUs. If you intended to use GPU acceleration in a parallel run, you "
                                        "can either fix the affected nodes, avoid using the nodes that don't have "
                                        "GPUs or place ranks with only CPU work on these nodes.", host);
            GMX_THROW(InconsistentInputError(message));
        }
        else
        {
            GpuTaskAssignments gpuTaskAssignmentOnRanksOfThisNode(gpuTasksOnRanks.size());
            return gpuTaskAssignmentOnRanksOfThisNode;
        }
    }
    return buildTaskAssignment(gpuTasksOnRanks, makeFlatGpuTaskAssignment(numGpuTasksOnThisNode, hardwareInfo_.compatibleGpus));
}

bool
NonbondedOnGpuAuto::didUserAssignGpus() const
{
    return false;
}

void
NonbondedOnGpuAuto::logPerformanceHints(const MDLogger           & /*mdlog*/,
                                        const GpuTaskAssignments & /*gpuTaskAssignments*/) const
{
}

} // namespace
