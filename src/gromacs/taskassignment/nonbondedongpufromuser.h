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
 * Declares gmx::NonbondedOnGpuFromUser class
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_taskassignment
 */
#ifndef GMX_TASKASSIGNMENT_NONBONDEDONGPUFROMUSER_H
#define GMX_TASKASSIGNMENT_NONBONDEDONGPUFROMUSER_H

#include <string>
#include <vector>

#include "gromacs/hardware/hw_info.h"
#include "gromacs/taskassignment/inodetaskassigner.h"
#include "gromacs/utility/basedefinitions.h"

#include "nonbondedongpu-impl.h"

struct gmx_hw_info_t;

namespace gmx
{

class MDLogger;

/*! \internal
 * \brief Handles assigning tasks to run on compatible GPUs designated by the user.
 *
 * Such as assignment occurs for a Verlet-scheme tpr when compatible
 * GPUs are detected and indexed in e.g. mdrun -gpu_id . */
class NonbondedOnGpuFromUser : public INodeTaskAssigner
{
    //! User-supplied GPU devices ID for nonbondeds.
    std::string          gpuIdString_;
    //! Information detected about hardware, including GPU compatibility.
    const gmx_hw_info_t  hardwareInfo_;
    public:
        /*! \brief Constructor.
         *
         * \param[in]   gpuIdString      The user-supplied string mapping the tasks on ranks
         *                               of this node to the devices on this node.
         * \param[in]   hardwareInfo     The detected hardware
         * \throws      std::bad_alloc   If out of memory. */
        NonbondedOnGpuFromUser(const std::string   &gpuIdString,
                               const gmx_hw_info_t &hardwareInfo);
        //! \copydoc INodeTaskAssigner::chooseWhetherToUseGpus(const bool, const bool)
        void
        chooseWhetherToUseGpus(const bool usingVerletScheme,
                               const bool gpuAccelerationIsUseful) override;
        //! \copydoc INodeTaskAssigner::areNonbondedOnGpu()
        bool
        areNonbondedOnGpu() const override;
        //! \copydoc INodeTaskAssigner::assignGpuTasksToDeviceIds
        GpuTaskAssignments
        assignGpuTasksToDeviceIds(const GpuTasksOnRanks &gpuTasksOnRanks) const override;
        //! \copydoc INodeTaskAssigner::didUserAssignGpus()
        bool
        didUserAssignGpus() const override;
        //! \copydoc INodeTaskAssigner::logPerformanceHints(const MDLogger &, const GpuTaskAssignments &)
        void
        logPerformanceHints(const MDLogger           &mdlog,
                            const GpuTaskAssignments &gpuTaskAssignments) const override;
};

} // namespace

#endif
