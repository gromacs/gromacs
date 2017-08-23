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
 * Declares gmx::NonbondedOnGpuAuto class
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_taskassignment
 */
#ifndef GMX_TASKASSIGNMENT_NONBONDEDONGPUAUTO_H
#define GMX_TASKASSIGNMENT_NONBONDEDONGPUAUTO_H

#include <vector>

#include "gromacs/hardware/hw_info.h"
#include "gromacs/taskassignment/inodetaskassigner.h"
#include "gromacs/utility/basedefinitions.h"

namespace gmx
{

class MDLogger;

/*! \internal
 * \brief Handles automatically assigning tasks to run on compatible GPUs.
 *
 * Such as assignment occurs by default for a Verlet-scheme tpr when
 * compatible GPUs are detected. However, if the user required the use
 * of GPUs for short-ranged interactions, then the handling of cases
 * where no compatible GPUs exist must differ, and this is handled by
 * class NonbondedOnGpuFromUser. */
class NonbondedOnGpuAuto : public INodeTaskAssigner
{
    //! Tracks whether the user forced the nonbondeds to run on a GPU.
    bool                 forceNonbondedOnPhysicalGpu_;
    //! Information detected about hardware, including GPU compatibility.
    const gmx_hw_info_t  hardwareInfo_;
    //! Whether the nonbondeds will run on a GPU.
    bool                 nonbondedOnGpu_;
    public:
        /*! \brief Constructor.
         *
         * \param[in]   forceNonbondedOnPhysicalGpu  Tracks whether the user forced the nonbondeds to run on a GPU.
         * \param[in]   hardwareInfo                 The detected hardware
         * \throws      std::bad_alloc               If out of memory. */
        NonbondedOnGpuAuto(bool                 forceNonbondedOnPhysicalGpu,
                           const gmx_hw_info_t &hardwareInfo);
        //! \copydoc INodeTaskAssigner::decideWhetherToUseGpus(const bool, const bool)
        void
        decideWhetherToUseGpus(const bool usingVerletScheme,
                               const bool gpuAccelerationIsUseful) override;
        //! \copydoc INodeTaskAssigner::areNonbondedOnGpu()
        bool
        areNonbondedOnGpu() const override;
        //! \copydoc INodeTaskAssigner::assignGpuTasksToDeviceIds()
        GpuTaskAssignments
        assignGpuTasksToDeviceIds(const GpuTasksOnRanks &gpuTasksOnRanks) const override;
        //! \copydoc INodeTaskAssigner::didUserAssignGpus()
        bool
        didUserAssignGpus() const override;
        //! \copydoc INodeTaskAssigner::logPerformanceHints()
        void
        logPerformanceHints(const MDLogger           &mdlog,
                            const GpuTaskAssignments &gpuTaskAssignments) const override;
};

} // namespace

#endif
