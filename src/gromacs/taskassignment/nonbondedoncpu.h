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
 * Declares gmx::NonbondedOnCpu class
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_taskassignment
 */
#ifndef GMX_TASKASSIGNMENT_NONBONDEDONCPU_H
#define GMX_TASKASSIGNMENT_NONBONDEDONCPU_H

#include <vector>

#include "gromacs/taskassignment/inodetaskassigner.h"
#include "gromacs/utility/basedefinitions.h"

namespace gmx
{

class MDLogger;

/*! \internal
 * \brief Handles assigning nonbonded tasks to run solely on CPUs.
 *
 * Such an assignment might occur when the user required such, e.g with
 * mdrun -nb cpu, or a request for emulation.
 *
 * Currently the implementation is trivial, but later we expect to
 * have to grapple also with CPU-side assignment issues, e.g. for
 * locality.
 *
 * \todo Later, this assigner might also handle preparing the implementation
 * of emulation, rather than the Verlet-scheme checking GMX_EMULATE_GPU
 * independently. */
class NonbondedOnCpu : public INodeTaskAssigner
{
    public:
        NonbondedOnCpu();
        //! \copydoc INodeTaskAssigner::decideWhetherToUseGpus(const bool, const bool)
        void
        decideWhetherToUseGpus(const bool usingVerletScheme,
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
