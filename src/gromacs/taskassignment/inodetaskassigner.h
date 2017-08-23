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
/*! \libinternal \file
 * \brief Declares the gmx::INodeTaskAssigner interface that is needed
 * to implement the various flavours of task assignment.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_taskassignment
 * \inlibraryapi
 */
#ifndef GMX_TASKASSIGNMENT_INODETASKASSIGNER_H
#define GMX_TASKASSIGNMENT_INODETASKASSIGNER_H

#include <vector>

#include "gromacs/taskassignment/taskassignment.h"

namespace gmx
{

class MDLogger;

/*! \libinternal
 * \brief Interface for handling assigning GPU tasks to devices
 * according to user choices.
 *
 * Types that implement this interface permit us to separate the
 * handling of different user choices, in the context of hardware
 * found on each node, in ways that permit related functionality to be
 * consistent at different points when a simulation runner needs to
 * call it. Generally, such task assigners will have state only
 * related to user-supplied data necessary for that particular kind of
 * task assigner to function.
 *
 * Three implementations exist for respectively handling scenarios
 * where the user explicitly required assignment to CPUs, explictily
 * assigned specific GPU IDs, or allowed automatic assignment.
 *
 * The different implementations throw gmx::InconsistentInputError
 * when the user input and/or available hardware do not permit
 * a valid task assignment to be made. One method also handles
 * the logging of performance hints.
 * */
class INodeTaskAssigner
{
    public:
        /*! \brief Given the information from the inputrec, choose whether GPUs will be used for non-bonded interactions.
         *
         * \param[in]  usingVerletScheme        Whether the nonbondeds are using the Verlet scheme.
         * \param[in]  gpuAccelerationIsUseful  Whether GPU acceleration is useful for this calculation.
         * \throws     InconsistentInputError   If the assigner's choice to not use GPUs is inconsistent with the user's requirement to do so. */
        virtual void
        chooseWhetherToUseGpus(const bool usingVerletScheme,
                               const bool gpuAccelerationIsUseful) = 0;
        //! Return whether nonbonded tasks will run on GPUs.
        virtual bool
        areNonbondedOnGpu() const = 0;
        /*! \brief Returns a task assignment for each rank on this node.
         *
         * \param[in]  gpuTasksOnRanks         Contains the descriptions of any tasks on ranks of this node eligible to be run on GPUs.
         * \returns    Task assignment of GPU device IDs for each task on each rank of this node.
         * \throws     std::bad_alloc          If out of memory
         *             InconsistentInputError  If the assigned GPUs are not valid. */
        virtual GpuTaskAssignments
        assignGpuTasksToDeviceIds(const GpuTasksOnRanks &gpuTasksOnRanks) const = 0;
        //! Return whether the user specified the GPUs for this task assignment.
        virtual bool
        didUserAssignGpus() const = 0;
        //! Logs to \c mdlog information that may help a user learn how to let mdrun make a task assignment that runs faster.
        virtual void
        logPerformanceHints(const MDLogger           &mdlog,
                            const GpuTaskAssignments &gpuTaskAssignments) const = 0;
        //! Virtual destructor required for interfaces with virtual functions.
        virtual ~INodeTaskAssigner() {};
};

} // namespace

#endif
