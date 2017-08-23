/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015,2016,2017, by the GROMACS development team, led by
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
 * \brief Defines functionality for deciding whether tasks will run on GPUs.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_taskassignment
 */

#include "gmxpre.h"

#include "decidegpuusage.h"

#include "config.h"

#include <stdlib.h>
#include <string.h>

#include <algorithm>

#include "gromacs/hardware/cpuinfo.h"
#include "gromacs/hardware/detecthardware.h"
#include "gromacs/hardware/hardwaretopology.h"
#include "gromacs/hardware/hw_info.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdlib/nb_verlet.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/taskassignment/taskassignment.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/baseversion.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/stringutil.h"


namespace gmx
{

namespace
{

//! Record where a compute task is targetted.
enum class TaskTarget : int
{
    Auto,
    Cpu,
    Gpu
};

//! Make a TaskTarget from an mdrun argument string.
static TaskTarget findTaskTarget(const char *optionString)
{
    if (strncmp(optionString, "auto", 3) == 0)
    {
        return TaskTarget::Auto;
    }
    else if (strncmp(optionString, "cpu", 3) == 0)
    {
        return TaskTarget::Cpu;
    }
    else if (strncmp(optionString, "gpu", 3) == 0)
    {
        return TaskTarget::Gpu;
    }
    else
    {
        GMX_ASSERT(false, "Option string should have been checked for sanity already");
    }
}

}   // namespace

bool
decideWhetherToUseGpusForNonbondedWithThreadMpi(const char               *nbOptionString,
                                                const std::vector<int>   &gpuIds,
                                                const EmulateGpuNonbonded emulateGpuNonbonded,
                                                const gmx_hw_info_t      &hardwareInfo,
                                                const bool                usingVerletScheme,
                                                const bool                nonbondedOnGpuIsUseful,
                                                const int                 numRanksPerSimulation)
{
    auto nonbondedTarget = findTaskTarget(nbOptionString);

    // First, exclude all cases where we can't run NB on GPUs.

    if (nonbondedTarget == TaskTarget::Cpu ||
        (emulateGpuNonbonded == EmulateGpuNonbonded::Yes) ||
        !usingVerletScheme ||
        !nonbondedOnGpuIsUseful)
    {
        // If the user required NB on GPUs, we issue an error later.
        return false;
    }

    // We now know that NB on GPUs makes sense, if we have any.

    if (!gpuIds.empty())
    {
        // Force the user to be completely explicit when they have
        // supplied GPU IDs.
        if (nonbondedTarget == TaskTarget::Auto ||
            numRanksPerSimulation < 1)
        {
            GMX_THROW(InconsistentInputError
                          ("When you run mdrun -gpu_id, all of -nb, -pme and -ntmpi must be set to "
                          "non-default values, so that the device IDs can be interpreted correctly."));
        }
        return true;
    }

    auto haveGpus = !hardwareInfo.compatibleGpus.empty();

    // If we have GPUs, we run nonbonded tasks on them.
    return haveGpus;
}

// TODO should this be returning a PmeRunMode?
bool
decideWhetherToUseGpusForPmeWithThreadMpi(bool                    useGpuForNonbonded,
                                          const char             *pmeOptionString,
                                          const char             *pmeFftOptionString,
                                          const std::vector<int> &gpuIds,
                                          const gmx_hw_info_t    &hardwareInfo,
                                          const bool              usingAnyPme,
                                          const int               numRanksPerSimulation)
{
    auto pmeTarget = findTaskTarget(pmeOptionString);
    // TODO find a sensible home and behaviour for this
    auto pmeFftTarget    = findTaskTarget(pmeFftOptionString);
    GMX_UNUSED_VALUE(pmeFftTarget);

    // First, exclude all cases where we can't run PME on GPUs.

    if ((pmeTarget == TaskTarget::Cpu) ||
        !useGpuForNonbonded ||
        !usingAnyPme)
    {
        // Thus PME can't run on a GPU either. If the user required
        // that, we issue an error later.
        return false;
    }

    // We now know that PME on GPUs might make sense, if we have any.

    if (!gpuIds.empty())
    {
        // Follow the user's choice of GPU IDs, if we can. Checking
        // that their IDs are for compatible GPUs comes later.

        // Force the user to be completely explicit when they have
        // supplied GPU IDs.
        if (pmeTarget == TaskTarget::Auto ||
            numRanksPerSimulation < 1)
        {
            GMX_THROW(InconsistentInputError
                          ("When you run mdrun -gpu_id, all of -nb, -pme and -ntmpi must be set to "
                          "non-default values, so that the device IDs can be interpreted correctly."));
        }

        // PME on GPUs is only supported in a single case
        if (pmeTarget == TaskTarget::Gpu)
        {
            if (numRanksPerSimulation > 1)
            {
                GMX_THROW(InconsistentInputError
                              ("When you run mdrun -nb pme -gpu_id, you must supply a PME .tpr file and -ntmpi 1"));
            }
            return true;
        }

        // pmeTarget == TaskTarget::Auto
        return numRanksPerSimulation == 1;
    }

    if (numRanksPerSimulation < 1)
    {
        // Full automated mode for thread-MPI (the default).  PME can
        // run well on a GPU shared with NB, and we permit mdrun to
        // default to it in only that case.
        return (hardwareInfo.compatibleGpus.size() == 1);
    }

    if (numRanksPerSimulation == 1 && hardwareInfo.compatibleGpus.size() >= 1)
    {
        // PME can run well on a GPU shared with NB, and we permit
        // mdrun to default to that.
        return true;
    }

    // Not enough support for PME on GPUs for anything else
    return false;
}

// TODO should this function be split like the above two?
// TODO should this be returning a PmeRunMode?
std::pair<bool, bool> decideWhetherToUseGpus(const char                *nbOptionString,
                                             const char                *pmeOptionString,
                                             const std::vector<int>    &gpuIds,
                                             const EmulateGpuNonbonded  emulateGpuNonbonded,
                                             const gmx_hw_info_t       &hardwareInfo,
                                             const bool                 usingVerletScheme,
                                             const bool                 nonbondedOnGpuIsUseful,
                                             const int                  numRanksPerSimulation)
{
    auto nonbondedTarget = findTaskTarget(nbOptionString);
    auto pmeTarget       = findTaskTarget(pmeOptionString);

    if (nonbondedTarget == TaskTarget::Cpu)
    {
        if (pmeTarget == TaskTarget::Gpu)
        {
            GMX_THROW(NotImplementedError
                          ("The combination of nonbonded interactions on the CPU, and PME "
                          "on the GPU is not supported."));
        }
        if (!gpuIds.empty())
        {
            GMX_THROW(InconsistentInputError
                          ("GPU IDs were specified, but both nonbonded interactions were "
                          "assigned to the CPU. Make no more than one of these choices."));
        }

        return std::make_pair(false, false);
    }

    // TODO refactor all these TaskTarget::Gpu checks into one place?
    // e.g. use a subfunction that handles only the cases where
    // TaskTargets are not Cpu?
    if (emulateGpuNonbonded == EmulateGpuNonbonded::Yes)
    {
        if (nonbondedTarget == TaskTarget::Gpu)
        {
            GMX_THROW(InconsistentInputError
                          ("Nonbonded interactions on the GPU were required, which is inconsistent "
                          "with choosing emulation. Make no more than one of these choices."));
        }
        if (pmeTarget == TaskTarget::Gpu)
        {
            GMX_THROW(NotImplementedError
                          ("PME interactions on the GPU were required, which is not supported "
                          "with emulation of nonbonded interactions. Either permit PME to run "
                          "on the CPU, or do not use emulation."));
        }
        if (!gpuIds.empty())
        {
            auto message = "GPU IDs were specified, as was GPU emulation. Make no more than one of these choices.";
            GMX_THROW(InconsistentInputError(message));
        }

        return std::make_pair(false, false);
    }

    if (!usingVerletScheme)
    {
        if (nonbondedTarget == TaskTarget::Gpu)
        {
            GMX_THROW(InconsistentInputError
                          ("Nonbonded interactions on the GPU were required, which requires using "
                          "the Verlet scheme. Either use the Verlet scheme, or do not require using GPUs."));
        }
        if (pmeTarget == TaskTarget::Gpu)
        {
            GMX_THROW(NotImplementedError
                          ("PME interactions on the GPU were required, which requires using "
                          "the Verlet scheme. Either use the Verlet scheme, or do not require using GPUs."));
        }

        return std::make_pair(false, false);
    }

    if (!nonbondedOnGpuIsUseful)
    {
        if (nonbondedTarget == TaskTarget::Gpu)
        {
            GMX_THROW(InconsistentInputError
                          ("Nonbonded interactions on the GPU were required, but this would not be "
                          "useful. Probably you should not require using GPUs."));
        }
        if (pmeTarget == TaskTarget::Gpu)
        {
            GMX_THROW(NotImplementedError
                          ("PME interactions on the GPU were required, but that requires that "
                          "nonbonded interactions also run on the GPU, and the latter would not be "
                          "useful. Probably you should not require using GPUs."));
        }

        return std::make_pair(false, false);
    }

    // Note that we can't (yet) give an error if compatible GPUs are
    // missing, e.g. because that is fine on a node with only PME
    // ranks, which ignore -nb gpu and/or the GPU IDs if they are only
    // for nonbonded tasks. So we can't handle such aspects in this
    // function.
    if (!gpuIds.empty())
    {
        // TODO This means that mpirun -np 2 gmx_mpi mdrun -gpu_id 01
        // continues to work fine. And now mpirun -np 1 gmx_mpi mdrun
        // -gpu_id 1 would be rejected. Is this necessary and OK?
        if (numRanksPerSimulation == 1 &&
            (nonbondedTarget == TaskTarget::Auto ||
             pmeTarget == TaskTarget::Auto))
        {
            GMX_THROW(InconsistentInputError
                          ("When you set mdrun -gpu_id, the number of ranks and both of -nb and -pme must be set to "
                          "non-default values, so that the device IDs can be interpreted correctly."));
        }

        return std::make_pair(true, pmeTarget == TaskTarget::Gpu);
    }

    // If we get here, then the user permitted some degree of
    // automated assignment.
    auto haveGpus = !hardwareInfo.compatibleGpus.empty();
    if (numRanksPerSimulation == 1)
    {
        // This is a special case where we know that any tasks that
        // can run must do so on this rank, so we know already whether
        // we can honour -nb gpu, etc.
        if (!haveGpus)
        {
            // If there are no GPUs detected, then no tasks can run on them.
            if (nonbondedTarget == TaskTarget::Gpu)
            {
                GMX_THROW(InconsistentInputError
                              ("Nonbonded interactions on the GPU were required, but no GPUs were "
                              "detected."));
            }
            if (pmeTarget == TaskTarget::Gpu)
            {
                GMX_THROW(NotImplementedError
                              ("PME interactions on the GPU were required, but no GPUs were "
                              "detected."));
            }
            return std::make_pair(false, false);
        }

        // PME can run well on a single GPU shared with NB when
        // there is one rank, so we permit mdrun to default to
        // that.
        return std::make_pair(true, true);
    }

    if (pmeTarget == TaskTarget::Gpu)
    {
        GMX_THROW(NotImplementedError
                      ("PME tasks were required to run on GPUs, but that is not implemented with "
                      "more than one rank. Use a single rank, or permit PME tasks to be assigned "
                      "to the CPU."));
    }

    // If there are no GPUs detected, then no tasks can run on
    // them. If we have at least one GPU available, NB will run on
    // GPUs. We do not permit the default to choose to run a
    // separate PME rank (because it is unclear whether that is
    // useful), and DD with PME on GPUs is not implemented, so PME
    // will not run on GPUs at all.
    return std::make_pair(haveGpus, false);
}

} // namespace
