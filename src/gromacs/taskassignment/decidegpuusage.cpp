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
#include <string>

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

//! Helper variable to localise the text of an often repeated message.
const char * g_specifyEverythingFormatString =
    "When you use mdrun -gputasks, %s must be set to non-default "
    "values, so that the device IDs can be interpreted correctly."
#if GMX_GPU != GMX_GPU_NONE
    " If you simply want to restrict which GPUs are used, then it is "
    "better to use mdrun -gpu_id. Otherwise, setting the "
#  if GMX_GPU == GMX_GPU_CUDA
    "CUDA_VISIBLE_DEVICES"
#  elif GMX_GPU == GMX_GPU_OPENCL
    // Technically there is no portable way to do this offered by the
    // OpenCL standard, but the only current relevant case for GROMACS
    // is AMD OpenCL, which offers this variable.
    "GPU_DEVICE_ORDINAL"
#  else
#  error "Unreachable branch"
#  endif
    " environment variable in your bash profile or job "
    "script may be more convenient."
#endif
;

}   // namespace

bool
decideWhetherToUseGpusForNonbondedWithThreadMpi(const TaskTarget          nonbondedTarget,
                                                const std::vector<int>   &gpuIdsToUse,
                                                const std::vector<int>   &userGpuTaskAssignment,
                                                const EmulateGpuNonbonded emulateGpuNonbonded,
                                                const bool                usingVerletScheme,
                                                const bool                nonbondedOnGpuIsUseful,
                                                const int                 numRanksPerSimulation)
{
    // First, exclude all cases where we can't run NB on GPUs.
    if (nonbondedTarget == TaskTarget::Cpu ||
        emulateGpuNonbonded == EmulateGpuNonbonded::Yes ||
        !usingVerletScheme ||
        !nonbondedOnGpuIsUseful)
    {
        // If the user required NB on GPUs, we issue an error later.
        return false;
    }

    // We now know that NB on GPUs makes sense, if we have any.

    if (!userGpuTaskAssignment.empty())
    {
        // Specifying -gputasks requires specifying everything.
        if (nonbondedTarget == TaskTarget::Auto ||
            numRanksPerSimulation < 1)
        {
            GMX_THROW(InconsistentInputError(formatString(g_specifyEverythingFormatString, "-nb and -ntmpi")));
        }
        return true;
    }

    if (nonbondedTarget == TaskTarget::Gpu)
    {
        return true;
    }

    // Because this is thread-MPI, we already know about the GPUs that
    // all potential ranks can use, and can use that in a global
    // decision that will later be consistent.
    auto haveGpus = !gpuIdsToUse.empty();

    // If we get here, then the user permitted or required GPUs.
    return haveGpus;
}

bool
decideWhetherToUseGpusForPmeWithThreadMpi(const bool              useGpuForNonbonded,
                                          const TaskTarget        pmeTarget,
                                          const std::vector<int> &gpuIdsToUse,
                                          const std::vector<int> &userGpuTaskAssignment,
                                          const bool              canUseGpuForPme,
                                          const int               numRanksPerSimulation,
                                          const int               numPmeRanksPerSimulation)
{
    // First, exclude all cases where we can't run PME on GPUs.
    if ((pmeTarget == TaskTarget::Cpu) ||
        !useGpuForNonbonded ||
        !canUseGpuForPme)
    {
        // PME can't run on a GPU. If the user required that, we issue
        // an error later.
        return false;
    }

    // We now know that PME on GPUs might make sense, if we have any.

    if (!userGpuTaskAssignment.empty())
    {
        // Follow the user's choice of GPU task assignment, if we
        // can. Checking that their IDs are for compatible GPUs comes
        // later.

        // Specifying -gputasks requires specifying everything.
        if (pmeTarget == TaskTarget::Auto ||
            numRanksPerSimulation < 1)
        {
            GMX_THROW(InconsistentInputError(formatString(g_specifyEverythingFormatString, "all of -nb, -pme, and -ntmpi")));
        }

        // PME on GPUs is only supported in a single case
        if (pmeTarget == TaskTarget::Gpu)
        {
            if (((numRanksPerSimulation > 1) && (numPmeRanksPerSimulation == 0)) ||
                (numPmeRanksPerSimulation > 1))
            {
                GMX_THROW(InconsistentInputError
                              ("When you run mdrun -pme gpu -gputasks, you must supply a PME-enabled .tpr file and use a single PME rank."));
            }
            return true;
        }

        // pmeTarget == TaskTarget::Auto
        return numRanksPerSimulation == 1;
    }

    // Because this is thread-MPI, we already know about the GPUs that
    // all potential ranks can use, and can use that in a global
    // decision that will later be consistent.

    if (pmeTarget == TaskTarget::Gpu)
    {
        if (((numRanksPerSimulation > 1) && (numPmeRanksPerSimulation == 0)) ||
            (numPmeRanksPerSimulation > 1))
        {
            GMX_THROW(NotImplementedError
                          ("PME tasks were required to run on GPUs, but that is not implemented with "
                          "more than one PME rank. Use a single rank simulation, or a separate PME rank, "
                          "or permit PME tasks to be assigned to the CPU."));
        }
        return true;
    }

    if (numRanksPerSimulation == 1)
    {
        // PME can run well on a GPU shared with NB, and we permit
        // mdrun to default to try that.
        return gpuIdsToUse.size() >= 1;
    }

    if (numRanksPerSimulation < 1)
    {
        // Full automated mode for thread-MPI (the default). PME can
        // run well on a GPU shared with NB, and we permit mdrun to
        // default to it if there is only one GPU available.
        return (gpuIdsToUse.size() == 1);
    }

    // Not enough support for PME on GPUs for anything else
    return false;
}

bool decideWhetherToUseGpusForNonbonded(const TaskTarget           nonbondedTarget,
                                        const std::vector<int>    &userGpuTaskAssignment,
                                        const EmulateGpuNonbonded  emulateGpuNonbonded,
                                        const bool                 usingVerletScheme,
                                        const bool                 nonbondedOnGpuIsUseful,
                                        const bool                 gpusWereDetected)
{
    if (nonbondedTarget == TaskTarget::Cpu)
    {
        if (!userGpuTaskAssignment.empty())
        {
            GMX_THROW(InconsistentInputError
                          ("A GPU task assignment was specified, but nonbonded interactions were "
                          "assigned to the CPU. Make no more than one of these choices."));
        }

        return false;
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
        if (!userGpuTaskAssignment.empty())
        {
            GMX_THROW(InconsistentInputError
                          ("GPU ID usage was specified, as was GPU emulation. Make no more than one of these choices."));
        }

        return false;
    }

    if (!usingVerletScheme)
    {
        if (nonbondedTarget == TaskTarget::Gpu)
        {
            GMX_THROW(InconsistentInputError
                          ("Nonbonded interactions on the GPU were required, which requires using "
                          "the Verlet scheme. Either use the Verlet scheme, or do not require using GPUs."));
        }

        return false;
    }

    if (!nonbondedOnGpuIsUseful)
    {
        if (nonbondedTarget == TaskTarget::Gpu)
        {
            GMX_THROW(InconsistentInputError
                          ("Nonbonded interactions on the GPU were required, but this would not be "
                          "useful. Probably you should not require using GPUs."));
        }

        return false;
    }

    if (!userGpuTaskAssignment.empty())
    {
        // Specifying -gputasks requires specifying everything.
        if (nonbondedTarget == TaskTarget::Auto)
        {
            GMX_THROW(InconsistentInputError(formatString(g_specifyEverythingFormatString, "-nb and -ntmpi")));
        }

        return true;
    }

    if (nonbondedTarget == TaskTarget::Gpu)
    {
        // We still don't know whether it is an error if no GPUs are found
        // because we don't know the duty of this rank, yet. For example,
        // a node with only PME ranks and -pme cpu is OK if there are not
        // GPUs.
        return true;
    }

    // If we get here, then the user permitted GPUs, which we should
    // use for nonbonded interactions.
    return gpusWereDetected;
}

bool decideWhetherToUseGpusForPme(const bool              useGpuForNonbonded,
                                  const TaskTarget        pmeTarget,
                                  const std::vector<int> &userGpuTaskAssignment,
                                  const bool              canUseGpuForPme,
                                  const int               numRanksPerSimulation,
                                  const int               numPmeRanksPerSimulation,
                                  const bool              gpusWereDetected)
{
    if (pmeTarget == TaskTarget::Cpu)
    {
        return false;
    }

    if (!useGpuForNonbonded)
    {
        if (pmeTarget == TaskTarget::Gpu)
        {
            GMX_THROW(NotImplementedError
                          ("The PME on the GPU is only supported when nonbonded interactions run on GPUs also."));
        }
        return false;
    }

    if (!canUseGpuForPme)
    {
        if (pmeTarget == TaskTarget::Gpu)
        {
            // TODO Pass in the inputrec so we can give more help here?
            GMX_THROW(NotImplementedError
                          ("The input simulation did not use PME in a way that is supported on the GPU."));
        }
        return false;
    }

    if (pmeTarget == TaskTarget::Cpu)
    {
        if (!userGpuTaskAssignment.empty())
        {
            GMX_THROW(InconsistentInputError
                          ("A GPU task assignment was specified, but PME interactions were "
                          "assigned to the CPU. Make no more than one of these choices."));
        }

        return false;
    }

    if (!userGpuTaskAssignment.empty())
    {
        // Specifying -gputasks requires specifying everything.
        if (pmeTarget == TaskTarget::Auto)
        {
            GMX_THROW(InconsistentInputError(formatString(g_specifyEverythingFormatString, "all of -nb, -pme, and -ntmpi")));
        }

        return true;
    }

    // We still don't know whether it is an error if no GPUs are found
    // because we don't know the duty of this rank, yet. For example,
    // a node with only PME ranks and -pme cpu is OK if there are not
    // GPUs.

    if (pmeTarget == TaskTarget::Gpu)
    {
        if (((numRanksPerSimulation > 1) && (numPmeRanksPerSimulation == 0)) ||
            (numPmeRanksPerSimulation > 1))
        {
            GMX_THROW(NotImplementedError
                          ("PME tasks were required to run on GPUs, but that is not implemented with "
                          "more than one PME rank. Use a single rank simulation, or a separate PME rank, "
                          "or permit PME tasks to be assigned to the CPU."));
        }
        return true;
    }

    // If we get here, then the user permitted GPUs.
    if (numRanksPerSimulation == 1)
    {
        // PME can run well on a single GPU shared with NB when there
        // is one rank, so we permit mdrun to try that if we have
        // detected GPUs.
        return gpusWereDetected;
    }

    // Not enough support for PME on GPUs for anything else
    return false;
}

} // namespace
