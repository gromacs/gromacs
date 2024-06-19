/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
/*! \internal \file
 * \brief Declares utility functions to manage step, domain-lifetime, and run workload data structures.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \author Szilárd Páll <pall.szilard@gmail.com>
 * \ingroup module_taskassignment
 */
#include "gmxpre.h"

#include "gromacs/taskassignment/decidesimulationworkload.h"

#include "config.h"

#include <cstdint>

#include <bitset>
#include <memory>
#include <vector>

#include "gromacs/essentialdynamics/edsam.h"
#include "gromacs/ewald/pme.h"
#include "gromacs/listed_forces/listed_forces.h"
#include "gromacs/listed_forces/listed_forces_gpu.h"
#include "gromacs/mdlib/force_flags.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/iforceprovider.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/multipletimestepping.h"
#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/pulling/pull.h"
#include "gromacs/taskassignment/decidegpuusage.h"
#include "gromacs/taskassignment/taskassignment.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/gmxassert.h"

struct gmx_edsam;
struct pull_t;

namespace gmx
{

SimulationWorkload createSimulationWorkload(const t_inputrec& inputrec,
                                            const bool        disableNonbondedCalculation,
                                            const DevelopmentFeatureFlags& devFlags,
                                            bool       havePpDomainDecomposition,
                                            bool       haveSeparatePmeRank,
                                            bool       useGpuForNonbonded,
                                            PmeRunMode pmeRunMode,
                                            bool       useGpuForBonded,
                                            bool       useGpuForUpdate,
                                            bool       useGpuDirectHalo,
                                            bool       canUseDirectGpuComm,
                                            bool       useGpuPmeDecomposition)
{
    SimulationWorkload simulationWorkload;
    simulationWorkload.computeNonbonded = !disableNonbondedCalculation;
    simulationWorkload.computeNonbondedAtMtsLevel1 =
            simulationWorkload.computeNonbonded && inputrec.useMts
            && inputrec.mtsLevels.back().forceGroups[static_cast<int>(MtsForceGroups::Nonbonded)];
    simulationWorkload.computeMuTot    = inputrecNeedMutot(&inputrec);
    simulationWorkload.useCpuNonbonded = !useGpuForNonbonded;
    simulationWorkload.useGpuNonbonded = useGpuForNonbonded;
    simulationWorkload.useCpuPme       = (pmeRunMode == PmeRunMode::CPU);
    simulationWorkload.useGpuPme = (pmeRunMode == PmeRunMode::GPU || pmeRunMode == PmeRunMode::Mixed);
    simulationWorkload.useGpuPmeFft              = (pmeRunMode == PmeRunMode::GPU);
    simulationWorkload.useGpuBonded              = useGpuForBonded;
    simulationWorkload.useGpuUpdate              = useGpuForUpdate;
    simulationWorkload.havePpDomainDecomposition = havePpDomainDecomposition;
    simulationWorkload.useCpuHaloExchange        = havePpDomainDecomposition && !useGpuDirectHalo;
    simulationWorkload.useGpuHaloExchange        = useGpuDirectHalo;
    if (pmeRunMode == PmeRunMode::None)
    {
        GMX_RELEASE_ASSERT(!haveSeparatePmeRank, "Can not have separate PME rank(s) without PME.");
    }
    simulationWorkload.haveSeparatePmeRank = haveSeparatePmeRank;
    simulationWorkload.useGpuPmePpCommunication =
            haveSeparatePmeRank && canUseDirectGpuComm
            && (pmeRunMode == PmeRunMode::GPU || pmeRunMode == PmeRunMode::Mixed);
    simulationWorkload.useCpuPmePpCommunication =
            haveSeparatePmeRank && !simulationWorkload.useGpuPmePpCommunication;
    GMX_RELEASE_ASSERT(!(simulationWorkload.useGpuPmePpCommunication
                         && simulationWorkload.useCpuPmePpCommunication),
                       "Cannot do PME-PP communication on both CPU and GPU");
    simulationWorkload.useGpuDirectCommunication =
            simulationWorkload.useGpuHaloExchange || simulationWorkload.useGpuPmePpCommunication;
    simulationWorkload.useGpuPmeDecomposition       = useGpuPmeDecomposition;
    simulationWorkload.haveEwaldSurfaceContribution = haveEwaldSurfaceContribution(inputrec);
    simulationWorkload.useMts                       = inputrec.useMts;
    const bool featuresRequireGpuBufferOps = useGpuForUpdate || simulationWorkload.useGpuDirectCommunication;
    simulationWorkload.useGpuXBufferOpsWhenAllowed =
            (devFlags.enableGpuBufferOps || featuresRequireGpuBufferOps) && !inputrec.useMts;
    simulationWorkload.useGpuFBufferOpsWhenAllowed =
            (devFlags.enableGpuBufferOps || featuresRequireGpuBufferOps) && !inputrec.useMts;
    if (simulationWorkload.useGpuXBufferOpsWhenAllowed || simulationWorkload.useGpuFBufferOpsWhenAllowed)
    {
        GMX_ASSERT(simulationWorkload.useGpuNonbonded,
                   "Can only offload X/F buffer ops if nonbonded computation is also offloaded");
    }
    simulationWorkload.useMdGpuGraph =
            devFlags.enableCudaGraphs && useGpuForUpdate
            && (simulationWorkload.haveSeparatePmeRank ? simulationWorkload.useGpuPmePpCommunication : true)
            && (havePpDomainDecomposition ? simulationWorkload.useGpuHaloExchange : true)
            && (havePpDomainDecomposition ? (GMX_THREAD_MPI > 0) : true);

    simulationWorkload.useNvshmem = devFlags.enableNvshmem && simulationWorkload.useGpuPmePpCommunication;
    return simulationWorkload;
}


/*! \brief Return true if there are special forces computed.
 *
 * The conditionals exactly correspond to those in sim_util.cpp:computeSpecialForces().
 */
static bool haveSpecialForces(const t_inputrec&          inputrec,
                              const gmx::ForceProviders& forceProviders,
                              const pull_t*              pull_work,
                              const gmx_edsam*           ed)
{

    return ((forceProviders.hasForceProvider()) ||                 // forceProviders
            (inputrec.bPull && pull_have_potential(*pull_work)) || // pull
            inputrec.bRot ||                                       // enforced rotation
            (ed != nullptr) ||                                     // flooding
            (inputrec.bIMD));                                      // IMD
}

DomainLifetimeWorkload setupDomainLifetimeWorkload(const t_inputrec&         inputrec,
                                                   const t_forcerec&         fr,
                                                   const pull_t*             pull_work,
                                                   const gmx_edsam*          ed,
                                                   const t_mdatoms&          mdatoms,
                                                   const SimulationWorkload& simulationWork)
{
    DomainLifetimeWorkload domainWork;
    // Note that haveSpecialForces is constant over the whole run
    domainWork.haveSpecialForces = haveSpecialForces(inputrec, *fr.forceProviders, pull_work, ed);
    domainWork.haveCpuListedForceWork = false;
    domainWork.haveCpuBondedWork      = false;
    for (const auto& listedForces : fr.listedForces)
    {
        if (listedForces.haveCpuListedForces(*fr.fcdata))
        {
            domainWork.haveCpuListedForceWork = true;
        }
        if (listedForces.haveCpuBondeds())
        {
            domainWork.haveCpuBondedWork = true;
        }
    }
    domainWork.haveGpuBondedWork =
            ((fr.listedForcesGpu != nullptr) && fr.listedForcesGpu->haveInteractions());
    // Note that haveFreeEnergyWork is constant over the whole run
    domainWork.haveFreeEnergyWork =
            (fr.efep != FreeEnergyPerturbationType::No && mdatoms.nPerturbed != 0);
    // We assume we have local force work if there are CPU
    // force tasks including PME or nonbondeds.
    domainWork.haveCpuLocalForceWork =
            domainWork.haveSpecialForces || domainWork.haveCpuListedForceWork
            || domainWork.haveFreeEnergyWork || simulationWork.useCpuNonbonded || simulationWork.useCpuPme
            || simulationWork.haveEwaldSurfaceContribution || inputrec.nwall > 0;
    domainWork.haveCpuNonLocalForceWork = domainWork.haveCpuBondedWork || domainWork.haveFreeEnergyWork;
    domainWork.haveLocalForceContribInCpuBuffer =
            domainWork.haveCpuLocalForceWork || simulationWork.havePpDomainDecomposition;

    return domainWork;
}

/*! \brief Set up force flag struct from the force bitmask.
 *
 * \param[in]      legacyFlags          Force bitmask flags used to construct the new flags
 * \param[in]      mtsLevels            The multiple time-stepping levels, either empty or 2 levels
 * \param[in]      step                 The current MD step
 * \param[in]      domainWork           Domain lifetime workload description.
 * \param[in]      simulationWork       Simulation workload description.
 *
 * \returns New Stepworkload description.
 */
StepWorkload setupStepWorkload(const int                     legacyFlags,
                               ArrayRef<const gmx::MtsLevel> mtsLevels,
                               const int64_t                 step,
                               const DomainLifetimeWorkload& domainWork,
                               const SimulationWorkload&     simulationWork)
{
    GMX_ASSERT(mtsLevels.empty() || mtsLevels.size() == 2, "Expect 0 or 2 MTS levels");
    const bool computeSlowForces = (mtsLevels.empty() || step % mtsLevels[1].stepFactor == 0);

    StepWorkload flags;
    flags.stateChanged                  = ((legacyFlags & GMX_FORCE_STATECHANGED) != 0);
    flags.haveDynamicBox                = ((legacyFlags & GMX_FORCE_DYNAMICBOX) != 0);
    flags.doNeighborSearch              = ((legacyFlags & GMX_FORCE_NS) != 0);
    flags.computeSlowForces             = computeSlowForces;
    flags.computeVirial                 = ((legacyFlags & GMX_FORCE_VIRIAL) != 0);
    flags.computeEnergy                 = ((legacyFlags & GMX_FORCE_ENERGY) != 0);
    flags.computeForces                 = ((legacyFlags & GMX_FORCE_FORCES) != 0);
    flags.useOnlyMtsCombinedForceBuffer = ((legacyFlags & GMX_FORCE_DO_NOT_NEED_NORMAL_FORCE) != 0);
    flags.computeListedForces           = ((legacyFlags & GMX_FORCE_LISTED) != 0);
    flags.computeNonbondedForces =
            ((legacyFlags & GMX_FORCE_NONBONDED) != 0) && simulationWork.computeNonbonded
            && !(simulationWork.computeNonbondedAtMtsLevel1 && !computeSlowForces);
    flags.computeDhdl = ((legacyFlags & GMX_FORCE_DHDL) != 0);

    if (simulationWork.useGpuXBufferOpsWhenAllowed || simulationWork.useGpuFBufferOpsWhenAllowed)
    {
        GMX_ASSERT(simulationWork.useGpuNonbonded,
                   "Can only offload buffer ops if nonbonded computation is also offloaded");
    }
    flags.useGpuXBufferOps = simulationWork.useGpuXBufferOpsWhenAllowed && !flags.doNeighborSearch;
    // on virial steps the CPU reduction path is taken
    flags.useGpuFBufferOps = simulationWork.useGpuFBufferOpsWhenAllowed && !flags.computeVirial;
    flags.useGpuPmeFReduction =
            flags.computeSlowForces && flags.useGpuFBufferOps
            && (simulationWork.haveGpuPmeOnPpRank() || simulationWork.useGpuPmePpCommunication);
    flags.useGpuXHalo              = simulationWork.useGpuHaloExchange && !flags.doNeighborSearch;
    flags.useGpuFHalo              = simulationWork.useGpuHaloExchange && flags.useGpuFBufferOps;
    flags.haveGpuPmeOnThisRank     = simulationWork.haveGpuPmeOnPpRank() && flags.computeSlowForces;
    flags.computePmeOnSeparateRank = simulationWork.haveSeparatePmeRank && flags.computeSlowForces;
    flags.combineMtsForcesBeforeHaloExchange =
            (flags.computeForces && simulationWork.useMts && flags.computeSlowForces
             && flags.useOnlyMtsCombinedForceBuffer
             && !(flags.computeVirial || simulationWork.useGpuNonbonded || flags.haveGpuPmeOnThisRank));
    // On NS steps, the buffer is cleared in stateGpu->reinit, no need to clear it twice.
    flags.clearGpuFBufferEarly =
            flags.useGpuFHalo && !domainWork.haveCpuLocalForceWork && !flags.doNeighborSearch;

    return flags;
}

} // namespace gmx
