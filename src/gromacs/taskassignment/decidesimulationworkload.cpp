/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019,2020,2021, by the GROMACS development team, led by
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
 * Defines routine for building simulation workload task description.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \ingroup module_taskassignment
 */
#include "gmxpre.h"

#include "decidesimulationworkload.h"

#include "gromacs/ewald/pme.h"
#include "gromacs/mdtypes/multipletimestepping.h"
#include "gromacs/taskassignment/decidegpuusage.h"
#include "gromacs/taskassignment/taskassignment.h"
#include "gromacs/utility/arrayref.h"

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
                                            bool       useGpuDirectHalo)
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
    simulationWorkload.useGpuPmeFft = (pmeRunMode == PmeRunMode::Mixed);
    simulationWorkload.useGpuBonded = useGpuForBonded;
    simulationWorkload.useGpuUpdate = useGpuForUpdate;
    simulationWorkload.useGpuXBufferOps =
            (devFlags.enableGpuBufferOps || useGpuForUpdate) && !inputrec.useMts;
    simulationWorkload.useGpuFBufferOps =
            (devFlags.enableGpuBufferOps || useGpuForUpdate) && !inputrec.useMts;
    if (simulationWorkload.useGpuXBufferOps || simulationWorkload.useGpuFBufferOps)
    {
        GMX_ASSERT(simulationWorkload.useGpuNonbonded,
                   "Can only offload X/F buffer ops if nonbonded computation is also offloaded");
    }
    simulationWorkload.havePpDomainDecomposition = havePpDomainDecomposition;
    simulationWorkload.useCpuHaloExchange        = havePpDomainDecomposition && !useGpuDirectHalo;
    simulationWorkload.useGpuHaloExchange        = useGpuDirectHalo;
    if (pmeRunMode == PmeRunMode::None)
    {
        GMX_RELEASE_ASSERT(!haveSeparatePmeRank, "Can not have separate PME rank(s) without PME.");
    }
    simulationWorkload.haveSeparatePmeRank = haveSeparatePmeRank;
    simulationWorkload.useGpuPmePpCommunication =
            haveSeparatePmeRank && devFlags.enableGpuPmePPComm && (pmeRunMode == PmeRunMode::GPU);
    simulationWorkload.useCpuPmePpCommunication =
            haveSeparatePmeRank && !simulationWorkload.useGpuPmePpCommunication;
    GMX_RELEASE_ASSERT(!(simulationWorkload.useGpuPmePpCommunication
                         && simulationWorkload.useCpuPmePpCommunication),
                       "Cannot do PME-PP communication on both CPU and GPU");
    simulationWorkload.useGpuDirectCommunication =
            devFlags.enableGpuHaloExchange || devFlags.enableGpuPmePPComm;
    simulationWorkload.haveEwaldSurfaceContribution = haveEwaldSurfaceContribution(inputrec);
    simulationWorkload.useMts                       = inputrec.useMts;

    return simulationWorkload;
}

} // namespace gmx
