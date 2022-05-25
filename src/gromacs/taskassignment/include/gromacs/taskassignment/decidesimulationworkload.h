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
/*! \internal
 * \file
 * \brief Declares routine for deciding simulation workload based on GPU tasks.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \ingroup module_taskassignment
 */
#ifndef GMX_TASKASSIGNMENT_DECIDESIMULATIONWORKLOAD_H
#define GMX_TASKASSIGNMENT_DECIDESIMULATIONWORKLOAD_H

#include <vector>

#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/simulation_workload.h"

enum class PmeRunMode;

namespace gmx
{

struct DevelopmentFeatureFlags;

/*! \brief
 * Build datastructure that contains decisions whether to run different workload
 * task on GPUs.
 *
 * \param[in] inputrec           The input record
 * \param[in] disableNonbondedCalculation  Disable calculation of nonbonded forces
 * \param[in] devFlags           The development feature flags
 * \param[in] havePpDomainDecomposition Whether PP domain decomposition is used in this run.
 * \param[in] haveSeparatePmeRank Whether separate PME rank(s) are used in this run.
 * \param[in] useGpuForNonbonded Whether we have short-range nonbonded interactions
 *                               calculations on GPU(s).
 * \param[in] pmeRunMode         Run mode indicating what resource is PME execured on.
 * \param[in] useGpuForBonded    Whether bonded interactions are calculated on GPU(s).
 * \param[in] useGpuForUpdate    Whether coordinate update and constraint solving is performed on
 *                               GPU(s).
 * \param[in] useGpuDirectHalo   Whether halo exchange is performed directly between GPUs.
 * \param[in] canUseDirectGpuComm Whether direct GPU communication can be used in this run.
 * \param[in] useGpuPmeDecomposition GPU based PME decomposition used.
 *
 * \returns Simulation lifetime constant workload description.
 */
SimulationWorkload createSimulationWorkload(const t_inputrec& inputrec,
                                            bool              disableNonbondedCalculation,
                                            const DevelopmentFeatureFlags& devFlags,
                                            bool       havePpDomainDecomposition,
                                            bool       haveSeparatePmeRank,
                                            bool       useGpuForNonbonded,
                                            PmeRunMode pmeRunMode,
                                            bool       useGpuForBonded,
                                            bool       useGpuForUpdate,
                                            bool       useGpuDirectHalo,
                                            bool       canUseDirectGpuComm,
                                            bool       useGpuPmeDecomposition);

} // namespace gmx

#endif
