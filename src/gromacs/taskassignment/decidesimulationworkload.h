/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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

#include "gromacs/mdtypes/simulation_workload.h"

namespace gmx
{

/*! \brief
 * Build datastructure that contains decisions whether to run different workload
 * task on GPUs.
 *
 * \param[in] useGpuForNonbonded If we have short-range nonbonded interactions
 *                               calculations on GPU(s).
 * \param[in] useGpuForPme       If long range PME interactions are calculated on GPU(s).
 * \param[in] useGpuForPmeFft    If FFT solving for PME is done on the GPU.
 * \param[in] useGpuForBonded    If bonded interactions are calculated on GPU(s).
 * \param[in] useGpuForUpdateConstraints If coordinate update and constraint solving is performed on
 *                                       GPU(s).
 * \param[in] useGpuForBufferOps If buffer ops / reduction are calculated on GPU(s).
 * \param[in] useGpuHaloExchange If GPU direct communication is used in halo exchange.
 * \param[in] useGpuPmePpComm    If GPu direct communication is used in PME-PP communication.
 * \returns Simulation lifetime constant workload description.
 */
SimulationWorkload createSimulationWorkload(bool useGpuForNonbonded,
                                            bool useGpuForPme,
                                            bool useGpuForPmeFft,
                                            bool useGpuForBonded,
                                            bool useGpuForUpdateConstraints,
                                            bool useGpuForBufferOps,
                                            bool useGpuHaloExchange,
                                            bool useGpuPmePpComm);


}  // namespace gmx

#endif
