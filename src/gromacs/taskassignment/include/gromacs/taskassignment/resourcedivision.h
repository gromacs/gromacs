/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2015- The GROMACS Authors
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
/*! \libinternal \file
 * \brief Declares utility functionality for dividing resources and
 * checking for consistency and usefulness.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_taskassignment
 * \inlibraryapi
 */

#ifndef GMX_TASKASSIGNMENT_RESOURCEDIVISION_H
#define GMX_TASKASSIGNMENT_RESOURCEDIVISION_H

#include <cstdio>

#include <vector>

#include "gromacs/ewald/pme.h"
#include "gromacs/utility/basedefinitions.h"

struct gmx_hw_info_t;
struct gmx_hw_opt_t;
struct gmx_mtop_t;
struct gmx_multisim_t;
struct t_inputrec;

namespace gmx
{
class MpiComm;
class HardwareTopology;
class MDLogger;
class PhysicalNodeCommunicator;
} // namespace gmx

/*! \brief Return the number of threads to use for thread-MPI based on how many
 * were requested, which algorithms we're using,
 * and how many particles there are.
 * At the point we have already called check_and_update_hw_opt.
 * Thus all options should be internally consistent and consistent
 * with the hardware, except that ntmpi could be larger than number of GPUs.
 * If necessary, this function will modify hw_opt->nthreads_omp.
 */
int get_nthreads_mpi(const gmx_hw_info_t* hwinfo,
                     gmx_hw_opt_t*        hw_opt,
                     int                  numDevicesToUse,
                     bool                 nonbondedOnGpu,
                     bool                 pmeOnGpu,
                     const t_inputrec*    inputrec,
                     const gmx_mtop_t&    mtop,
                     const gmx::MDLogger& mdlog,
                     bool                 doMembed);

/*! \brief Check if the number of OpenMP threads is within reasonable range
 * considering the hardware used. This is a crude check, but mainly
 * intended to catch cases where the user starts 1 MPI rank per hardware
 * thread or 1 rank per physical node.
 * With a sub-optimal setup a note is printed to fplog and stderr.
 * This function should be called after thread-MPI and OpenMP are set up.
 */
void check_resource_division_efficiency(const gmx_hw_info_t* hwinfo,
                                        bool                 willUsePhysicalGpu,
                                        const gmx::MpiComm*  mpiCommMySim,
                                        const gmx::MDLogger& mdlog);

/*! \brief Checks what our hardware options are based on how Gromacs was compiled
 *  and user-set options
 *
 *  \param[in]      mdlog                     Logger
 *  \param[in, out] hw_opt                    Hardware-related and threading options
 *  \param[in]      isSimulationMainRank
 *  \param[in]      nPmeRanks                 Number of PME ranks
 *  \param[in]      inputrec                  The input record, should point to a valid object when \p isSimulationMainRank = true
 *  */
void checkAndUpdateHardwareOptions(const gmx::MDLogger& mdlog,
                                   gmx_hw_opt_t*        hw_opt,
                                   bool                 isSimulationMainRank,
                                   int                  nPmeRanks,
                                   const t_inputrec*    inputrec);

/*! \brief Check, and if necessary update, the number of OpenMP threads requested
 *
 * Should be called when we know the MPI rank count and PME run mode.
 */
void checkAndUpdateRequestedNumOpenmpThreads(gmx_hw_opt_t*         hw_opt,
                                             const gmx_hw_info_t&  hwinfo,
                                             const gmx::MpiComm&   mpiCommMySim,
                                             const gmx_multisim_t* ms,
                                             int                   numRanksOnThisNode,
                                             PmeRunMode            pmeRunMode,
                                             const gmx_mtop_t&     mtop,
                                             const t_inputrec&     inputrec);

namespace gmx
{

/*! \brief Warns for oversubscribing the hardware threads, when that is the case
 */
void checkHardwareOversubscription(int                             numThreadsOnThisRank,
                                   int                             rank,
                                   const HardwareTopology&         hwTop,
                                   const PhysicalNodeCommunicator& comm,
                                   const MDLogger&                 mdlog);

} // namespace gmx

#endif
