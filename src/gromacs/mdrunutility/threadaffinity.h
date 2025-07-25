/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2012- The GROMACS Authors
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
 * \brief
 * Declares functions for managing mdrun thread affinity.
 *
 * \inlibraryapi
 * \ingroup module_mdrunutility
 */
#ifndef GMX_MDRUNUTILITY_THREADAFFINITY_H
#define GMX_MDRUNUTILITY_THREADAFFINITY_H

#include <cstdio>

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/gmxmpi.h"

struct gmx_hw_opt_t;
struct gmx_multisim_t;

namespace gmx
{

class MpiComm;
class HardwareTopology;
class MDLogger;
class PhysicalNodeCommunicator;

class IThreadAffinityAccess
{
public:
    virtual bool isThreadAffinitySupported() const        = 0;
    virtual bool setCurrentThreadAffinityToCore(int core) = 0;

protected:
    virtual ~IThreadAffinityAccess();
};

} // namespace gmx

/*! \brief Communicates within physical nodes to discover the
 * distribution of threads over ranks. */
void analyzeThreadsOnThisNode(const gmx::PhysicalNodeCommunicator& physicalNodeComm,
                              int                                  numThreadsOnThisRank,
                              int*                                 numThreadsOnThisNode,
                              int*                                 intraNodeThreadOffset);

/*! \brief
 * Sets the thread affinity using the requested setting stored in hw_opt.
 *
 * See analyzeThreadsOnThisNode(), which prepares some of the input.
 *
 * \param[out] mdlog                  Logger.
 * \param[in]  mpiCommMySim           Communication handler for the simulation.
 * \param[in]  hw_opt                 Accesses user choices for thread affinity handling.
 * \param[in]  hwTop                  Detected hardware topology.
 * \param[in]  numThreadsOnThisRank   The number of threads on this rank.
 * \param[in]  numThreadsOnThisNode   The number of threads on all ranks of this node.
 * \param[in]  intraNodeThreadOffset  The index of the first hardware thread of this rank
 *   in the set of all the threads of all MPI ranks within a node (ordered by MPI rank ID).
 * \param[in]  affinityAccess         Interface for low-level access to affinity details.
 */
void gmx_set_thread_affinity(const gmx::MDLogger&         mdlog,
                             const gmx::MpiComm&          mpiCommMySim,
                             const gmx_hw_opt_t*          hw_opt,
                             const gmx::HardwareTopology& hwTop,
                             int                          numThreadsOnThisRank,
                             int                          numThreadsOnThisNode,
                             int                          intraNodeThreadOffset,
                             gmx::IThreadAffinityAccess*  affinityAccess);

/*! \brief
 * Checks the process affinity mask and if it is found to be non-zero,
 * will honor it and disable mdrun internal affinity setting.
 *
 * This function should be called first before the OpenMP library gets
 * initialized with the last argument FALSE (which will detect affinity
 * set by external tools like taskset), and later, after the OpenMP
 * initialization, with the last argument TRUE to detect affinity changes
 * made by the OpenMP library.
 *
 * Note that this will only work on Linux as we use a GNU feature.
 * With bAfterOpenmpInit false, it will also detect whether OpenMP environment
 * variables for setting the affinity are set.
 */
void gmx_check_thread_affinity_set(const gmx::MDLogger& mdlog,
                                   gmx_hw_opt_t*        hw_opt,
                                   int                  ncpus,
                                   gmx_bool             bAfterOpenmpInit,
                                   MPI_Comm             world);

#endif
