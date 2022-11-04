/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2022- The GROMACS Authors
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
 *
 * \brief Declares the MD Graph class
 *
 * \author Alan Gray <alang@nvidia.com>
 *
 *
 * \ingroup module_mdlib
 */
#ifndef GMX_MDRUN_MDGRAPH_H
#define GMX_MDRUN_MDGRAPH_H

#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/gmxmpi.h"

class GpuEventSynchronizer;

namespace gmx
{

class MdGpuGraph
{
public:
    /*! \brief Create MD graph object
     *
     *  This object manages the creation and execution of a graph of
     *  GPU activities corresponding to a MD step, allowing these
     *  multiple activities to be launched as a single entity,
     *  which reduces launch overheads. The graph is: defined by
     *  "recording" the pre-existing streams of activities (including
     *  across multiple thread-MPI tasks); converted to an executable
     *  format; and launched to execute the MD step. To record the
     *  graph, we insert some extra stream synchronization such that
     *  we can associate the start and end points of the recording
     *  with a single stream, while allowing all participating streams
     *  to fork and join. We also associate a separate graph with even
     *  and odd MD steps for 2 reasons: to cater for the the alternate
     *  non-bonded prune scheme (where local and non-local buffers
     *  are pruned on alternate steps); and to also allow the extra
     *  artificial GPU-side synchronizations to be overlapped with
     *  real work across alternate steps.
     *
     * \param [in] deviceStreamManager  Device stream manager object
     * \param [in] simulationWork       Simulation workload structure
     * \param [in] mpiComm              MPI communicator for PP domain decomposition
     * \param [in] evenOrOddStep        Whether this graph corresponds to even or odd step
     * \param [in] wcycle               Wall cycle timer object
     */
    MdGpuGraph(const DeviceStreamManager& deviceStreamManager,
               SimulationWorkload         simulationWork,
               MPI_Comm                   mpiComm,
               MdGraphEvenOrOddStep       evenOrOddStep,
               gmx_wallcycle*             wcycle);

    ~MdGpuGraph();

    /*! \brief Reset graph */
    void reset();

    /*! \brief Disable graph across all PP ranks for this domain,
     * due to presence of CPU forces on any PP rank.
     * \param [in] disableGraphAcrossAllPpRanks  Whether graph usage should be disabled across all PP ranks
     */
    void disableForDomainIfAnyPpRankHasCpuForces(bool disableGraphAcrossAllPpRanks);

    /*! \brief Decide whether graph will be captured this step
     * \param [in] canUseGraphThisStep   Whether graph can be used this step
     * \returns Whether graph will be captured this step
     */
    bool captureThisStep(bool canUseGraphThisStep);

    /*! \brief Set whether graph was used in the previous step.
     * \param [in] usedGraphLastStep     Whether graph was used in the last step
     */
    void setUsedGraphLastStep(bool usedGraphLastStep);

    /*! \brief Denote start of graph region, starts graph recording,
     * and forks those streams involved in graph.
     * \param [in] xReadyOnDeviceEvent   Event marked when coordinates are ready on device
     */
    void startRecord(GpuEventSynchronizer* xReadyOnDeviceEvent);

    /*! \brief Join streams involved in graph, and end recording */
    void endRecord();

    /*! \brief Create executable graph from recorded graph
     * \param [in] forceGraphReinstantiation Whether graph reinstantiation should be used instead of graph exec update
     */
    void createExecutableGraph(bool forceGraphReinstantiation);

    /*! \brief Launch graph corresponding to MD step
     * \param [inout] xUpdatedOnDeviceEvent  Event marked when coordinates have been updated on device
     */
    void launchGraphMdStep(GpuEventSynchronizer* xUpdatedOnDeviceEvent);

    /*! \brief Whether graph is in use this step */
    bool useGraphThisStep() const;

    /*! \brief Whether graph is capturing */
    bool graphIsCapturingThisStep() const;

    /*! \brief Set PP task completion event for graph on alternate step */
    void setAlternateStepPpTaskCompletionEvent(GpuEventSynchronizer* event);

    /*! \brief Getter for task completion event for this graph
     * \returns ppTaskCompletionEvent_
     */
    GpuEventSynchronizer* getPpTaskCompletionEvent();

private:
    class Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace gmx
#endif
