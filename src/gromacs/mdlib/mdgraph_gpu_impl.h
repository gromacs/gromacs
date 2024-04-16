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
#ifndef GMX_MDRUN_MDGRAPH_IMPL_H
#define GMX_MDRUN_MDGRAPH_IMPL_H

#include "gromacs/gpu_utils/device_stream_manager.h"

#include "mdgraph_gpu.h"

namespace gmx
{

/*! \brief State of graph */
enum class GraphState : int
{
    Invalid,      //!< Invalid, i.e before recording has started (including after reset)
    Recording,    //!< Recording is underway and has not yet finished
    Recorded,     //!< Recording has finished, but graph is not yet instantiated
    Instantiated, //!< Instantiated and ready to launch
    Count         //!< Number of valid values
};

class MdGpuGraph::Impl
{
public:
    /*! \brief Create MD graph object
     * \param [in] deviceStreamManager  Device stream manager object
     * \param [in] simulationWork       Simulation workload structure
     * \param [in] mpiComm              MPI communicator for PP domain decomposition
     * \param [in] evenOrOddStep        Whether this graph corresponds to even or odd step
     * \param [in] wcycle               Wall cycle timer object
     */
    Impl(const DeviceStreamManager& deviceStreamManager,
         SimulationWorkload         simulationWork,
         MPI_Comm                   mpiComm,
         MdGraphEvenOrOddStep       evenOrOddStep,
         gmx_wallcycle*             wcycle);
    // NOLINTNEXTLINE(performance-trivially-destructible)
    ~Impl();

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
     * \param [inout] xUpdatedOnDeviceEvent  Event marked when coordinates have been updated o\
n device
     */
    void launchGraphMdStep(GpuEventSynchronizer* xUpdatedOnDeviceEvent);

    /*! \brief Whether graph is in use this step */
    bool useGraphThisStep() const { return useGraphThisStep_; }

    /*! \brief Whether graph is capturing */
    bool graphIsCapturingThisStep() const { return graphIsCapturingThisStep_; }

    /*! \brief Set PP task completion event for graph on alternate step */
    void setAlternateStepPpTaskCompletionEvent(GpuEventSynchronizer* event);

    /*! \brief Getter for task completion event for this graph
     * \returns ppTaskCompletionEvent_
     */
    GpuEventSynchronizer* getPpTaskCompletionEvent();

#if GMX_GPU_CUDA
    using Graph         = cudaGraph_t;
    using GraphInstance = cudaGraphExec_t;
#else
    using Graph         = void*;
    using GraphInstance = void*;
#endif

private:
    /*! \brief Collective operation to enqueue events from all PP ranks to a stream on PP rank 0
     * \param [in] event   Event to enqueue, valid on all PP ranks
     * \param [in] stream  Stream to enqueue events, valid on PP rank 0
     */
    void enqueueEventFromAllPpRanksToRank0Stream(GpuEventSynchronizer* event, const DeviceStream& stream);


    /*! \brief Collective operation to enqueue an event from PP rank 0 to streams on all PP ranks
     * \param [in] event   Event to enqueue, valid on PP rank 0
     * \param [in] stream  Stream to enqueue events, valid on all PP ranks
     */
    void enqueueRank0EventToAllPpStreams(GpuEventSynchronizer* event, const DeviceStream& stream);

    //! Captured graph object
    Graph graph_;
    //! Instantiated graph object
    GraphInstance instance_;
    //! Whether graph has already been created
    bool graphCreated_ = false;
    //! Whether graph is capturing in this step
    bool graphIsCapturingThisStep_ = false;
    //! Whether graph should be used this step
    bool useGraphThisStep_ = false;
    //! Whether graph was used in the last step
    bool usedGraphLastStep_ = false;
    //! Whether graph should be disabled for this domain (including if due to disablement on other rank.)
    bool disableGraphAcrossAllPpRanks_ = false;
    //! Device stream manager object
    const DeviceStreamManager& deviceStreamManager_;
    //! Device stream for launching graph
    std::unique_ptr<DeviceStream> launchStream_;
    //! Alternate Device stream for launching graph
    std::unique_ptr<DeviceStream> launchStreamAlternate_;
    //! Whether PP domain decomposition is in use
    bool havePPDomainDecomposition_;
    //! Whether there is PME GPU task on this PP rank
    bool haveGpuPmeOnThisPpRank_;
    //! Whether PME is handled on a separate rank
    bool haveSeparatePmeRank_;
    //! MPI communicator for PP domain decomposition
    MPI_Comm mpiComm_;
    //! PP Rank for this MD graph object
    int ppRank_ = 0;
    //! Number of PP ranks in use
    int ppSize_ = 1;
    //! Temporary event used for forking and joining streams in graph
    std::unique_ptr<GpuEventSynchronizer> helperEvent_;
    //! Whether step is even or odd, where different graphs are used for each
    MdGraphEvenOrOddStep evenOrOddStep_;
    //! event marked on this step when this PP task has completed its tasks
    std::unique_ptr<GpuEventSynchronizer> ppTaskCompletionEvent_;
    //! event marked on alternate step when PP task has completed its tasks
    GpuEventSynchronizer* alternateStepPpTaskCompletionEvent_;
    //! Wall cycle timer object
    gmx_wallcycle* wcycle_;
    //! Whether the graph object has been allocated
    bool graphAllocated_ = false;
    //! Whether the graph instance object has been allocated
    bool graphInstanceAllocated_ = false;
    //! State of graph
    GraphState graphState_ = GraphState::Invalid;
    //! Whether a perormance bug workaround is needed in graph update/reinstantiation
    bool needOldDriverTransferWorkaround_ = false;
};

} // namespace gmx
#endif
