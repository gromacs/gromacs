/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2024- The GROMACS Authors
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
 * \brief Defines the MD Graph class
 *
 * \author Alan Gray <alang@nvidia.com>
 * \author Andrey Alekseenko <al42and@gmail.com>
 *
 *
 * \ingroup module_mdlib
 */

#include "gmxpre.h"

#include "gromacs/gpu_utils/device_context.h"
#include "gromacs/gpu_utils/device_stream.h"
#include "gromacs/gpu_utils/gmxsycl.h"
#include "gromacs/gpu_utils/gpueventsynchronizer.h"
#include "gromacs/hardware/device_information.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxmpi.h"

#include "mdgraph_gpu_impl.h"

#if GMX_HAVE_GPU_GRAPH_SUPPORT

namespace syclex = sycl::ext::oneapi::experimental;

namespace gmx
{

MdGpuGraph::Impl::Impl(const DeviceStreamManager& deviceStreamManager,
                       SimulationWorkload         simulationWork,
                       MPI_Comm                   mpiComm,
                       MdGraphEvenOrOddStep       evenOrOddStep,
                       gmx_wallcycle*             wcycle) :
    deviceStreamManager_(deviceStreamManager),
    launchStream_(new DeviceStream(deviceStreamManager.context(), DeviceStreamPriority::Normal, false)),
    launchStreamAlternate_(nullptr),
    havePPDomainDecomposition_(simulationWork.havePpDomainDecomposition),
    haveGpuPmeOnThisPpRank_(simulationWork.haveGpuPmeOnPpRank()),
    haveSeparatePmeRank_(simulationWork.haveSeparatePmeRank),
    mpiComm_(mpiComm),
    evenOrOddStep_(evenOrOddStep),
    wcycle_(wcycle)
{
    helperEvent_           = std::make_unique<GpuEventSynchronizer>();
    ppTaskCompletionEvent_ = std::make_unique<GpuEventSynchronizer>();

    if (!deviceStreamManager.deviceInfo().supportsSyclGraph)
    {
        GMX_THROW(gmx::InconsistentInputError(
                "SYCL Graphs requested on a device that does not support them."));
    }

    GMX_RELEASE_ASSERT(!havePPDomainDecomposition_,
                       "GPU Graphs with multiple ranks require threadMPI with GPU-direct "
                       "communication, but it is not supported in SYCL");
    GMX_RELEASE_ASSERT(!haveSeparatePmeRank_,
                       "GPU Graphs with separate PME rank require threadMPI with GPU-direct "
                       "communication, but it is not supported in SYCL");

    // Not used in SYCL since threadMPI direct comm is not supported
    GMX_UNUSED_VALUE(mpiComm_);
    GMX_UNUSED_VALUE(ppSize_);
    GMX_UNUSED_VALUE(evenOrOddStep_);
    // NVIDIA-specific, not used in SYCL
    GMX_UNUSED_VALUE(needOldDriverTransferWorkaround_);
}

MdGpuGraph::Impl::~Impl()
{
    launchStream_->synchronize();
}


void MdGpuGraph::Impl::enqueueEventFromAllPpRanksToRank0Stream(GpuEventSynchronizer*, const DeviceStream&)
{
    GMX_RELEASE_ASSERT(false, "enqueueEventFromAllPpRanksToRank0Stream not supported in SYCL");
}

void MdGpuGraph::Impl::enqueueRank0EventToAllPpStreams(GpuEventSynchronizer*, const DeviceStream&)
{
    GMX_RELEASE_ASSERT(false, "enqueueRank0EventToAllPpStreams not supported in SYCL");
}

void MdGpuGraph::Impl::reset()
{
    graphCaptureStarted_      = false;
    useGraphThisStep_         = false;
    graphIsCapturingThisStep_ = false;
    graphState_               = GraphState::Invalid;
}

void MdGpuGraph::Impl::disableForDomainIfAnyPpRankHasCpuForces(bool disableGraphAcrossAllPpRanks)
{
    disableGraphAcrossAllPpRanks_ = disableGraphAcrossAllPpRanks;
}

bool MdGpuGraph::Impl::captureThisStep(bool canUseGraphThisStep)
{
    useGraphThisStep_         = canUseGraphThisStep && !disableGraphAcrossAllPpRanks_;
    graphIsCapturingThisStep_ = useGraphThisStep_ && !graphCaptureStarted_;
    return graphIsCapturingThisStep_;
}

void MdGpuGraph::Impl::setUsedGraphLastStep(bool usedGraphLastStep)
{
    usedGraphLastStep_ = usedGraphLastStep;
}

void MdGpuGraph::Impl::startRecord(GpuEventSynchronizer* xReadyOnDeviceEvent)
{
    GMX_ASSERT(useGraphThisStep_,
               "startRecord should not have been called if graph is not in use this step");
    GMX_ASSERT(graphIsCapturingThisStep_,
               "startRecord should not have been called if graph is not capturing this step");
    GMX_ASSERT(graphState_ == GraphState::Invalid,
               "Graph should be in an invalid state before recording");
    GMX_RELEASE_ASSERT(ppRank_ == 0, "SYCL Graph does not support recording with PP decomposition");

    wallcycle_start(wcycle_, WallCycleCounter::MdGpuGraph);
    wallcycle_sub_start(wcycle_, WallCycleSubCounter::MdGpuGraphCapture);

    // Disable check for cycles in graph
    static const sycl::property_list s_propList{ syclex::property::graph::no_cycle_check() };
    graph_ = std::make_unique<syclex::command_graph<syclex::graph_state::modifiable>>(
            deviceStreamManager_.context().context(), deviceStreamManager_.deviceInfo().syclDevice, s_propList);
    graphCaptureStarted_ = true;

    std::vector<sycl::queue> queuesToRecord;
    queuesToRecord.emplace_back(deviceStreamManager_.stream(gmx::DeviceStreamType::NonBondedLocal).stream());
    queuesToRecord.emplace_back(
            deviceStreamManager_.stream(gmx::DeviceStreamType::UpdateAndConstraints).stream());
    if (haveGpuPmeOnThisPpRank_)
    {
        queuesToRecord.emplace_back(deviceStreamManager_.stream(gmx::DeviceStreamType::Pme).stream());
    }

    graph_->begin_recording(queuesToRecord);

    // Re-mark xReadyOnDeviceEvent to allow full isolation within graph capture
    // We explicitly want to replace an existing event, so we call the reset here
    if (xReadyOnDeviceEvent->isMarked())
    {
        xReadyOnDeviceEvent->reset();
    }
    xReadyOnDeviceEvent->markEvent(deviceStreamManager_.stream(gmx::DeviceStreamType::UpdateAndConstraints));

    graphState_ = GraphState::Recording;
};


void MdGpuGraph::Impl::endRecord()
{

    GMX_ASSERT(useGraphThisStep_,
               "endRecord should not have been called if graph is not in use this step");
    GMX_ASSERT(graphIsCapturingThisStep_,
               "endRecord should not have been called if graph is not capturing this step");
    GMX_ASSERT(graphState_ == GraphState::Recording,
               "Graph should be in a recording state before recording is ended");


    graph_->end_recording();

    graphState_ = GraphState::Recorded;

    // Sync all tasks before closing timing region, since the graph capture should be treated as a collective operation for timing purposes.
    wallcycle_sub_stop(wcycle_, WallCycleSubCounter::MdGpuGraphCapture);
    wallcycle_stop(wcycle_, WallCycleCounter::MdGpuGraph);
};

void MdGpuGraph::Impl::createExecutableGraph(bool forceGraphReinstantiation)
{

    GMX_ASSERT(
            useGraphThisStep_,
            "createExecutableGraph should not have been called if graph is not in use this step");
    GMX_ASSERT(graphIsCapturingThisStep_,
               "createExecutableGraph should not have been called if graph is not capturing this "
               "step");
    GMX_ASSERT(graphState_ == GraphState::Recorded,
               "Graph should be in a recorded state before instantiation");

    // graph::update  API exists, but it is not supported, so we always re-instantiate
    GMX_UNUSED_VALUE(forceGraphReinstantiation);

    wallcycle_start(wcycle_, WallCycleCounter::MdGpuGraph);
    wallcycle_sub_start(wcycle_, WallCycleSubCounter::MdGpuGraphInstantiateOrUpdate);

    // Instantiate graph
    const auto instance = graph_->finalize();
    instance_ = std::make_unique<syclex::command_graph<syclex::graph_state::executable>>(std::move(instance));

    graphState_ = GraphState::Instantiated;

    wallcycle_sub_stop(wcycle_, WallCycleSubCounter::MdGpuGraphInstantiateOrUpdate);
    wallcycle_stop(wcycle_, WallCycleCounter::MdGpuGraph);
};

void MdGpuGraph::Impl::launchGraphMdStep(GpuEventSynchronizer* xUpdatedOnDeviceEvent)
{

    GMX_ASSERT(useGraphThisStep_,
               "launchGraphMdStep should not have been called if graph is not in use this step");
    GMX_ASSERT(graphState_ == GraphState::Instantiated,
               "Graph should be in an instantiated state before launching");

    wallcycle_start(wcycle_, WallCycleCounter::MdGpuGraph);
    wallcycle_sub_start(wcycle_, WallCycleSubCounter::MdGpuGraphLaunch);

    const DeviceStream* thisLaunchStream = launchStream_.get();

    if (!usedGraphLastStep_)
    {

        // Sync update and constraints
        helperEvent_->markEvent(deviceStreamManager_.stream(gmx::DeviceStreamType::UpdateAndConstraints));
        helperEvent_->enqueueWaitEvent(*thisLaunchStream);

        // Sync NB local and non-local (to ensure no race condition with pruning)
        helperEvent_->markEvent(deviceStreamManager_.stream(gmx::DeviceStreamType::NonBondedLocal));
        helperEvent_->enqueueWaitEvent(*thisLaunchStream);

        // If PME on same rank, sync PME (to ensure no race condition with clearing)
        // Note that separate rank PME has implicit sync, including clearing.
        if (haveGpuPmeOnThisPpRank_)
        {
            helperEvent_->markEvent(deviceStreamManager_.stream(gmx::DeviceStreamType::Pme));
            helperEvent_->enqueueWaitEvent(*thisLaunchStream);
        }
    }


    if (alternateStepPpTaskCompletionEvent_->isMarked())
    {
        alternateStepPpTaskCompletionEvent_->enqueueWaitEvent(
                *thisLaunchStream); // TODO: Suboptimal, but we can't have external nodes yet (2024-09-23).
    }
    thisLaunchStream->stream().ext_oneapi_graph(*instance_);
    if (ppTaskCompletionEvent_->isMarked())
    {
        ppTaskCompletionEvent_->reset();
    }
    ppTaskCompletionEvent_->markEvent(
            *thisLaunchStream); // TODO: Suboptimal, but we can't have external nodes yet (2024-09-23).

    if (xUpdatedOnDeviceEvent->isMarked())
    {
        xUpdatedOnDeviceEvent->reset();
    }
    xUpdatedOnDeviceEvent->markEvent(*thisLaunchStream);

    wallcycle_sub_stop(wcycle_, WallCycleSubCounter::MdGpuGraphLaunch);
    wallcycle_stop(wcycle_, WallCycleCounter::MdGpuGraph);
}

void MdGpuGraph::Impl::setAlternateStepPpTaskCompletionEvent(GpuEventSynchronizer* event)
{
    alternateStepPpTaskCompletionEvent_ = event;
}

GpuEventSynchronizer* MdGpuGraph::Impl::getPpTaskCompletionEvent()
{
    return ppTaskCompletionEvent_.get();
}

MdGpuGraph::MdGpuGraph(const DeviceStreamManager& deviceStreamManager,
                       SimulationWorkload         simulationWork,
                       MPI_Comm                   mpiComm,
                       MdGraphEvenOrOddStep       evenOrOddStep,
                       gmx_wallcycle*             wcycle) :
    impl_(new Impl(deviceStreamManager, simulationWork, mpiComm, evenOrOddStep, wcycle))
{
}

MdGpuGraph::~MdGpuGraph() = default;

void MdGpuGraph::reset()
{
    impl_->reset();
}

void MdGpuGraph::disableForDomainIfAnyPpRankHasCpuForces(bool disableGraphAcrossAllPpRanks)
{
    impl_->disableForDomainIfAnyPpRankHasCpuForces(disableGraphAcrossAllPpRanks);
}

bool MdGpuGraph::captureThisStep(bool canUseGraphThisStep)
{
    return impl_->captureThisStep(canUseGraphThisStep);
}

void MdGpuGraph::setUsedGraphLastStep(bool usedGraphLastStep)
{
    impl_->setUsedGraphLastStep(usedGraphLastStep);
}

void MdGpuGraph::startRecord(GpuEventSynchronizer* xReadyOnDeviceEvent)
{
    impl_->startRecord(xReadyOnDeviceEvent);
}

void MdGpuGraph::endRecord()
{
    impl_->endRecord();
}

void MdGpuGraph::createExecutableGraph(bool forceGraphReinstantiation)
{
    impl_->createExecutableGraph(forceGraphReinstantiation);
}

void MdGpuGraph::launchGraphMdStep(GpuEventSynchronizer* xUpdatedOnDeviceEvent)
{
    impl_->launchGraphMdStep(xUpdatedOnDeviceEvent);
}

bool MdGpuGraph::useGraphThisStep() const
{
    return impl_->useGraphThisStep();
}

bool MdGpuGraph::graphIsCapturingThisStep() const
{
    return impl_->graphIsCapturingThisStep();
}

void MdGpuGraph::setAlternateStepPpTaskCompletionEvent(GpuEventSynchronizer* event)
{
    impl_->setAlternateStepPpTaskCompletionEvent(event);
}

GpuEventSynchronizer* MdGpuGraph::getPpTaskCompletionEvent()
{
    return impl_->getPpTaskCompletionEvent();
}

} // namespace gmx

#endif // GMX_HAVE_GPU_GRAPH_SUPPORT
