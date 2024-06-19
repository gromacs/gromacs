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
 * \brief May be used to implement MD graph CUDA interfaces for non-GPU builds.
 *
 * Currently, reports and exits if any of the interfaces are called.
 * Needed to satisfy compiler on systems, where CUDA is not available.
 *
 * \author Alan Gray <alang@nvidia.com>
 *
 *
 * \ingroup module_mdlib
 */
#include "gmxpre.h"

#include "config.h"

#include <utility>

#include "gromacs/gpu_utils/device_stream_manager.h"
#include "gromacs/mdlib/mdgraph_gpu.h"
#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/gmxmpi.h"

class GpuEventSynchronizer;
enum class MdGraphEvenOrOddStep : int;
struct gmx_wallcycle;

#if !GMX_HAVE_GPU_GRAPH_SUPPORT

namespace gmx
{

/*!\brief Impl class stub. */
class MdGpuGraph::Impl
{
};

MdGpuGraph::MdGpuGraph(const DeviceStreamManager& /* deviceStreamManager */,
                       SimulationWorkload /* simulationWork */,
                       MPI_Comm /* mpiComm */,
                       MdGraphEvenOrOddStep /* evenOrOddStep */,
                       gmx_wallcycle* /* wcycle */) :
    impl_(nullptr)
{
    GMX_ASSERT(!impl_, "A CPU stub for MD Graph was called instead of the correct implementation.");
}

MdGpuGraph::~MdGpuGraph() = default;

void MdGpuGraph::reset()
{
    GMX_ASSERT(!impl_, "A CPU stub for MD Graph was called instead of the correct implementation.");
}

void MdGpuGraph::disableForDomainIfAnyPpRankHasCpuForces(bool /* disableGraphAcrossAllPpRanks */)
{
    GMX_ASSERT(!impl_, "A CPU stub for MD Graph was called instead of the correct implementation.");
}

bool MdGpuGraph::captureThisStep(bool /* canUseGraphThisStep */)
{
    GMX_ASSERT(!impl_, "A CPU stub for MD Graph was called instead of the correct implementation.");
    return false;
}

void MdGpuGraph::setUsedGraphLastStep(bool /* usedGraphLastStep */)
{
    GMX_ASSERT(!impl_, "A CPU stub for MD Graph was called instead of the correct implementation.");
}

void MdGpuGraph::startRecord(GpuEventSynchronizer* /* xReadyOnDeviceEvent */)
{
    GMX_ASSERT(!impl_, "A CPU stub for MD Graph was called instead of the correct implementation.");
}

void MdGpuGraph::endRecord()
{
    GMX_ASSERT(!impl_, "A CPU stub for MD Graph was called instead of the correct implementation.");
}

void MdGpuGraph::createExecutableGraph(bool /* forceGraphReinstantiation */)
{
    GMX_ASSERT(!impl_, "A CPU stub for MD Graph was called instead of the correct implementation.");
}

void MdGpuGraph::launchGraphMdStep(GpuEventSynchronizer* /* xUpdatedOnDeviceEvent */)
{
    GMX_ASSERT(!impl_, "A CPU stub for MD Graph was called instead of the correct implementation.");
}

bool MdGpuGraph::useGraphThisStep() const
{
    GMX_ASSERT(!impl_, "A CPU stub for MD Graph was called instead of the correct implementation.");
    return false;
}

bool MdGpuGraph::graphIsCapturingThisStep() const
{
    GMX_ASSERT(!impl_, "A CPU stub for MD Graph was called instead of the correct implementation.");
    return false;
}

void MdGpuGraph::setAlternateStepPpTaskCompletionEvent(GpuEventSynchronizer* /* event */)
{
    GMX_ASSERT(!impl_, "A CPU stub for MD Graph was called instead of the correct implementation.");
}

GpuEventSynchronizer* MdGpuGraph::getPpTaskCompletionEvent()
{
    GMX_ASSERT(!impl_, "A CPU stub for MD Graph was called instead of the correct implementation.");
    return nullptr;
}

} // namespace gmx

#endif // !GMX_HAVE_GPU_GRAPH_SUPPORT
