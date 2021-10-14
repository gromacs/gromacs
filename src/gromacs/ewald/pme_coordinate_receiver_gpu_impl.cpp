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
 *
 * \brief May be used to implement PME-PP GPU comm interfaces for non-GPU builds.
 *
 * Currently, reports and exits if any of the interfaces are called.
 * Needed to satisfy compiler on systems, where CUDA is not available.
 *
 * \author Alan Gray <alang@nvidia.com>
 *
 * \ingroup module_ewald
 */
#include "gmxpre.h"

#include "config.h"

#include "gromacs/ewald/pme_coordinate_receiver_gpu.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/gmxassert.h"

#if !GMX_GPU_CUDA

namespace gmx
{

/*!\brief \internal Impl class stub. */
class PmeCoordinateReceiverGpu::Impl
{
};

/*!\brief Constructor stub. */
PmeCoordinateReceiverGpu::PmeCoordinateReceiverGpu(MPI_Comm /* comm */,
                                                   const DeviceContext& /* deviceContext */,
                                                   gmx::ArrayRef<PpRanks> /* ppRanks */) :
    impl_(nullptr)
{
    GMX_ASSERT(!impl_,
               "A CPU stub for PME-PP GPU communication was called instead of the correct "
               "implementation.");
}

PmeCoordinateReceiverGpu::~PmeCoordinateReceiverGpu() = default;

/*!\brief init PME-PP GPU communication stub */
void PmeCoordinateReceiverGpu::reinitCoordinateReceiver(DeviceBuffer<RVec> /* d_x */)
{
    GMX_ASSERT(!impl_,
               "A CPU stub for PME-PP GPU communication initialization was called instead of the "
               "correct implementation.");
}

void PmeCoordinateReceiverGpu::receiveCoordinatesSynchronizerFromPpCudaDirect(int /* ppRank */)
{
    GMX_ASSERT(!impl_,
               "A CPU stub for PME-PP GPU communication was called instead of the correct "
               "implementation.");
}

void PmeCoordinateReceiverGpu::launchReceiveCoordinatesFromPpCudaMpi(DeviceBuffer<RVec> /* recvbuf */,
                                                                     int /* numAtoms */,
                                                                     int /* numBytes */,
                                                                     int /* ppRank */)
{
    GMX_ASSERT(!impl_,
               "A CPU stub for PME-PP GPU communication was called instead of the correct "
               "implementation.");
}

int PmeCoordinateReceiverGpu::synchronizeOnCoordinatesFromPpRank(int /* pipelineStage */,
                                                                 const DeviceStream& /* deviceStream */)
{
    GMX_ASSERT(!impl_,
               "A CPU stub for PME-PP GPU communication was called instead of the correct "
               "implementation.");
    return 0;
}

void PmeCoordinateReceiverGpu::synchronizeOnCoordinatesFromAllPpRanks(const DeviceStream& /* deviceStream */)
{
    GMX_ASSERT(!impl_,
               "A CPU stub for PME-PP GPU communication was called instead of the correct "
               "implementation.");
}

DeviceStream* PmeCoordinateReceiverGpu::ppCommStream(int /* senderIndex */)
{
    GMX_ASSERT(!impl_,
               "A CPU stub for PME-PP GPU communication was called instead of the correct "
               "implementation.");
    return nullptr;
}

std::tuple<int, int> PmeCoordinateReceiverGpu::ppCommAtomRange(int /* senderIndex */)
{
    GMX_ASSERT(!impl_,
               "A CPU stub for PME-PP GPU communication was called instead of the correct "
               "implementation.");
    return std::make_tuple(0, 0);
}

int PmeCoordinateReceiverGpu::ppCommNumSenderRanks()
{
    GMX_ASSERT(!impl_,
               "A CPU stub for PME-PP GPU communication was called instead of the correct "
               "implementation.");
    return 0;
}


} // namespace gmx

#endif // !GMX_GPU_CUDA
