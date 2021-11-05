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
 * \brief Implements GPU stream manager.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 *
 * \ingroup module_gpu_utils
 */
#include "gmxpre.h"

#include "device_stream_manager.h"

#include "gromacs/gpu_utils/device_context.h"
#include "gromacs/gpu_utils/device_stream.h"
#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

/*! \libinternal
 * \brief Impl class to manages the lifetime of the GPU streams.
 *
 * If supported by the GPU API, the available runtime and the
 * indicated device, some streams will be configured at high
 * priority. Otherwise, all streams will share the default priority
 * appropriate to the situation.
 */
class DeviceStreamManager::Impl
{
public:
    /*! \brief Constructor.
     *
     * \throws InternalError  If any of the required resources could not be initialized.
     */
    Impl(const DeviceInformation& deviceInfo,
         bool                     havePpDomainDecomposition,
         SimulationWorkload       simulationWork,
         bool                     useTiming);
    ~Impl();

    //! Device context.
    DeviceContext context_;
    //! GPU command streams.
    EnumerationArray<DeviceStreamType, std::unique_ptr<DeviceStream>> streams_;
};

// DeviceStreamManager::Impl
DeviceStreamManager::Impl::Impl(const DeviceInformation& deviceInfo,
                                const bool               havePpDomainDecomposition,
                                const SimulationWorkload simulationWork,
                                const bool               useTiming) :
    context_(deviceInfo)
{
    try
    {
        streams_[DeviceStreamType::NonBondedLocal] =
                std::make_unique<DeviceStream>(context_, DeviceStreamPriority::Normal, useTiming);

        if (simulationWork.useGpuPme)
        {
            /* Creating a PME GPU stream:
             * - default high priority with CUDA
             * - no priorities implemented yet with OpenCL; see #2532
             */
            streams_[DeviceStreamType::Pme] =
                    std::make_unique<DeviceStream>(context_, DeviceStreamPriority::High, useTiming);
        }

        if (havePpDomainDecomposition)
        {
            streams_[DeviceStreamType::NonBondedNonLocal] =
                    std::make_unique<DeviceStream>(context_, DeviceStreamPriority::High, useTiming);
        }
        // Update stream is used both for coordinates transfers and for GPU update/constraints
        if (simulationWork.useGpuPme || simulationWork.useGpuUpdate || simulationWork.useGpuXBufferOps)
        {
            streams_[DeviceStreamType::UpdateAndConstraints] =
                    std::make_unique<DeviceStream>(context_, DeviceStreamPriority::Normal, useTiming);
        }
        if (simulationWork.useGpuPmePpCommunication)
        {
            streams_[DeviceStreamType::PmePpTransfer] =
                    std::make_unique<DeviceStream>(context_, DeviceStreamPriority::Normal, useTiming);
        }
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
}

DeviceStreamManager::Impl::~Impl() = default;

// DeviceStreamManager
DeviceStreamManager::DeviceStreamManager(const DeviceInformation& deviceInfo,
                                         const bool               havePpDomainDecomposition,
                                         const SimulationWorkload simulationWork,
                                         const bool               useTiming) :
    impl_(new Impl(deviceInfo, havePpDomainDecomposition, simulationWork, useTiming))
{
}

DeviceStreamManager::~DeviceStreamManager() = default;

const DeviceInformation& DeviceStreamManager::deviceInfo() const
{
    return impl_->context_.deviceInfo();
}

const DeviceContext& DeviceStreamManager::context() const
{
    return impl_->context_;
}

const DeviceStream& DeviceStreamManager::stream(DeviceStreamType streamToGet) const
{
    return *impl_->streams_[streamToGet];
}

const DeviceStream& DeviceStreamManager::bondedStream(bool hasPPDomainDecomposition) const
{
    if (hasPPDomainDecomposition)
    {
        GMX_RELEASE_ASSERT(stream(DeviceStreamType::NonBondedNonLocal).isValid(),
                           "GPU non-bonded non-local stream should be valid in order to use GPU "
                           "version of bonded forces with domain decomposition.");
        return stream(DeviceStreamType::NonBondedNonLocal);
    }
    else
    {
        GMX_RELEASE_ASSERT(stream(DeviceStreamType::NonBondedLocal).isValid(),
                           "GPU non-bonded local stream should be valid in order to use GPU "
                           "version of bonded forces without domain decomposition.");
        return stream(DeviceStreamType::NonBondedLocal);
    }
}

bool DeviceStreamManager::streamIsValid(DeviceStreamType streamToCheck) const
{
    return impl_->streams_[streamToCheck] != nullptr && impl_->streams_[streamToCheck]->isValid();
}

} // namespace gmx
