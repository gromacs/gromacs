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
/*! \internal \file
 *
 * \brief Definitions of interfaces for GPU state data propagator object.
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 *
 * \ingroup module_mdtypes
 */
#include "gmxpre.h"

#include "config.h"

#if GMX_GPU != GMX_GPU_NONE

#if GMX_GPU == GMX_GPU_CUDA
#include "gromacs/gpu_utils/cudautils.cuh"
#endif
#include "gromacs/gpu_utils/devicebuffer.h"
#if GMX_GPU == GMX_GPU_OPENCL
#include "gromacs/gpu_utils/oclutils.h"
#endif
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/state_propagator_data_gpu.h"
#include "gromacs/utility/classhelpers.h"

#include "state_propagator_data_gpu_impl.h"

namespace gmx
{

StatePropagatorDataGpu::Impl::Impl(gmx_unused const void *commandStream,
                                   gmx_unused const void *deviceContext,
                                   GpuApiCallBehavior     transferKind,
                                   int                    paddingSize) :
    transferKind_(transferKind),
    paddingSize_(paddingSize)
{

    GMX_RELEASE_ASSERT(getenv("GMX_USE_GPU_BUFFER_OPS") == nullptr, "GPU buffer ops are not supported in this build.");

    // Set the stream-context pair for the OpenCL builds,
    // use the nullptr stream for CUDA builds
#if GMX_GPU == GMX_GPU_OPENCL
    if (commandStream != nullptr)
    {
        commandStream_ = *static_cast<const CommandStream*>(commandStream);
    }
    if (deviceContext != nullptr)
    {
        deviceContext_ = *static_cast<const DeviceContext*>(deviceContext);
    }
#endif

}

StatePropagatorDataGpu::Impl::~Impl()
{
}

void StatePropagatorDataGpu::Impl::reinit(int numAtomsLocal, int numAtomsAll)
{
#if GMX_GPU == GMX_GPU_OPENCL
    GMX_ASSERT(deviceContext_ != nullptr, "GPU context should be set in OpenCL builds.");
#endif
    numAtomsLocal_ = numAtomsLocal;
    numAtomsAll_   = numAtomsAll;

    int numAtomsPadded;
    if (paddingSize_ > 0)
    {
        numAtomsPadded = ((numAtomsAll_ + paddingSize_ - 1 ) / paddingSize_ )*paddingSize_;
    }
    else
    {
        numAtomsPadded = numAtomsAll_;
    }

    reallocateDeviceBuffer(&d_x_, DIM*numAtomsPadded, &d_xSize_, &d_xCapacity_, deviceContext_);

    const size_t paddingAllocationSize = numAtomsPadded - numAtomsAll_;
    if (paddingAllocationSize > 0)
    {
        clearDeviceBufferAsync(&d_x_, DIM*numAtomsAll_, DIM*paddingAllocationSize, commandStream_);
    }

    reallocateDeviceBuffer(&d_v_, DIM*numAtomsAll_, &d_vSize_, &d_vCapacity_, deviceContext_);
    reallocateDeviceBuffer(&d_f_, DIM*numAtomsAll_, &d_fSize_, &d_fCapacity_, deviceContext_);

}

std::tuple<int, int> StatePropagatorDataGpu::Impl::getAtomRangesFromAtomLocality(AtomLocality  atomLocality)
{
    int atomsStartAt   = 0;
    int numAtomsToCopy = 0;
    switch (atomLocality)
    {
        case AtomLocality::All:
            atomsStartAt    = 0;
            numAtomsToCopy  = numAtomsAll_;
            break;
        case AtomLocality::Local:
            atomsStartAt    = 0;
            numAtomsToCopy  = numAtomsLocal_;
            break;
        case AtomLocality::NonLocal:
            atomsStartAt    = numAtomsLocal_;
            numAtomsToCopy  = numAtomsAll_ - numAtomsLocal_;
            break;
        default:
            GMX_RELEASE_ASSERT(false, "Wrong range of atoms requested in GPU state data manager. Should be All, Local or NonLocal.");
    }
    GMX_ASSERT(atomsStartAt   >= 0, "The first elemtnt to copy has negative index. Probably, the GPU propagator state was not initialized.");
    GMX_ASSERT(numAtomsToCopy >= 0, "Number of atoms to copy is negative. Probably, the GPU propagator state was not initialized.");
    return std::make_tuple(atomsStartAt, numAtomsToCopy);
}

void StatePropagatorDataGpu::Impl::copyToDevice(DeviceBuffer<float>                   d_data,
                                                const gmx::ArrayRef<const gmx::RVec>  h_data,
                                                int                                   dataSize,
                                                AtomLocality                          atomLocality)
{

#if GMX_GPU == GMX_GPU_OPENCL
    GMX_ASSERT(deviceContext_ != nullptr, "GPU context should be set in OpenCL builds.");
#endif

    GMX_UNUSED_VALUE(dataSize);

    GMX_ASSERT(dataSize >= 0, "Trying to copy to device buffer before it was allocated.");

    int atomsStartAt, numAtomsToCopy;
    std::tie(atomsStartAt, numAtomsToCopy) = getAtomRangesFromAtomLocality(atomLocality);

    int elementsStartAt   = atomsStartAt*DIM;
    int numElementsToCopy = numAtomsToCopy*DIM;

    if (numAtomsToCopy != 0)
    {
        GMX_ASSERT(elementsStartAt + numElementsToCopy <= dataSize, "The device allocation is smaller than requested copy range.");
        GMX_ASSERT(atomsStartAt + numAtomsToCopy <= h_data.ssize(), "The host buffer is smaller than the requested copy range.");

        // TODO: Use the proper stream
        copyToDeviceBuffer(&d_data, reinterpret_cast<const float *>(&h_data.data()[atomsStartAt]),
                           elementsStartAt, numElementsToCopy,
                           commandStream_, transferKind_, nullptr);
    }
}

void StatePropagatorDataGpu::Impl::copyFromDevice(gmx::ArrayRef<gmx::RVec>  h_data,
                                                  DeviceBuffer<float>       d_data,
                                                  int                       dataSize,
                                                  AtomLocality              atomLocality)
{

#if GMX_GPU == GMX_GPU_OPENCL
    GMX_ASSERT(deviceContext_ != nullptr, "GPU context should be set in OpenCL builds.");
#endif

    GMX_UNUSED_VALUE(dataSize);

    GMX_ASSERT(dataSize >= 0, "Trying to copy from device buffer before it was allocated.");

    int atomsStartAt, numAtomsToCopy;
    std::tie(atomsStartAt, numAtomsToCopy) = getAtomRangesFromAtomLocality(atomLocality);

    int elementsStartAt   = atomsStartAt*DIM;
    int numElementsToCopy = numAtomsToCopy*DIM;

    if (numAtomsToCopy != 0)
    {
        GMX_ASSERT(elementsStartAt + numElementsToCopy <= dataSize, "The device allocation is smaller than requested copy range.");
        GMX_ASSERT(atomsStartAt + numAtomsToCopy <= h_data.ssize(), "The host buffer is smaller than the requested copy range.");

        // TODO: Use the proper stream
        copyFromDeviceBuffer(reinterpret_cast<float*>(&h_data.data()[atomsStartAt]), &d_data,
                             elementsStartAt, numElementsToCopy,
                             commandStream_, transferKind_, nullptr);

    }
}

DeviceBuffer<float> StatePropagatorDataGpu::Impl::getCoordinates()
{
    return d_x_;
}

void StatePropagatorDataGpu::Impl::copyCoordinatesToGpu(const gmx::ArrayRef<const gmx::RVec>  h_x,
                                                        AtomLocality                          atomLocality)
{
    copyToDevice(d_x_, h_x, d_xSize_, atomLocality);
}

void StatePropagatorDataGpu::Impl::copyCoordinatesFromGpu(gmx::ArrayRef<gmx::RVec>  h_x,
                                                          AtomLocality              atomLocality)
{
    copyFromDevice(h_x, d_x_, d_xSize_, atomLocality);
}


DeviceBuffer<float> StatePropagatorDataGpu::Impl::getVelocities()
{
    return d_v_;
}

void StatePropagatorDataGpu::Impl::copyVelocitiesToGpu(const gmx::ArrayRef<const gmx::RVec>  h_v,
                                                       AtomLocality                          atomLocality)
{
    copyToDevice(d_v_, h_v, d_vSize_, atomLocality);
}

void StatePropagatorDataGpu::Impl::copyVelocitiesFromGpu(gmx::ArrayRef<gmx::RVec>  h_v,
                                                         AtomLocality              atomLocality)
{
    copyFromDevice(h_v, d_v_, d_vSize_, atomLocality);
}


DeviceBuffer<float> StatePropagatorDataGpu::Impl::getForces()
{
    return d_f_;
}

void StatePropagatorDataGpu::Impl::copyForcesToGpu(const gmx::ArrayRef<const gmx::RVec>  h_f,
                                                   AtomLocality                          atomLocality)
{
    copyToDevice(d_f_, h_f, d_fSize_, atomLocality);
}

void StatePropagatorDataGpu::Impl::copyForcesFromGpu(gmx::ArrayRef<gmx::RVec>  h_f,
                                                     AtomLocality              atomLocality)
{
    copyFromDevice(h_f, d_f_, d_fSize_, atomLocality);
}

void StatePropagatorDataGpu::Impl::synchronizeStream()
{
    gpuStreamSynchronize(commandStream_);
}

int StatePropagatorDataGpu::Impl::numAtomsLocal()
{
    return numAtomsLocal_;
}

int StatePropagatorDataGpu::Impl::numAtomsAll()
{
    return numAtomsAll_;
}



StatePropagatorDataGpu::StatePropagatorDataGpu(const void        *commandStream,
                                               const void        *deviceContext,
                                               GpuApiCallBehavior transferKind,
                                               int                paddingSize)
    : impl_(new Impl(commandStream,
                     deviceContext,
                     transferKind,
                     paddingSize))
{
}

StatePropagatorDataGpu::StatePropagatorDataGpu(StatePropagatorDataGpu && /* other */) noexcept = default;

StatePropagatorDataGpu &StatePropagatorDataGpu::operator=(StatePropagatorDataGpu && /* other */) noexcept = default;

StatePropagatorDataGpu::~StatePropagatorDataGpu() = default;


void StatePropagatorDataGpu::reinit(int numAtomsLocal, int numAtomsAll)
{
    return impl_->reinit(numAtomsLocal, numAtomsAll);
}

std::tuple<int, int> StatePropagatorDataGpu::getAtomRangesFromAtomLocality(AtomLocality  atomLocality)
{
    return impl_->getAtomRangesFromAtomLocality(atomLocality);
}


DeviceBuffer<float> StatePropagatorDataGpu::getCoordinates()
{
    return impl_->getCoordinates();
}

void StatePropagatorDataGpu::copyCoordinatesToGpu(const gmx::ArrayRef<const gmx::RVec>  h_x,
                                                  AtomLocality                          atomLocality)
{
    return impl_->copyCoordinatesToGpu(h_x, atomLocality);
}

void StatePropagatorDataGpu::copyCoordinatesFromGpu(gmx::ArrayRef<RVec>  h_x,
                                                    AtomLocality         atomLocality)
{
    return impl_->copyCoordinatesFromGpu(h_x, atomLocality);
}


DeviceBuffer<float> StatePropagatorDataGpu::getVelocities()
{
    return impl_->getVelocities();
}

void StatePropagatorDataGpu::copyVelocitiesToGpu(const gmx::ArrayRef<const gmx::RVec>  h_v,
                                                 AtomLocality                          atomLocality)
{
    return impl_->copyVelocitiesToGpu(h_v, atomLocality);
}

void StatePropagatorDataGpu::copyVelocitiesFromGpu(gmx::ArrayRef<RVec>  h_v,
                                                   AtomLocality         atomLocality)
{
    return impl_->copyVelocitiesFromGpu(h_v, atomLocality);
}


DeviceBuffer<float> StatePropagatorDataGpu::getForces()
{
    return impl_->getForces();
}

void StatePropagatorDataGpu::copyForcesToGpu(const gmx::ArrayRef<const gmx::RVec>  h_f,
                                             AtomLocality                          atomLocality)
{
    return impl_->copyForcesToGpu(h_f, atomLocality);
}

void StatePropagatorDataGpu::copyForcesFromGpu(gmx::ArrayRef<RVec>  h_f,
                                               AtomLocality         atomLocality)
{
    return impl_->copyForcesFromGpu(h_f, atomLocality);
}

void StatePropagatorDataGpu::synchronizeStream()
{
    return impl_->synchronizeStream();
}

int StatePropagatorDataGpu::numAtomsLocal()
{
    return impl_->numAtomsLocal();
}

int StatePropagatorDataGpu::numAtomsAll()
{
    return impl_->numAtomsAll();
}

}      // namespace gmx

#endif // GMX_GPU == GMX_GPU_NONE
