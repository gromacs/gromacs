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
 * \brief Implements update and constraints class using CUDA.
 *
 * The class combines Leap-Frog integrator with LINCS and SETTLE constraints.
 *
 * \todo The computational procedures in members should be integrated to improve
 *       computational performance.
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 *
 * \ingroup module_mdlib
 */
#include "gmxpre.h"

#include "update_constrain_cuda_impl.h"

#include <assert.h>
#include <stdio.h>

#include <cmath>

#include <algorithm>

#include "gromacs/gpu_utils/cudautils.cuh"
#include "gromacs/gpu_utils/devicebuffer.h"
#include "gromacs/gpu_utils/gputraits.cuh"
#include "gromacs/gpu_utils/vectype_ops.cuh"
#include "gromacs/mdlib/leapfrog_cuda.cuh"
#include "gromacs/mdlib/lincs_cuda.cuh"
#include "gromacs/mdlib/settle_cuda.cuh"
#include "gromacs/mdlib/update_constrain_cuda.h"

namespace gmx
{
/*!\brief Number of CUDA threads in a block
 *
 * \todo Check if using smaller block size will lead to better prformance.
 */
constexpr static int c_threadsPerBlock = 256;
//! Maximum number of threads in a block (for __launch_bounds__)
constexpr static int c_maxThreadsPerBlock = c_threadsPerBlock;

/*! \brief Scaling matrix struct.
 *
 * \todo Should be generalized.
 */
struct ScalingMatrix
{
    float xx, yy, zz, yx, zx, zy;
};

__launch_bounds__(c_maxThreadsPerBlock) __global__
        static void scaleCoordinates_kernel(const int numAtoms,
                                            float3* __restrict__ gm_x,
                                            const ScalingMatrix scalingMatrix)
{
    int threadIndex = blockIdx.x * blockDim.x + threadIdx.x;
    if (threadIndex < numAtoms)
    {
        float3 x = gm_x[threadIndex];

        x.x = scalingMatrix.xx * x.x + scalingMatrix.yx * x.y + scalingMatrix.zx * x.z;
        x.y = scalingMatrix.yy * x.y + scalingMatrix.zy * x.z;
        x.z = scalingMatrix.zz * x.z;

        gm_x[threadIndex] = x;
    }
}

void UpdateConstrainCuda::Impl::integrate(GpuEventSynchronizer*             fReadyOnDevice,
                                          const real                        dt,
                                          const bool                        updateVelocities,
                                          const bool                        computeVirial,
                                          tensor                            virial,
                                          const bool                        doTemperatureScaling,
                                          gmx::ArrayRef<const t_grp_tcstat> tcstat,
                                          const bool                        doParrinelloRahman,
                                          const float                       dtPressureCouple,
                                          const matrix                      prVelocityScalingMatrix)
{
    // Clearing virial matrix
    // TODO There is no point in having separate virial matrix for constraints
    clear_mat(virial);

    // Make sure that the forces are ready on device before proceeding with the update.
    fReadyOnDevice->enqueueWaitEvent(commandStream_);

    // The integrate should save a copy of the current coordinates in d_xp_ and write updated once
    // into d_x_. The d_xp_ is only needed by constraints.
    integrator_->integrate(d_x_, d_xp_, d_v_, d_f_, dt, doTemperatureScaling, tcstat,
                           doParrinelloRahman, dtPressureCouple, prVelocityScalingMatrix);
    // Constraints need both coordinates before (d_x_) and after (d_xp_) update. However, after constraints
    // are applied, the d_x_ can be discarded. So we intentionally swap the d_x_ and d_xp_ here to avoid the
    // d_xp_ -> d_x_ copy after constraints. Note that the integrate saves them in the wrong order as well.
    lincsCuda_->apply(d_xp_, d_x_, updateVelocities, d_v_, 1.0 / dt, computeVirial, virial);
    settleCuda_->apply(d_xp_, d_x_, updateVelocities, d_v_, 1.0 / dt, computeVirial, virial);

    // scaledVirial -> virial (methods above returns scaled values)
    float scaleFactor = 0.5f / (dt * dt);
    for (int i = 0; i < DIM; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            virial[i][j] = scaleFactor * virial[i][j];
        }
    }

    coordinatesReady_->markEvent(commandStream_);

    return;
}

void UpdateConstrainCuda::Impl::scaleCoordinates(const matrix scalingMatrix)
{
    ScalingMatrix mu;
    mu.xx = scalingMatrix[XX][XX];
    mu.yy = scalingMatrix[YY][YY];
    mu.zz = scalingMatrix[ZZ][ZZ];
    mu.yx = scalingMatrix[YY][XX];
    mu.zx = scalingMatrix[ZZ][XX];
    mu.zy = scalingMatrix[ZZ][YY];

    const auto kernelArgs = prepareGpuKernelArguments(
            scaleCoordinates_kernel, coordinateScalingKernelLaunchConfig_, &numAtoms_, &d_x_, &mu);
    launchGpuKernel(scaleCoordinates_kernel, coordinateScalingKernelLaunchConfig_, nullptr,
                    "scaleCoordinates_kernel", kernelArgs);
    // TODO: Although this only happens on the pressure coupling steps, this synchronization
    //       can affect the perfornamce if nstpcouple is small.
    gpuStreamSynchronize(commandStream_);
}

UpdateConstrainCuda::Impl::Impl(const t_inputrec&     ir,
                                const gmx_mtop_t&     mtop,
                                const void*           commandStream,
                                GpuEventSynchronizer* xUpdatedOnDevice) :
    coordinatesReady_(xUpdatedOnDevice)
{
    GMX_ASSERT(xUpdatedOnDevice != nullptr, "The event synchronizer can not be nullptr.");
    commandStream != nullptr ? commandStream_ = *static_cast<const CommandStream*>(commandStream)
                             : commandStream_ = nullptr;


    integrator_ = std::make_unique<LeapFrogCuda>(commandStream_);
    lincsCuda_  = std::make_unique<LincsCuda>(ir.nLincsIter, ir.nProjOrder, commandStream_);
    settleCuda_ = std::make_unique<SettleCuda>(mtop, commandStream_);

    coordinateScalingKernelLaunchConfig_.blockSize[0]     = c_threadsPerBlock;
    coordinateScalingKernelLaunchConfig_.blockSize[1]     = 1;
    coordinateScalingKernelLaunchConfig_.blockSize[2]     = 1;
    coordinateScalingKernelLaunchConfig_.sharedMemorySize = 0;
    coordinateScalingKernelLaunchConfig_.stream           = commandStream_;
}

UpdateConstrainCuda::Impl::~Impl() {}

void UpdateConstrainCuda::Impl::set(DeviceBuffer<float>       d_x,
                                    DeviceBuffer<float>       d_v,
                                    const DeviceBuffer<float> d_f,
                                    const t_idef&             idef,
                                    const t_mdatoms&          md,
                                    const int                 numTempScaleValues)
{
    GMX_ASSERT(d_x != nullptr, "Coordinates device buffer should not be null.");
    GMX_ASSERT(d_v != nullptr, "Velocities device buffer should not be null.");
    GMX_ASSERT(d_f != nullptr, "Forces device buffer should not be null.");

    d_x_ = reinterpret_cast<float3*>(d_x);
    d_v_ = reinterpret_cast<float3*>(d_v);
    d_f_ = reinterpret_cast<float3*>(d_f);

    numAtoms_ = md.nr;

    reallocateDeviceBuffer(&d_xp_, numAtoms_, &numXp_, &numXpAlloc_, nullptr);

    reallocateDeviceBuffer(&d_inverseMasses_, numAtoms_, &numInverseMasses_,
                           &numInverseMassesAlloc_, nullptr);

    // Integrator should also update something, but it does not even have a method yet
    integrator_->set(md, numTempScaleValues, md.cTC);
    lincsCuda_->set(idef, md);
    settleCuda_->set(idef, md);

    coordinateScalingKernelLaunchConfig_.gridSize[0] =
            (numAtoms_ + c_threadsPerBlock - 1) / c_threadsPerBlock;
}

void UpdateConstrainCuda::Impl::setPbc(const t_pbc* pbc)
{
    setPbcAiuc(pbc->ndim_ePBC, pbc->box, &pbcAiuc_);
    integrator_->setPbc(pbc);
    lincsCuda_->setPbc(pbc);
    settleCuda_->setPbc(pbc);
}

GpuEventSynchronizer* UpdateConstrainCuda::Impl::getCoordinatesReadySync()
{
    return coordinatesReady_;
}

UpdateConstrainCuda::UpdateConstrainCuda(const t_inputrec&     ir,
                                         const gmx_mtop_t&     mtop,
                                         const void*           commandStream,
                                         GpuEventSynchronizer* xUpdatedOnDevice) :
    impl_(new Impl(ir, mtop, commandStream, xUpdatedOnDevice))
{
}

UpdateConstrainCuda::~UpdateConstrainCuda() = default;

void UpdateConstrainCuda::integrate(GpuEventSynchronizer*             fReadyOnDevice,
                                    const real                        dt,
                                    const bool                        updateVelocities,
                                    const bool                        computeVirial,
                                    tensor                            virialScaled,
                                    const bool                        doTemperatureScaling,
                                    gmx::ArrayRef<const t_grp_tcstat> tcstat,
                                    const bool                        doParrinelloRahman,
                                    const float                       dtPressureCouple,
                                    const matrix                      prVelocityScalingMatrix)
{
    impl_->integrate(fReadyOnDevice, dt, updateVelocities, computeVirial, virialScaled, doTemperatureScaling,
                     tcstat, doParrinelloRahman, dtPressureCouple, prVelocityScalingMatrix);
}

void UpdateConstrainCuda::scaleCoordinates(const matrix scalingMatrix)
{
    impl_->scaleCoordinates(scalingMatrix);
}

void UpdateConstrainCuda::set(DeviceBuffer<float>       d_x,
                              DeviceBuffer<float>       d_v,
                              const DeviceBuffer<float> d_f,
                              const t_idef&             idef,
                              const t_mdatoms&          md,
                              const int                 numTempScaleValues)
{
    impl_->set(d_x, d_v, d_f, idef, md, numTempScaleValues);
}

void UpdateConstrainCuda::setPbc(const t_pbc* pbc)
{
    impl_->setPbc(pbc);
}

GpuEventSynchronizer* UpdateConstrainCuda::getCoordinatesReadySync()
{
    return impl_->getCoordinatesReadySync();
}

bool UpdateConstrainCuda::isNumCoupledConstraintsSupported(const gmx_mtop_t& mtop)
{
    return LincsCuda::isNumCoupledConstraintsSupported(mtop);
}

} // namespace gmx
