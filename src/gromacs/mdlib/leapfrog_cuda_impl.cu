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
 * \brief Implements Leap-Frog using CUDA
 *
 * This file contains implementation of basic Leap-Frog integrator
 * using CUDA, including class initialization, data-structures management
 * and GPU kernel.
 *
 * \todo Reconsider naming towards using "gpu" suffix instead of "cuda".
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 *
 * \ingroup module_mdlib
 */
#include "gmxpre.h"

#include "leapfrog_cuda_impl.h"

#include <assert.h>
#include <stdio.h>

#include <cmath>

#include <algorithm>

#include "gromacs/gpu_utils/cudautils.cuh"
#include "gromacs/gpu_utils/devicebuffer.cuh"
#include "gromacs/gpu_utils/gputraits.cuh"
#include "gromacs/gpu_utils/vectype_ops.cuh"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/leapfrog_cuda.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pbcutil/pbc_aiuc_cuda.cuh"

namespace gmx
{

//! Number of CUDA threads in a block
constexpr static int c_threadsPerBlock = 256;
//! Maximum number of threads in a block (for __launch_bounds__)
constexpr static int c_maxThreadsPerBlock = c_threadsPerBlock;

/*! \brief Main kernel for Leap-Frog integrator.
 *
 *  Each GPU thread works with a single particle. Empty declaration is needed to
 *  avoid "no previous prototype for function" clang warning.
 *
 *  \todo Check if the force should be set to zero here.
 *  \todo This kernel can also accumulate incidental temperatures for each atom.
 *
 * \param[in]     numAtoms                  Total number of atoms.
 * \param[in]     gm_x                      Coordinates before the timestep
 * \param[out]    gm_xp                     Coordinates after the timestep.
 * \param[in,out] gm_v                      Velocities to update.
 * \param[in]     gm_f                      Atomic forces.
 * \param[in]     gm_inverseMasses          Reciprocal masses.
 * \param[in]     dt                        Timestep.
 */
__launch_bounds__(c_maxThreadsPerBlock)
__global__ void leapfrog_kernel(const int                  numAtoms,
                                const float3* __restrict__ gm_x,
                                float3* __restrict__       gm_xp,
                                float3* __restrict__       gm_v,
                                const float3* __restrict__ gm_f,
                                const float*  __restrict__ gm_inverseMasses,
                                const float                dt);

__launch_bounds__(c_maxThreadsPerBlock)
__global__ void leapfrog_kernel(const int                  numAtoms,
                                const float3* __restrict__ gm_x,
                                float3* __restrict__       gm_xp,
                                float3* __restrict__       gm_v,
                                const float3* __restrict__ gm_f,
                                const float*  __restrict__ gm_inverseMasses,
                                const float                dt)
{
    int threadIndex = blockIdx.x*blockDim.x + threadIdx.x;
    if (threadIndex < numAtoms)
    {
        float3 xi           = gm_x[threadIndex];
        float3 vi           = gm_v[threadIndex];
        float3 fi           = gm_f[threadIndex];
        float  imi          = gm_inverseMasses[threadIndex];
        float  imidt        = imi*dt;
        vi                 += fi*imidt;
        xi                 += vi*dt;
        gm_v[threadIndex]   = vi;
        gm_xp[threadIndex]  = xi;
    }
    return;
}

void LeapFrogCuda::Impl::integrate(const float3 *d_x,
                                   float3       *d_xp,
                                   float3       *d_v,
                                   const float3 *d_f,
                                   const real    dt)
{

    ensureNoPendingCudaError("In CUDA version of Leap-Frog integrator");

    KernelLaunchConfig config;
    config.blockSize[0]     = c_threadsPerBlock;
    config.blockSize[1]     = 1;
    config.blockSize[2]     = 1;
    config.gridSize[0]      = (numAtoms_ + c_threadsPerBlock - 1)/c_threadsPerBlock;
    config.sharedMemorySize = 0;
    config.stream           = stream_;

    auto          kernelPtr         = leapfrog_kernel;
    const float3 *gm_x              = d_x;
    float3       *gm_xp             = d_xp;
    float3       *gm_v              = d_v;
    const float3 *gm_f              = d_f;
    const float  *gm_inverseMasses  = d_inverseMasses_;

    const auto    kernelArgs = prepareGpuKernelArguments(kernelPtr, config,
                                                         &numAtoms_,
                                                         &gm_x, &gm_xp,
                                                         &gm_v,
                                                         &gm_f,
                                                         &gm_inverseMasses, &dt);
    launchGpuKernel(kernelPtr, config, nullptr, "leapfrog_kernel", kernelArgs);

    return;
}

void LeapFrogCuda::Impl::copyIntegrateCopy(const int   numAtoms,
                                           const rvec *h_x,
                                           rvec       *h_xp,
                                           rvec       *h_v,
                                           const rvec *h_f,
                                           const real  dt)
{
    float3 *d_x, *d_xp, *d_v, *d_f;

    allocateDeviceBuffer(&d_x,  numAtoms, nullptr);
    allocateDeviceBuffer(&d_xp, numAtoms, nullptr);
    allocateDeviceBuffer(&d_v,  numAtoms, nullptr);
    allocateDeviceBuffer(&d_f,  numAtoms, nullptr);

    copyToDeviceBuffer(&d_x,  (float3*)h_x,  0, numAtoms, stream_, GpuApiCallBehavior::Sync, nullptr);
    copyToDeviceBuffer(&d_xp, (float3*)h_xp, 0, numAtoms, stream_, GpuApiCallBehavior::Sync, nullptr);
    copyToDeviceBuffer(&d_v,  (float3*)h_v,  0, numAtoms, stream_, GpuApiCallBehavior::Sync, nullptr);
    copyToDeviceBuffer(&d_f,  (float3*)h_f,  0, numAtoms, stream_, GpuApiCallBehavior::Sync, nullptr);

    integrate(d_x, d_xp, d_v, d_f, dt);

    copyFromDeviceBuffer((float3*)h_xp, &d_xp, 0, numAtoms, stream_, GpuApiCallBehavior::Sync, nullptr);
    copyFromDeviceBuffer((float3*)h_v, &d_v, 0, numAtoms, stream_, GpuApiCallBehavior::Sync, nullptr);

    freeDeviceBuffer(&d_x);
    freeDeviceBuffer(&d_xp);
    freeDeviceBuffer(&d_v);
    freeDeviceBuffer(&d_f);

}

LeapFrogCuda::Impl::Impl()
{
    numAtoms_ = 0;

    // TODO When the code will be integrated into the schedule, it will be assigned non-default stream.
    stream_ = nullptr;
}

LeapFrogCuda::Impl::~Impl()
{
    freeDeviceBuffer(&d_inverseMasses_);

}

void LeapFrogCuda::Impl::setPbc(const t_pbc *pbc)
{
    setPbcAiuc(pbc->ndim_ePBC, pbc->box, &pbcAiuc_);
}

void LeapFrogCuda::Impl::set(const t_mdatoms &md)
{
    if (md.nr > numAtoms_)
    {
        if (numAtoms_ > 0)
        {
            freeDeviceBuffer(&d_inverseMasses_);
        }
        numAtoms_ = md.nr;
        allocateDeviceBuffer(&d_inverseMasses_,  numAtoms_, nullptr);
    }
    copyToDeviceBuffer(&d_inverseMasses_, (float*)md.invmass,
                       0, numAtoms_, stream_, GpuApiCallBehavior::Sync, nullptr);
}


LeapFrogCuda::LeapFrogCuda()
    : impl_(new Impl())
{
}

LeapFrogCuda::~LeapFrogCuda() = default;

void LeapFrogCuda::copyIntegrateCopy(const int   numAtoms,
                                     const rvec *h_x,
                                     rvec       *h_xp,
                                     rvec       *h_v,
                                     const rvec *h_f,
                                     const real  dt)
{
    impl_->copyIntegrateCopy(numAtoms, h_x, h_xp, h_v, h_f, dt);
}

void LeapFrogCuda::setPbc(const t_pbc *pbc)
{
    impl_->setPbc(pbc);
}

void LeapFrogCuda::set(const t_mdatoms &md)
{
    impl_->set(md);
}

} //namespace gmx
