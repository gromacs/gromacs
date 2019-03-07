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

constexpr static int c_threadsPerBlock = 256;

/*! \brief Main kernel for Leap-Frog integrator.
 *
 *  Each GPU thread works with a single particle.
 *  \todo Check if the force should be set to zero here.
 *  \todo This kernel can also accumulate incindental temperatures for each atom.
 *
 * \param[in]     nAtom                     Total number atoms.
 * \param[in]     gm_x                      Coordinates before the timestep
 * \param[out]    gm_xp                     Coordinates after the timestep.
 * \param[in,out] gm_v                      Velocities to update.
 * \param[in]     pbcAiuc                   Periodic boundary information.
 * \param[in]     gm_inverseMasses          Reciprocal masses.
 * \param[in]     dt                        Timestep.
 */
__global__ void leapfrog_kernel(const int                  nAtom,
                                const float3* __restrict__ gm_x,
                                float3* __restrict__       gm_xp,
                                float3* __restrict__       gm_v,
                                const float3* __restrict__ gm_f,
                                const float*  __restrict__ gm_inverseMasses,
                                const float                dt)
{
    int threadIndex = blockIdx.x*blockDim.x + threadIdx.x;
    if (threadIndex < nAtom)
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

/*! \brief Integrate
 *
 * Integrates the equation of motion using Leap-Frog algorithm.
 * Updates d_xp_ and d_v_ fields of this object.
 *
 * \param[in] dt             Timestep
 */
void LeapFrogCuda::Impl::integrate(const real  dt)
{

    ensureNoPendingCudaError("In CUDA version of Leap-Frog integrator");

    KernelLaunchConfig config;
    config.blockSize[0]     = c_threadsPerBlock;
    config.blockSize[1]     = 1;
    config.blockSize[2]     = 1;
    config.gridSize[0]      = (nAtom_ + c_threadsPerBlock - 1)/c_threadsPerBlock;
    config.sharedMemorySize = 1;
    config.stream           = stream_;

    auto          kernelPtr        = leapfrog_kernel;
    const float3 *d_x              = d_x_;
    float3       *d_xp             = d_xp_;
    float3       *d_v              = d_v_;
    const float3 *d_f              = d_f_;
    const float  *d_inverseMasses  = d_inverseMasses_;

    const auto    kernelArgs = prepareGpuKernelArguments(kernelPtr, config,
                                                         &nAtom_,
                                                         &d_x, &d_xp,
                                                         &d_v,
                                                         &d_f,
                                                         &d_inverseMasses, &dt);
    launchGpuKernel(kernelPtr, config, nullptr, "leapfrog_kernel", kernelArgs);

    return;
}

/*! \brief Create Leap-Frog object
 *
 * \param [in] nAtom  Number of atoms.
 */
LeapFrogCuda::Impl::Impl(int nAtom)
    : nAtom_(nAtom)
{
    allocateDeviceBuffer(&d_x_,              nAtom, nullptr);
    allocateDeviceBuffer(&d_xp_,             nAtom, nullptr);
    allocateDeviceBuffer(&d_v_,              nAtom, nullptr);
    allocateDeviceBuffer(&d_f_,              nAtom, nullptr);
    allocateDeviceBuffer(&d_inverseMasses_,  nAtom, nullptr);

    // \todo When the code will be integrated into the schedule, it will be assigned non-default stream.
    stream_ = nullptr;
}

LeapFrogCuda::Impl::~Impl()
{
}

/*! \brief
 * Update PBC data.
 *
 * Converts pbc data from t_pbc into the PbcAiuc format and stores the latter.
 *
 * \param[in] *pbc The PBC data in t_pbc format.
 */
void LeapFrogCuda::Impl::setPbc(const t_pbc *pbc)
{
    setPbcAiuc(pbc->ndim_ePBC, pbc->box, &pbcAiuc_);
}

/*! \brief
 * Copy inverse masses from CPU to GPU.
 *
 * The data are assumed to be in float format (single precision).
 *
 * \param[in] *h_inverseMasses  CPU pointer to inverse masses to copy from.
 */
void LeapFrogCuda::Impl::copyInverseMassesToGpu(const real *h_inverseMasses)
{
    copyToDeviceBuffer(&d_inverseMasses_, (float*)h_inverseMasses,
                       0, nAtom_, stream_, GpuApiCallBehavior::Sync, nullptr);
}

/*! \brief
 * Copy coordinates from CPU to GPU.
 *
 * The data are assumed to be in float3/fvec format (single precision).
 *
 * \param[in] *h_x  CPU pointer where coordinates should be copied from.
 */
void LeapFrogCuda::Impl::copyCoordinatesToGpu(const rvec *h_x)
{
    copyToDeviceBuffer(&d_x_, (float3*)h_x, 0, nAtom_, stream_, GpuApiCallBehavior::Sync, nullptr);
}

/*! \brief
 * Copy velocities from CPU to GPU.
 *
 * The data are assumed to be in float3/fvec format (single precision).
 *
 * \param[in] *h_v  CPU pointer where velocities should be copied from.
 */
void LeapFrogCuda::Impl::copyVelocitiesToGpu(const rvec *h_v)
{
    copyToDeviceBuffer(&d_v_, (float3*)h_v, 0, nAtom_, stream_, GpuApiCallBehavior::Sync, nullptr);
}

/*! \brief
 * Copy forces from CPU to GPU.
 *
 * The data are assumed to be in float3/fvec format (single precision).
 *
 * \param[in] *h_f  CPU pointer where forces should be copied from.
 */
void LeapFrogCuda::Impl::copyForcesToGpu(const rvec *h_f)
{
    copyToDeviceBuffer(&d_f_, (float3*)h_f, 0, nAtom_, stream_, GpuApiCallBehavior::Sync, nullptr);
}

/*! \brief
 * Copy coordinates from GPU to CPU.
 *
 * The data are assumed to be in float3/fvec format (single precision).
 *
 * \param[out] *h_xp CPU pointer where coordinates should be copied to.
 */
void LeapFrogCuda::Impl::copyCoordinatesFromGpu(rvec *h_xp)
{
    copyFromDeviceBuffer((float3*)h_xp, &d_xp_, 0, nAtom_, stream_, GpuApiCallBehavior::Sync, nullptr);
}

/*! \brief
 * Copy velocities from GPU to CPU.
 *
 * The velocities are assumed to be in float3/fvec format (single precision).
 *
 * \param[in] *h_v  Pointer to velocities data.
 */
void LeapFrogCuda::Impl::copyVelocitiesFromGpu(rvec *h_v)
{
    copyFromDeviceBuffer((float3*)h_v, &d_v_, 0, nAtom_, stream_, GpuApiCallBehavior::Sync, nullptr);
}

/*! \brief
 * Copy forces from GPU to CPU.
 *
 * The forces are assumed to be in float3/fvec format (single precision).
 *
 * \param[in] *h_f  Pointer to forces data.
 */
void LeapFrogCuda::Impl::copyForcesFromGpu(rvec *h_f)
{
    copyFromDeviceBuffer((float3*)h_f, &d_f_, 0, nAtom_, stream_, GpuApiCallBehavior::Sync, nullptr);
}

/*! \brief
 * Set the internal GPU-memory x, xprime and v pointers.
 *
 * Data is not copied. The data are assumed to be in float3/fvec format
 * (float3 is used internally, but the data layout should be identical).
 *
 * \param[in] *d_x  Pointer to the coordinates for the input (on GPU)
 * \param[in] *d_xp Pointer to the coordinates for the output (on GPU)
 * \param[in] *d_v  Pointer to the velocities (on GPU)
 * \param[in] *d_f  Pointer to the forces (on GPU)
 */
void LeapFrogCuda::Impl::setXVFPointers(rvec *d_x, rvec *d_xp, rvec *d_v, rvec *d_f)
{
    d_x_  = (float3*)d_x;
    d_xp_ = (float3*)d_xp;
    d_v_  = (float3*)d_v;
    d_f_  = (float3*)d_f;
}


LeapFrogCuda::LeapFrogCuda(const int nAtom)
    : impl_(new Impl(nAtom))
{
}

LeapFrogCuda::~LeapFrogCuda() = default;

void LeapFrogCuda::integrate(const real  dt)
{
    impl_->integrate(dt);
}

void LeapFrogCuda::setPbc(const t_pbc *pbc)
{
    impl_->setPbc(pbc);
}

void LeapFrogCuda::copyInverseMassesToGpu(const real *h_inverseMasses)
{
    impl_->copyInverseMassesToGpu(h_inverseMasses);
}

void LeapFrogCuda::copyCoordinatesToGpu(const rvec *h_x)
{
    impl_->copyCoordinatesToGpu(h_x);
}

void LeapFrogCuda::copyVelocitiesToGpu(const rvec *h_v)
{
    impl_->copyVelocitiesToGpu(h_v);
}

void LeapFrogCuda::copyForcesToGpu(const rvec *h_f)
{
    impl_->copyForcesToGpu(h_f);
}

void LeapFrogCuda::copyCoordinatesFromGpu(rvec *h_xp)
{
    impl_->copyCoordinatesFromGpu(h_xp);
}

void LeapFrogCuda::copyVelocitiesFromGpu(rvec *h_v)
{
    impl_->copyVelocitiesFromGpu(h_v);
}

void LeapFrogCuda::copyForcesFromGpu(rvec *h_f)
{
    impl_->copyForcesFromGpu(h_f);
}

void LeapFrogCuda::setXVFPointers(rvec *d_x, rvec *d_xp, rvec *d_v, rvec *d_f)
{
    impl_->setXVFPointers(d_x, d_xp, d_v, d_f);
}

} //namespace gmx
