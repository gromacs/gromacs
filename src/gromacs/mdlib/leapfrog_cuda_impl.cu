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
 * \param[in]     nAtom                     Total number atoms.
 * \param[in]     x                         Coordinates before the timestep
 * \param[out]    xp                        Coordinates after the timestep.
 * \param[in,out] v                         Velocities to update.
 * \param[in]     pbcAiuc                   Periodic boundary information.
 * \param[in]     inverseMasses             Reciprocal masses.
 * \param[in]     dt                        Timestep.
 */
__global__ void leapfrog_kernel(const int             nAtom,
                                const float3         *x,
                                float3               *xp,
                                float3               *v,
                                const float3         *f,
                                const float          *inverseMasses,
                                const float           dt)
{
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i < nAtom)
    {
        float3 xi    = x[i];
        float3 vi    = v[i];
        float3 fi    = f[i];
        float  imi   = inverseMasses[i];
        float  imidt = imi*dt;
        vi   += fi*imidt;
        xi   += vi*dt;
        v[i]  = vi;
        xp[i] = xi;
    }

    return;
}

/*! \brief Integrate
 *
 * Integrates the equation of motion using Leap-Frog algorithm.
 * Updates xpDevice_ and vDevice_ fields of this object.
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

    auto          kernelPtr               = leapfrog_kernel;
    const float3 *xDevicePtr              = xDevice_;
    float3       *xpDevicePtr             = xpDevice_;
    float3       *vDevicePtr              = xDevice_;
    const float3 *fDevicePtr              = xDevice_;
    const float  *inverseMassesDevicePtr  = inverseMassesDevice_;

    const auto    kernelArgs = prepareGpuKernelArguments(kernelPtr, config,
                                                         &nAtom_,
                                                         &xDevicePtr, &xpDevicePtr,
                                                         &vDevicePtr,
                                                         &fDevicePtr,
                                                         &inverseMassesDevicePtr, &dt);
    launchGpuKernel(kernelPtr, config, nullptr, "leapfrog_kernel", kernelArgs);

    cudaError_t stat = cudaGetLastError();
    CU_RET_ERR(stat, "In CUDA version of Leap-Frog integrator");

    return;
}

/*! \brief Create Leap-Frog object
 *
 * \param [in] nAtom  Number of atoms.
 */
LeapFrogCuda::Impl::Impl(int nAtom)
    : nAtom_(nAtom)
{
    allocateDeviceBuffer(&xDevice_,              nAtom, nullptr);
    allocateDeviceBuffer(&xpDevice_,             nAtom, nullptr);
    allocateDeviceBuffer(&vDevice_,              nAtom, nullptr);
    allocateDeviceBuffer(&fDevice_,              nAtom, nullptr);
    allocateDeviceBuffer(&inverseMassesDevice_,  nAtom, nullptr);

    stream_ = nullptr;
    cudaError_t stat = cudaGetLastError();
    CU_RET_ERR(stat, "While CUDA version of Leap-Frog was initialized.");
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
 * \param[in] *inverseMasses  CPU pointer to inverse masses to copy from.
 */
void LeapFrogCuda::Impl::copyInverseMassesToGpu(const real *inverseMasses)
{
    copyToDeviceBuffer(&inverseMassesDevice_, (float*)inverseMasses,
                       0, nAtom_, stream_, GpuApiCallBehavior::Sync, nullptr);
    cudaError_t stat = cudaGetLastError();
    CU_RET_ERR(stat, "While inverse masses were copied to GPU in CUDA version of Leap-Frog integrator.");
}

/*! \brief
 * Copy coordinates from CPU to GPU.
 *
 * The data are assumed to be in float3/fvec format (single precision).
 *
 * \param[in] *x  CPU pointer where coordinates should be copied from.
 */
void LeapFrogCuda::Impl::copyCoordinatesToGpu(const rvec *x)
{
    copyToDeviceBuffer(&xDevice_, (float3*)x, 0, nAtom_, stream_, GpuApiCallBehavior::Sync, nullptr);
    cudaError_t stat = cudaGetLastError();
    CU_RET_ERR(stat, "While coordinates were copied to GPU in CUDA version of Leap-Frog integrator.");
}

/*! \brief
 * Copy velocities from CPU to GPU.
 *
 * The data are assumed to be in float3/fvec format (single precision).
 *
 * \param[in] *v  CPU pointer where velocities should be copied from.
 */
void LeapFrogCuda::Impl::copyVelocitiesToGpu(const rvec *v)
{
    copyToDeviceBuffer(&vDevice_, (float3*)v, 0, nAtom_, stream_, GpuApiCallBehavior::Sync, nullptr);
    cudaError_t stat = cudaGetLastError();
    CU_RET_ERR(stat, "While velocities were copied to GPU in CUDA version of Leap-Frog integrator.");
}

/*! \brief
 * Copy forces from CPU to GPU.
 *
 * The data are assumed to be in float3/fvec format (single precision).
 *
 * \param[in] *f  CPU pointer where forces should be copied from.
 */
void LeapFrogCuda::Impl::copyForcesToGpu(const rvec *f)
{
    copyToDeviceBuffer(&fDevice_, (float3*)f, 0, nAtom_, stream_, GpuApiCallBehavior::Sync, nullptr);
    cudaError_t stat = cudaGetLastError();
    CU_RET_ERR(stat, "While forces were copied to GPU in CUDA version of Leap-Frog integrator.");
}

/*! \brief
 * Copy coordinates from GPU to CPU.
 *
 * The data are assumed to be in float3/fvec format (single precision).
 *
 * \param[out] *xp CPU pointer where coordinates should be copied to.
 */
void LeapFrogCuda::Impl::copyCoordinatesFromGpu(rvec *xp)
{
    copyFromDeviceBuffer((float3*)xp, &xpDevice_, 0, nAtom_, stream_, GpuApiCallBehavior::Sync, nullptr);
    cudaError_t stat = cudaGetLastError();
    CU_RET_ERR(stat, "While coordinates were copied from GPU in CUDA version of Leap-Frog integrator.");
}

/*! \brief
 * Copy velocities from GPU to CPU.
 *
 * The velocities are assumed to be in float3/fvec format (single precision).
 *
 * \param[in] *v  Pointer to velocities data.
 */
void LeapFrogCuda::Impl::copyVelocitiesFromGpu(rvec *v)
{
    copyFromDeviceBuffer((float3*)v, &vDevice_, 0, nAtom_, stream_, GpuApiCallBehavior::Sync, nullptr);
    cudaError_t stat = cudaGetLastError();
    CU_RET_ERR(stat, "While velocities were copied from GPU in CUDA version of Leap-Frog integrator.");
}

/*! \brief
 * Copy forces from GPU to CPU.
 *
 * The forces are assumed to be in float3/fvec format (single precision).
 *
 * \param[in] *f  Pointer to forces data.
 */
void LeapFrogCuda::Impl::copyForcesFromGpu(rvec *f)
{
    copyFromDeviceBuffer((float3*)f, &fDevice_, 0, nAtom_, stream_, GpuApiCallBehavior::Sync, nullptr);
    cudaError_t stat = cudaGetLastError();
    CU_RET_ERR(stat, "While forces were copied from GPU in CUDA version of Leap-Frog integrator.");
}

/*! \brief
 * Set the internal GPU-memory x, xprime and v pointers.
 *
 * Data is not copied. The data are assumed to be in float3/fvec format
 * (float3 is used internally, but the data layout should be identical).
 *
 * \param[in] *xDevice  Pointer to the coordinates for the input (on GPU)
 * \param[in] *xpDevice Pointer to the coordinates for the output (on GPU)
 * \param[in] *vDevice  Pointer to the velocities (on GPU)
 * \param[in] *fDevice  Pointer to the forces (on GPU)
 */
void LeapFrogCuda::Impl::setXVFPointers(rvec *xDevice, rvec *xpDevice, rvec *vDevice, rvec *fDevice)
{
    xDevice_  = (float3*)xDevice;
    xpDevice_ = (float3*)xpDevice;
    vDevice_  = (float3*)vDevice;
    fDevice_  = (float3*)fDevice;
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

void LeapFrogCuda::copyInverseMassesToGpu(const real *inverseMasses)
{
    impl_->copyInverseMassesToGpu(inverseMasses);
}

void LeapFrogCuda::copyCoordinatesToGpu(const rvec *x)
{
    impl_->copyCoordinatesToGpu(x);
}

void LeapFrogCuda::copyVelocitiesToGpu(const rvec *v)
{
    impl_->copyVelocitiesToGpu(v);
}

void LeapFrogCuda::copyForcesToGpu(const rvec *f)
{
    impl_->copyForcesToGpu(f);
}

void LeapFrogCuda::copyCoordinatesFromGpu(rvec *xp)
{
    impl_->copyCoordinatesFromGpu(xp);
}

void LeapFrogCuda::copyVelocitiesFromGpu(rvec *v)
{
    impl_->copyVelocitiesFromGpu(v);
}

void LeapFrogCuda::copyForcesFromGpu(rvec *f)
{
    impl_->copyForcesFromGpu(f);
}

void LeapFrogCuda::setXVFPointers(rvec *xDevice, rvec *xpDevice, rvec *vDevice, rvec *fDevice)
{
    impl_->setXVFPointers(xDevice, xpDevice, vDevice, fDevice);
}

} //namespace gmx
