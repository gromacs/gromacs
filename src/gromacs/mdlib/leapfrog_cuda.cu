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

#include "leapfrog_cuda.cuh"

#include <assert.h>
#include <stdio.h>

#include <cmath>

#include <algorithm>

#include "gromacs/gpu_utils/cudautils.cuh"
#include "gromacs/gpu_utils/devicebuffer.h"
#include "gromacs/gpu_utils/vectype_ops.cuh"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pbcutil/pbc_aiuc_cuda.cuh"
#include "gromacs/utility/arrayref.h"

namespace gmx
{

//! Number of CUDA threads in a block
constexpr static int c_threadsPerBlock = 256;
//! Maximum number of threads in a block (for __launch_bounds__)
constexpr static int c_maxThreadsPerBlock = c_threadsPerBlock;

/*! \brief Sets the number of different temperature coupling values
 *
 *  This is needed to template the kernel
 *  \todo Unify with similar enum in CPU update module
 */
enum class NumTempScaleValues
{
    None,      //!< No temperature coupling
    Single,    //!< Single T-scaling value (one group)
    Multiple   //!< Multiple T-scaling values, need to use T-group indices
};

/*! \brief Different variants of the Parrinello-Rahman velocity scaling
 *
 *  This is needed to template the kernel
 *  \todo Unify with similar enum in CPU update module
 */
enum class VelocityScalingType
{
    None,      //!< Do not apply velocity scaling (not a PR-coupling run or step)
    Diagonal,  //!< Apply velocity scaling using a diagonal matrix
    Full       //!< Apply velocity scaling using a full matrix
};

/*! \brief Main kernel for Leap-Frog integrator.
 *
 *  Each GPU thread works with a single particle. Empty declaration is needed to
 *  avoid "no previous prototype for function" clang warning.
 *
 *  \todo Check if the force should be set to zero here.
 *  \todo This kernel can also accumulate incidental temperatures for each atom.
 *
 * \tparam        numTempScaleValues             The number of different T-couple values.
 * \tparam        velocityScaling                Type of the Parrinello-Rahman velocity rescaling.
 * \param[in]     numAtoms                       Total number of atoms.
 * \param[in]     gm_x                           Coordinates before the timestep
 * \param[out]    gm_xp                          Coordinates after the timestep.
 * \param[in,out] gm_v                           Velocities to update.
 * \param[in]     gm_f                           Atomic forces.
 * \param[in]     gm_inverseMasses               Reciprocal masses.
 * \param[in]     dt                             Timestep.
 * \param[in]     gm_lambdas                     Temperature scaling factors (one per group)
 * \param[in]     gm_tempScaleGroups             Mapping of atoms into groups.
 * \param[in]     dtPressureCouple               Time step for pressure coupling
 * \param[in]     velocityScalingMatrixDiagonal  Diagonal elements of Parrinello-Rahman velocity scaling matrix
 */
template<NumTempScaleValues numTempScaleValues, VelocityScalingType velocityScaling>
__launch_bounds__(c_maxThreadsPerBlock)
__global__ void leapfrog_kernel(const int                             numAtoms,
                                const float3* __restrict__            gm_x,
                                float3* __restrict__                  gm_xp,
                                float3* __restrict__                  gm_v,
                                const float3* __restrict__            gm_f,
                                const float*  __restrict__            gm_inverseMasses,
                                const float                           dt,
                                const float*  __restrict__            gm_lambdas,
                                const unsigned short*    __restrict__ gm_tempScaleGroups,
                                const float3                          velocityScalingMatrixDiagonal);

template<NumTempScaleValues numTempScaleValues, VelocityScalingType velocityScaling>
__launch_bounds__(c_maxThreadsPerBlock)
__global__ void leapfrog_kernel(const int                             numAtoms,
                                const float3* __restrict__            gm_x,
                                float3* __restrict__                  gm_xp,
                                float3* __restrict__                  gm_v,
                                const float3* __restrict__            gm_f,
                                const float*  __restrict__            gm_inverseMasses,
                                const float                           dt,
                                const float*  __restrict__            gm_lambdas,
                                const unsigned short*    __restrict__ gm_tempScaleGroups,
                                const float3                          velocityScalingMatrixDiagonal)
{
    int threadIndex = blockIdx.x*blockDim.x + threadIdx.x;
    if (threadIndex < numAtoms)
    {
        float3 x           = gm_x[threadIndex];
        float3 v           = gm_v[threadIndex];
        float3 f           = gm_f[threadIndex];
        float  im          = gm_inverseMasses[threadIndex];
        float  imdt        = im*dt;

        if (numTempScaleValues != NumTempScaleValues::None || velocityScaling != VelocityScalingType::None)
        {
            float3 vp = v;

            if (numTempScaleValues != NumTempScaleValues::None)
            {
                float lambda = 1.0F;
                if (numTempScaleValues == NumTempScaleValues::Single)
                {
                    lambda    = gm_lambdas[0];
                }
                else if (numTempScaleValues == NumTempScaleValues::Multiple)
                {
                    int tempScaleGroup = gm_tempScaleGroups[threadIndex];
                    lambda             = gm_lambdas[tempScaleGroup];
                }
                vp *= lambda;
            }

            if (velocityScaling == VelocityScalingType::Diagonal)
            {
                vp.x -= velocityScalingMatrixDiagonal.x*v.x;
                vp.y -= velocityScalingMatrixDiagonal.y*v.y;
                vp.z -= velocityScalingMatrixDiagonal.z*v.z;
            }

            v = vp;
        }

        v += f*imdt;

        x                  += v*dt;
        gm_v[threadIndex]   = v;
        gm_xp[threadIndex]  = x;
    }
    return;
}

/*! \brief Select templated kernel.
 *
 * Returns pointer to a CUDA kernel based on the number of temperature coupling groups.
 * If zero is passed as an argument, it is assumed that no temperature coupling groups are used.
 *
 * \param[in]  numTempScaleValues  Numer of temperature coupling groups in the system
 * \param[in]  velocityScaling     Type of the Parrinello-Rahman velocity scaling
 *
 * \retrun                         Pointer to CUDA kernel
 */
inline auto selectLeapFrogKernelPtr(int numTempScaleValues, VelocityScalingType velocityScaling)
{
    auto kernelPtr = leapfrog_kernel<NumTempScaleValues::None, VelocityScalingType::None>;

    if (velocityScaling == VelocityScalingType::None)
    {
        if (numTempScaleValues == 0)
        {
            kernelPtr = leapfrog_kernel<NumTempScaleValues::None, VelocityScalingType::None>;
        }
        else if (numTempScaleValues == 1)
        {
            kernelPtr = leapfrog_kernel<NumTempScaleValues::Single, VelocityScalingType::None>;
        }
        else if (numTempScaleValues > 1)
        {
            kernelPtr = leapfrog_kernel<NumTempScaleValues::Multiple, VelocityScalingType::None>;
        }
        else
        {
            GMX_RELEASE_ASSERT(false, "Number of temperature coupling groups should be greater than zero (zero for no coupling).");
        }
    }
    else
    if (velocityScaling == VelocityScalingType::Diagonal)
    {
        if (numTempScaleValues == 0)
        {
            kernelPtr = leapfrog_kernel<NumTempScaleValues::None, VelocityScalingType::Diagonal>;
        }
        else if (numTempScaleValues == 1)
        {
            kernelPtr = leapfrog_kernel<NumTempScaleValues::Single, VelocityScalingType::Diagonal>;
        }
        else if (numTempScaleValues > 1)
        {
            kernelPtr = leapfrog_kernel<NumTempScaleValues::Multiple, VelocityScalingType::Diagonal>;
        }
        else
        {
            GMX_RELEASE_ASSERT(false, "Number of temperature coupling groups should be greater than zero (zero for no coupling).");
        }
    }
    else
    {
        GMX_RELEASE_ASSERT(false, "Only isotropic Parrinello-Rahman pressure coupling is supported.");
    }
    return kernelPtr;
}

void LeapFrogCuda::integrate(const float3                      *d_x,
                             float3                            *d_xp,
                             float3                            *d_v,
                             const float3                      *d_f,
                             const real                         dt,
                             const bool                         doTempCouple,
                             gmx::ArrayRef<const t_grp_tcstat>  tcstat,
                             const bool                         doPressureCouple,
                             const float                        dtPressureCouple,
                             const matrix                       velocityScalingMatrix)
{

    ensureNoPendingCudaError("In CUDA version of Leap-Frog integrator");

    auto kernelPtr          = leapfrog_kernel<NumTempScaleValues::None, VelocityScalingType::None>;
    if (doTempCouple || doPressureCouple)
    {
        if (doTempCouple)
        {
            GMX_ASSERT(numTempScaleValues_ == ssize(h_lambdas_), "Number of temperature scaling factors changed since it was set for the last time.");
            for (int i = 0; i < numTempScaleValues_; i++)
            {
                h_lambdas_[i] = tcstat[i].lambda;
            }
            copyToDeviceBuffer(&d_lambdas_, h_lambdas_.data(),
                               0, numTempScaleValues_, stream_, GpuApiCallBehavior::Async, nullptr);
        }
        VelocityScalingType velocityScaling = VelocityScalingType::None;
        if (doPressureCouple)
        {
            velocityScaling = VelocityScalingType::Diagonal;
            GMX_ASSERT(velocityScalingMatrix[YY][XX] == 0 && velocityScalingMatrix[ZZ][XX] == 0 && velocityScalingMatrix[ZZ][YY] == 0 &&
                       velocityScalingMatrix[XX][YY] == 0 && velocityScalingMatrix[XX][ZZ] == 0 && velocityScalingMatrix[YY][ZZ] == 0,
                       "Fully anisotropic Parrinello-Rahman pressure coupling is not yet supported in GPU version of Leap-Frog integrator.");
            velocityScalingMatrixDiagonal_ = make_float3(dtPressureCouple*velocityScalingMatrix[XX][XX],
                                                         dtPressureCouple*velocityScalingMatrix[YY][YY],
                                                         dtPressureCouple*velocityScalingMatrix[ZZ][ZZ]);
        }
        kernelPtr = selectLeapFrogKernelPtr(numTempScaleValues_, velocityScaling);
    }

    const auto            kernelArgs = prepareGpuKernelArguments(kernelPtr, kernelLaunchConfig_,
                                                                 &numAtoms_,
                                                                 &d_x, &d_xp,
                                                                 &d_v,
                                                                 &d_f,
                                                                 &d_inverseMasses_, &dt,
                                                                 &d_lambdas_, &d_tempScaleGroups_,
                                                                 &velocityScalingMatrixDiagonal_);
    launchGpuKernel(kernelPtr, kernelLaunchConfig_, nullptr, "leapfrog_kernel", kernelArgs);

    return;
}

LeapFrogCuda::LeapFrogCuda()
{
    numAtoms_ = 0;

    // TODO When the code is integrated into the schedule, it should be assigned non-default stream.
    stream_ = nullptr;

    changePinningPolicy(&h_lambdas_, gmx::PinningPolicy::PinnedIfSupported);

    kernelLaunchConfig_.blockSize[0]     = c_threadsPerBlock;
    kernelLaunchConfig_.blockSize[1]     = 1;
    kernelLaunchConfig_.blockSize[2]     = 1;
    kernelLaunchConfig_.sharedMemorySize = 0;
    kernelLaunchConfig_.stream           = stream_;
}

LeapFrogCuda::~LeapFrogCuda()
{
    freeDeviceBuffer(&d_inverseMasses_);

}

void LeapFrogCuda::setPbc(const t_pbc *pbc)
{
    setPbcAiuc(pbc->ndim_ePBC, pbc->box, &pbcAiuc_);
}

void LeapFrogCuda::set(const t_mdatoms      &md,
                       const int             numTempScaleValues,
                       const unsigned short *tempScaleGroups)
{
    numAtoms_                       = md.nr;
    kernelLaunchConfig_.gridSize[0] = (numAtoms_ + c_threadsPerBlock - 1)/c_threadsPerBlock;

    numTempScaleValues_ = numTempScaleValues;

    reallocateDeviceBuffer(&d_inverseMasses_, numAtoms_, &numInverseMasses_, &numInverseMassesAlloc_, nullptr);
    copyToDeviceBuffer(&d_inverseMasses_, (float*)md.invmass,
                       0, numAtoms_, stream_, GpuApiCallBehavior::Sync, nullptr);

    // Temperature scale group map only used if there are more then one group
    if (numTempScaleValues > 1)
    {
        reallocateDeviceBuffer(&d_tempScaleGroups_, numAtoms_, &numTempScaleGroups_, &numTempScaleGroupsAlloc_, nullptr);
        copyToDeviceBuffer(&d_tempScaleGroups_, tempScaleGroups,
                           0, numAtoms_, stream_, GpuApiCallBehavior::Sync, nullptr);
    }

    // If the temperature coupling is enabled, we need to make space for scaling factors
    if (numTempScaleValues_ > 0)
    {
        h_lambdas_.resize(numTempScaleValues);
        reallocateDeviceBuffer(&d_lambdas_, numTempScaleValues_, &numLambdas_, &numLambdasAlloc_, nullptr);
    }

}

} //namespace gmx
