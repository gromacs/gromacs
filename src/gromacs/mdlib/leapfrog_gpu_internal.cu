/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
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
 * \brief Implements Leap-Frog using CUDA
 *
 * This file contains CUDA implementation of back-end specific code for Leap-Frog.
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 *
 * \ingroup module_mdlib
 */
#include "gmxpre.h"

#include "leapfrog_gpu_internal.h"

#include "gromacs/gpu_utils/cudautils.cuh"
#include "gromacs/gpu_utils/devicebuffer.h"
#include "gromacs/gpu_utils/typecasts.cuh"
#include "gromacs/gpu_utils/vectype_ops.cuh"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/group.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pbcutil/pbc_aiuc_cuda.cuh"
#include "gromacs/utility/arrayref.h"

namespace gmx
{

/*!\brief Number of CUDA threads in a block
 *
 * \todo Check if using smaller block size will lead to better performance.
 */
constexpr static int c_threadsPerBlock = 256;
//! Maximum number of threads in a block (for __launch_bounds__)
constexpr static int c_maxThreadsPerBlock = c_threadsPerBlock;

/*! \brief Main kernel for Leap-Frog integrator.
 *
 *  The coordinates and velocities are updated on the GPU. Also saves the intermediate values of the coordinates for
 *   further use in constraints.
 *
 *  Each GPU thread works with a single particle. Empty declaration is needed to
 *  avoid "no previous prototype for function" clang warning.
 *
 *  \todo Check if the force should be set to zero here.
 *  \todo This kernel can also accumulate incidental temperatures for each atom.
 *
 * \tparam        numTempScaleValues               The number of different T-couple values.
 * \tparam        parrinelloRahmanVelocityScaling  The properties of the Parrinello-Rahman velocity scaling matrix.
 * \param[in]     numAtoms                         Total number of atoms.
 * \param[in,out] gm_x                             Coordinates to update upon integration.
 * \param[out]    gm_xp                            A copy of the coordinates before the integration (for constraints).
 * \param[in,out] gm_v                             Velocities to update.
 * \param[in]     gm_f                             Atomic forces.
 * \param[in]     gm_inverseMasses                 Reciprocal masses.
 * \param[in]     dt                               Timestep.
 * \param[in]     gm_lambdas                       Temperature scaling factors (one per group)
 * \param[in]     gm_tempScaleGroups               Mapping of atoms into groups.
 * \param[in]     prVelocityScalingMatrixDiagonal  Diagonal elements of Parrinello-Rahman velocity scaling matrix
 */
template<NumTempScaleValues numTempScaleValues, ParrinelloRahmanVelocityScaling parrinelloRahmanVelocityScaling>
__launch_bounds__(c_maxThreadsPerBlock) __global__
        void leapfrog_kernel(const int numAtoms,
                             float3* __restrict__ gm_x,
                             float3* __restrict__ gm_xp,
                             float3* __restrict__ gm_v,
                             const float3* __restrict__ gm_f,
                             const float* __restrict__ gm_inverseMasses,
                             const float dt,
                             const float* __restrict__ gm_lambdas,
                             const unsigned short* __restrict__ gm_tempScaleGroups,
                             const float3 prVelocityScalingMatrixDiagonal)
{
    int threadIndex = blockIdx.x * blockDim.x + threadIdx.x;
    if (threadIndex < numAtoms)
    {
        float3 x    = gm_x[threadIndex];
        float3 v    = gm_v[threadIndex];
        float3 f    = gm_f[threadIndex];
        float  im   = gm_inverseMasses[threadIndex];
        float  imdt = im * dt;

        // Swapping places for xp and x so that the x will contain the updated coordinates and xp - the
        // coordinates before update. This should be taken into account when (if) constraints are applied
        // after the update: x and xp have to be passed to constraints in the 'wrong' order.
        // TODO: Issue #3727
        gm_xp[threadIndex] = x;

        if (numTempScaleValues != NumTempScaleValues::None
            || parrinelloRahmanVelocityScaling != ParrinelloRahmanVelocityScaling::No)
        {
            float3 vp = v;

            if (numTempScaleValues != NumTempScaleValues::None)
            {
                float lambda = 1.0F;
                if (numTempScaleValues == NumTempScaleValues::Single)
                {
                    lambda = gm_lambdas[0];
                }
                else if (numTempScaleValues == NumTempScaleValues::Multiple)
                {
                    int tempScaleGroup = gm_tempScaleGroups[threadIndex];
                    lambda             = gm_lambdas[tempScaleGroup];
                }
                vp *= lambda;
            }

            if (parrinelloRahmanVelocityScaling == ParrinelloRahmanVelocityScaling::Diagonal)
            {
                vp.x -= prVelocityScalingMatrixDiagonal.x * v.x;
                vp.y -= prVelocityScalingMatrixDiagonal.y * v.y;
                vp.z -= prVelocityScalingMatrixDiagonal.z * v.z;
            }

            v = vp;
        }

        v += f * imdt;

        x += v * dt;
        gm_v[threadIndex] = v;
        gm_x[threadIndex] = x;
    }
}

/*! \brief Select templated kernel.
 *
 * Returns pointer to a CUDA kernel based on the number of temperature coupling groups and
 * whether or not the temperature and(or) pressure coupling is enabled.
 *
 * \param[in]  doTemperatureScaling   If the kernel with temperature coupling velocity scaling
 *                                    should be selected.
 * \param[in]  numTempScaleValues     Number of temperature coupling groups in the system.
 * \param[in]  parrinelloRahmanVelocityScaling  Type of the Parrinello-Rahman velocity scaling.
 *
 * \return                         Pointer to CUDA kernel
 */
inline auto selectLeapFrogKernelPtr(bool                            doTemperatureScaling,
                                    int                             numTempScaleValues,
                                    ParrinelloRahmanVelocityScaling parrinelloRahmanVelocityScaling)
{
    // Check input for consistency: if there is temperature coupling, at least one coupling group should be defined.
    GMX_ASSERT(!doTemperatureScaling || (numTempScaleValues > 0),
               "Temperature coupling was requested with no temperature coupling groups.");
    auto kernelPtr = leapfrog_kernel<NumTempScaleValues::None, ParrinelloRahmanVelocityScaling::No>;

    if (parrinelloRahmanVelocityScaling == ParrinelloRahmanVelocityScaling::No)
    {
        if (!doTemperatureScaling)
        {
            kernelPtr = leapfrog_kernel<NumTempScaleValues::None, ParrinelloRahmanVelocityScaling::No>;
        }
        else if (numTempScaleValues == 1)
        {
            kernelPtr = leapfrog_kernel<NumTempScaleValues::Single, ParrinelloRahmanVelocityScaling::No>;
        }
        else if (numTempScaleValues > 1)
        {
            kernelPtr = leapfrog_kernel<NumTempScaleValues::Multiple, ParrinelloRahmanVelocityScaling::No>;
        }
    }
    else if (parrinelloRahmanVelocityScaling == ParrinelloRahmanVelocityScaling::Diagonal)
    {
        if (!doTemperatureScaling)
        {
            kernelPtr = leapfrog_kernel<NumTempScaleValues::None, ParrinelloRahmanVelocityScaling::Diagonal>;
        }
        else if (numTempScaleValues == 1)
        {
            kernelPtr = leapfrog_kernel<NumTempScaleValues::Single, ParrinelloRahmanVelocityScaling::Diagonal>;
        }
        else if (numTempScaleValues > 1)
        {
            kernelPtr = leapfrog_kernel<NumTempScaleValues::Multiple, ParrinelloRahmanVelocityScaling::Diagonal>;
        }
    }
    else
    {
        GMX_RELEASE_ASSERT(false,
                           "Only isotropic Parrinello-Rahman pressure coupling is supported.");
    }
    return kernelPtr;
}


void launchLeapFrogKernel(const int                             numAtoms,
                          DeviceBuffer<Float3>                  d_x,
                          DeviceBuffer<Float3>                  d_xp,
                          DeviceBuffer<Float3>                  d_v,
                          const DeviceBuffer<Float3>            d_f,
                          const DeviceBuffer<float>             d_inverseMasses,
                          const float                           dt,
                          const bool                            doTemperatureScaling,
                          const int                             numTempScaleValues,
                          const DeviceBuffer<unsigned short>    d_tempScaleGroups,
                          const DeviceBuffer<float>             d_lambdas,
                          const ParrinelloRahmanVelocityScaling parrinelloRahmanVelocityScaling,
                          const Float3                          prVelocityScalingMatrixDiagonal,
                          const DeviceStream&                   deviceStream)
{
    // Checking the buffer types against the kernel argument types
    static_assert(sizeof(*d_inverseMasses) == sizeof(float), "Incompatible types");

    KernelLaunchConfig kernelLaunchConfig;

    kernelLaunchConfig.gridSize[0]      = (numAtoms + c_threadsPerBlock - 1) / c_threadsPerBlock;
    kernelLaunchConfig.blockSize[0]     = c_threadsPerBlock;
    kernelLaunchConfig.blockSize[1]     = 1;
    kernelLaunchConfig.blockSize[2]     = 1;
    kernelLaunchConfig.sharedMemorySize = 0;

    auto kernelPtr = selectLeapFrogKernelPtr(
            doTemperatureScaling, numTempScaleValues, parrinelloRahmanVelocityScaling);

    const auto kernelArgs = prepareGpuKernelArguments(kernelPtr,
                                                      kernelLaunchConfig,
                                                      &numAtoms,
                                                      asFloat3Pointer(&d_x),
                                                      asFloat3Pointer(&d_xp),
                                                      asFloat3Pointer(&d_v),
                                                      asFloat3Pointer(&d_f),
                                                      &d_inverseMasses,
                                                      &dt,
                                                      &d_lambdas,
                                                      &d_tempScaleGroups,
                                                      &prVelocityScalingMatrixDiagonal);
    launchGpuKernel(kernelPtr, kernelLaunchConfig, deviceStream, nullptr, "leapfrog_kernel", kernelArgs);
}

} // namespace gmx
