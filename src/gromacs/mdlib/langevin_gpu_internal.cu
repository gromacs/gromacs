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
 * \brief Implements Langevin (SD) integrator using CUDA
 *
 * This file contains implementation of the Langevin (SD) integrator
 * using CUDA, including class initialization, data-structures management
 * and GPU kernel.
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 * \author Magnus Lundborg <lundborg.magnus@gmail.com>
 *
 * \ingroup module_mdlib
 */
#include "gmxpre.h"

#include "langevin_gpu_internal.h"

#include "gromacs/gpu_utils/cudautils.cuh"
#include "gromacs/gpu_utils/devicebuffer.h"
#include "gromacs/gpu_utils/typecasts_cuda_hip.h"
#include "gromacs/pbcutil/pbc_aiuc_cuda.cuh"
#include "gromacs/random/tabulatednormaldistribution_gpu.h"
#include "gromacs/random/threefry.h"

#include "langevin_gpu.h"


namespace gmx
{

/*!\brief Number of CUDA threads in a block
 *
 * \todo Check if using smaller block size will lead to better performance.
 */
constexpr static int c_threadsPerBlock = 256;
//! Maximum number of threads in a block (for __launch_bounds__)
constexpr static int c_maxThreadsPerBlock = c_threadsPerBlock;

/*! \brief Main kernel for Langevin integrator.
 *
 *  The coordinates and velocities are updated on the GPU. Also saves the intermediate values of the coordinates for
 *  further use in constraints.
 *
 *  Each GPU thread works with a single particle. Empty declaration is needed to
 *  avoid "no previous prototype for function" clang warning.
 *
 * \tparam        updateType           Langevin integrator update type. The integration is divided in
 *                                     three steps if there are constraints:
 *                                     SDUpdate::ForcesOnly, Constraints, SDUpdate::FrictionAndNoiseOnly.
 *                                     If there are no constraints, SDUpdate::Combined can be used to do
 *                                     the full update in a single kernel launch.
 * \param[in]     numAtoms             Total number of atoms.
 * \param[in,out] gm_x                 Coordinates to update upon integration.
 * \param[out]    gm_xp                A copy of the coordinates before the integration (for constraints).
 * \param[in,out] gm_v                 Velocities to update.
 * \param[in]     gm_f                 Atomic forces.
 * \param[in]     gm_inverseMasses     Reciprocal masses.
 * \param[in]     dt                   Timestep.
 * \param[in]     seed                 Random seed for sd integrator.
 * \param[in]     step                 The step number in the simulation.
 * \param[in]     gm_tempCouplGroups   Mapping of atoms into temperate coupling groups.
 * \param[in]     gm_sdSigmaV          The sigma of the stochastic dynamics noise.
 * \param[in]     gm_sdConstEm         EM (alpha in the SD equation). We lose precision (compared to the CPU version) by using float.
 * \param[in]     gm_distributionTable The normal distribution table.
 */
template<SDUpdate updateType>
__launch_bounds__(c_maxThreadsPerBlock) __global__
        void langevin_kernel(const int numAtoms,
                             float3* __restrict__ gm_x,
                             float3* __restrict__ gm_xp,
                             float3* __restrict__ gm_v,
                             const float3* __restrict__ gm_f,
                             const float* __restrict__ gm_inverseMasses,
                             const float dt,
                             const int   seed,
                             const int   step,
                             const unsigned short* __restrict__ gm_tempCouplGroups,
                             const float* __restrict__ gm_sdSigmaV,
                             const float* __restrict__ gm_sdConstEm,
                             const float* __restrict__ gm_distributionTable)
{
    int threadIndex = blockIdx.x * blockDim.x + threadIdx.x;
    if (threadIndex < numAtoms)
    {
        // Even 0 bits internal counter gives 2x64 ints (more than enough for three table lookups)
        gmx::ThreeFry2x64<0> rng(seed, gmx::RandomDomain::UpdateCoordinates);
        gmx::TabulatedNormalDistributionGpu<sc_normalDistributionTableBits> dist(gm_distributionTable);
        rng.restart(step, threadIndex);
        dist.reset();

        float3 distByDim;
        if (updateType != SDUpdate::ForcesOnly)
        {
            distByDim.x = dist(rng);
            distByDim.y = dist(rng);
            distByDim.z = dist(rng);
        }

        float3 x                     = gm_x[threadIndex];
        float3 v                     = gm_v[threadIndex];
        float3 f                     = gm_f[threadIndex];
        float  inverseMass           = gm_inverseMasses[threadIndex];
        float  inverseSqrtMass       = sqrt(inverseMass);
        float  inverseMassDt         = inverseMass * dt;
        int    temperatureCouplGroup = gm_tempCouplGroups[threadIndex];
        float  sdSigmaV              = gm_sdSigmaV[temperatureCouplGroup];
        float  sdConstEm             = gm_sdConstEm[temperatureCouplGroup];


        if (updateType != SDUpdate::FrictionAndNoiseOnly)
        {
            // Swapping places for xp and x so that the x will contain the updated coordinates and xp - the
            // coordinates before update. This should be taken into account when (if) constraints are applied
            // after the update: x and xp have to be passed to constraints in the 'wrong' order.
            // TODO: Issue #3727
            gm_xp[threadIndex] = x;
        }

        if (updateType == SDUpdate::ForcesOnly)
        {
            float3 vn = v + f * inverseMassDt;
            v         = vn;
            x += v * dt;
        }
        else if (updateType == SDUpdate::FrictionAndNoiseOnly)
        {
            float3 vn = v;
            v         = vn * sdConstEm + inverseSqrtMass * sdSigmaV * distByDim;
            // The previous phase already updated the
            // positions with a full v*dt term that must
            // now be half removed.
            x += 0.5 * (v - vn) * dt;
        }
        else
        {
            // SDUpdate::Combined: forces, friction and noise in a single step.
            float3 vn = v + f * inverseMassDt;
            v         = vn * sdConstEm + inverseSqrtMass * sdSigmaV * distByDim;
            x += 0.5 * (vn + v) * dt;
        }
        gm_v[threadIndex] = v;
        gm_x[threadIndex] = x;
    }
    return;
}

/*! \brief Select templated kernel.
 *
 * Returns pointer to a CUDA kernel based on the type of SD integration.
 *
 * \param[in]  updateType   Langevin integrator update type. The integration is divided in
 *                          three steps if there are constraints:
 *                          SDUpdate::ForcesOnly, Constraints, SDUpdate::FrictionAndNoiseOnly
 *                          If there are no constraints it is only SDUpdate::Combined.
 *
 * \return                  Pointer to CUDA kernel
 */
inline auto selectLangevinKernelPtr(const SDUpdate updateType)
{
    GMX_ASSERT(updateType == SDUpdate::ForcesOnly || updateType == SDUpdate::FrictionAndNoiseOnly
                       || updateType == SDUpdate::Combined,
               "Unknown SD integrator update type.");

    auto kernelPtr = langevin_kernel<SDUpdate::ForcesOnly>;

    if (updateType == SDUpdate::FrictionAndNoiseOnly)
    {
        kernelPtr = langevin_kernel<SDUpdate::FrictionAndNoiseOnly>;
    }
    else if (updateType == SDUpdate::Combined)
    {
        kernelPtr = langevin_kernel<SDUpdate::Combined>;
    }
    return kernelPtr;
}


void launchLangevinKernel(const int                          numAtoms,
                          DeviceBuffer<Float3>               d_x,
                          DeviceBuffer<Float3>               d_xp,
                          DeviceBuffer<Float3>               d_v,
                          const DeviceBuffer<Float3>         d_f,
                          const DeviceBuffer<float>          d_inverseMasses,
                          const float                        dt,
                          const int                          seed,
                          const int                          step,
                          const DeviceBuffer<unsigned short> d_tempCouplGroups,
                          const DeviceBuffer<float>          d_sdSigmaV,
                          const DeviceBuffer<float>          d_sdConstEm,
                          const DeviceBuffer<float>          d_distributionTable,
                          const SDUpdate                     updateType,
                          const DeviceStream&                deviceStream)
{
    // Checking the buffer types against the kernel argument types
    static_assert(sizeof(*d_inverseMasses) == sizeof(float), "Incompatible types");

    KernelLaunchConfig kernelLaunchConfig;

    kernelLaunchConfig.gridSize[0]      = (numAtoms + c_threadsPerBlock - 1) / c_threadsPerBlock;
    kernelLaunchConfig.blockSize[0]     = c_threadsPerBlock;
    kernelLaunchConfig.blockSize[1]     = 1;
    kernelLaunchConfig.blockSize[2]     = 1;
    kernelLaunchConfig.sharedMemorySize = 0;

    auto kernelPtr = selectLangevinKernelPtr(updateType);

    const auto kernelArgs = prepareGpuKernelArguments(kernelPtr,
                                                      kernelLaunchConfig,
                                                      &numAtoms,
                                                      asFloat3Pointer(&d_x),
                                                      asFloat3Pointer(&d_xp),
                                                      asFloat3Pointer(&d_v),
                                                      asFloat3Pointer(&d_f),
                                                      &d_inverseMasses,
                                                      &dt,
                                                      &seed,
                                                      &step,
                                                      &d_tempCouplGroups,
                                                      &d_sdSigmaV,
                                                      &d_sdConstEm,
                                                      &d_distributionTable);
    launchGpuKernel(kernelPtr, kernelLaunchConfig, deviceStream, nullptr, "langevin_kernel", kernelArgs);
}

} // namespace gmx
