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
 * \brief Implements Langevin (SD) integrator using SYCL
 *
 * This file contains SYCL implementation of the back-end specific kernel
 * launch for the Langevin (SD) integrator. The actual kernel code is shared
 * with the CUDA/HIP backends through the portable random number facilities
 * in gmx::ThreeFry2x64 and gmx::TabulatedNormalDistributionGpu.
 *
 * \author Magnus Lundborg <lundborg.magnus@gmail.com>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 *
 * \ingroup module_mdlib
 */
#include "gmxpre.h"

#include "gromacs/gpu_utils/devicebuffer.h"
#include "gromacs/gpu_utils/gmxsycl.h"
#include "gromacs/random/seed.h"
#include "gromacs/random/tabulatednormaldistribution_gpu.h"
#include "gromacs/random/threefry.h"
#include "gromacs/utility/template_mp.h"

#include "langevin_gpu.h"
#include "langevin_gpu_internal.h"

//! \brief Class name for the Langevin SYCL kernel (used by some SYCL implementations).
template<gmx::SDUpdate updateType>
class LangevinKernel;

namespace gmx
{

/*! \brief Main kernel for the Langevin (SD) integrator.
 *
 * Mirrors the corresponding CUDA/HIP implementation; see langevin_gpu_internal.cu
 * for a detailed description of the algorithm and arguments.
 */
template<SDUpdate updateType>
auto langevinKernel(Float3* __restrict__ gm_x,
                    Float3* __restrict__ gm_xp,
                    Float3* __restrict__ gm_v,
                    const Float3* __restrict__ gm_f,
                    const float* __restrict__ gm_inverseMasses,
                    float dt,
                    int   seed,
                    int   step,
                    const unsigned short* __restrict__ gm_tempCouplGroups,
                    const float* __restrict__ gm_sdSigmaV,
                    const float* __restrict__ gm_sdConstEm,
                    const float* __restrict__ gm_distributionTable)
{
    return [=](sycl::id<1> itemIdx)
    {
        const int threadIndex = itemIdx;

        // Even 0 bits internal counter gives 2x64 ints (more than enough for three table lookups)
        gmx::ThreeFry2x64<0> rng(seed, gmx::RandomDomain::UpdateCoordinates);
        gmx::TabulatedNormalDistributionGpu<sc_normalDistributionTableBits> dist(gm_distributionTable);
        rng.restart(step, threadIndex);
        dist.reset();

        Float3 distByDim{ 0, 0, 0 };
        if constexpr (updateType != SDUpdate::ForcesOnly)
        {
            distByDim[0] = dist(rng);
            distByDim[1] = dist(rng);
            distByDim[2] = dist(rng);
        }

        Float3       x                     = gm_x[threadIndex];
        Float3       v                     = gm_v[threadIndex];
        const Float3 f                     = gm_f[threadIndex];
        const float  inverseMass           = gm_inverseMasses[threadIndex];
        const float  inverseSqrtMass       = sycl::sqrt(inverseMass);
        const float  inverseMassDt         = inverseMass * dt;
        const int    temperatureCouplGroup = gm_tempCouplGroups[threadIndex];
        const float  sdSigmaV              = gm_sdSigmaV[temperatureCouplGroup];
        const float  sdConstEm             = gm_sdConstEm[temperatureCouplGroup];

        if constexpr (updateType != SDUpdate::FrictionAndNoiseOnly)
        {
            // Save initial coordinates so constraints can use them.
            gm_xp[threadIndex] = x;
        }

        if constexpr (updateType == SDUpdate::ForcesOnly)
        {
            const Float3 vn = v + f * inverseMassDt;
            v               = vn;
            x += v * dt;
        }
        else if constexpr (updateType == SDUpdate::FrictionAndNoiseOnly)
        {
            const Float3 vn = v;
            v               = vn * sdConstEm + distByDim * (inverseSqrtMass * sdSigmaV);
            // The previous phase already updated the positions with a full
            // v*dt term that must now be half removed.
            x += (v - vn) * (0.5F * dt);
        }
        else
        {
            // SDUpdate::Combined: forces, friction and noise in a single step.
            const Float3 vn = v + f * inverseMassDt;
            v               = vn * sdConstEm + distByDim * (inverseSqrtMass * sdSigmaV);
            x += (vn + v) * (0.5F * dt);
        }
        gm_v[threadIndex] = v;
        gm_x[threadIndex] = x;
    };
}

//! \brief Langevin SYCL kernel launch code.
template<SDUpdate updateType, class... Args>
static void launchLangevinKernel(const DeviceStream& deviceStream, int numAtoms, Args&&... args)
{
    using kernelNameType = LangevinKernel<updateType>;

    const sycl::range<1> rangeAllAtoms(numAtoms);
    auto                 kernelFunctionBuilder = langevinKernel<updateType>;
    syclSubmitWithoutCghOrEvent<kernelNameType>(
            deviceStream.stream(), kernelFunctionBuilder, rangeAllAtoms, std::forward<Args>(args)...);
}

/*! \brief Select templated kernel and launch it. */
template<class... Args>
static inline void launchLangevinKernel(SDUpdate updateType, Args&&... args)
{
    GMX_ASSERT(updateType == SDUpdate::ForcesOnly || updateType == SDUpdate::FrictionAndNoiseOnly
                       || updateType == SDUpdate::Combined,
               "Unknown SD integrator update type.");

    dispatchTemplatedFunction(
            [&](auto updateType_)
            { return launchLangevinKernel<updateType_>(std::forward<Args>(args)...); },
            updateType);
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
    launchLangevinKernel(updateType,
                         deviceStream,
                         numAtoms,
                         d_x.get_pointer(),
                         d_xp.get_pointer(),
                         d_v.get_pointer(),
                         d_f.get_pointer(),
                         d_inverseMasses.get_pointer(),
                         dt,
                         seed,
                         step,
                         d_tempCouplGroups.get_pointer(),
                         d_sdSigmaV.get_pointer(),
                         d_sdConstEm.get_pointer(),
                         d_distributionTable.get_pointer());
}

} // namespace gmx
