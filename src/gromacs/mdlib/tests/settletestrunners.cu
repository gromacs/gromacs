/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2019,2020, by the GROMACS development team, led by
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
 * \brief Defines the runner for CUDA version of SETTLE.
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 * \ingroup module_mdlib
 */
#include "gmxpre.h"

#include "settletestrunners.h"

#include "config.h"

#include <assert.h>

#include <cmath>

#include <algorithm>
#include <vector>

#include "gromacs/gpu_utils/devicebuffer.cuh"
#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/mdlib/settle_gpu.cuh"
#include "gromacs/utility/unique_cptr.h"

namespace gmx
{
namespace test
{

/*! \brief Apply SETTLE using GPU version of the algorithm.
 *
 * Initializes SETTLE object, copied data to the GPU, applies algorithm, copies the data back,
 * destroys the object. The coordinates, velocities and virial are updated in the testData object.
 *
 * \param[in,out] testData          An object, containing all the data structures needed by SETTLE.
 * \param[in]     pbc               Periodic boundary setup.
 * \param[in]     updateVelocities  If the velocities should be updated.
 * \param[in]     calcVirial        If the virial should be computed.
 * \param[in]     testDescription   Brief description that will be printed in case of test failure.
 */
void applySettleGpu(SettleTestData*  testData,
                    const t_pbc      pbc,
                    const bool       updateVelocities,
                    const bool       calcVirial,
                    gmx_unused const std::string& testDescription)
{
    // These should never fail since this function should only be called if CUDA is enabled and
    // there is a CUDA-capable device available.
    GMX_RELEASE_ASSERT(GMX_GPU == GMX_GPU_CUDA,
                       "CUDA version of SETTLE was called from non-CUDA build.");

    // TODO: Here we should check that at least 1 suitable GPU is available
    GMX_RELEASE_ASSERT(canPerformGpuDetection(), "Can't detect CUDA-capable GPUs.");

    DeviceInformation   deviceInfo;
    const DeviceContext deviceContext(deviceInfo);
    const DeviceStream  deviceStream(deviceContext, DeviceStreamPriority::Normal, false);

    auto settleGpu = std::make_unique<SettleGpu>(testData->mtop_, deviceContext, deviceStream);

    settleGpu->set(*testData->idef_, testData->mdatoms_);
    PbcAiuc pbcAiuc;
    setPbcAiuc(pbc.ndim_ePBC, pbc.box, &pbcAiuc);

    int numAtoms = testData->mdatoms_.homenr;

    float3 *d_x, *d_xp, *d_v;

    float3* h_x  = (float3*)(as_rvec_array(testData->x_.data()));
    float3* h_xp = (float3*)(as_rvec_array(testData->xPrime_.data()));
    float3* h_v  = (float3*)(as_rvec_array(testData->v_.data()));

    allocateDeviceBuffer(&d_x, numAtoms, deviceContext);
    allocateDeviceBuffer(&d_xp, numAtoms, deviceContext);
    allocateDeviceBuffer(&d_v, numAtoms, deviceContext);

    copyToDeviceBuffer(&d_x, (float3*)h_x, 0, numAtoms, deviceStream, GpuApiCallBehavior::Sync, nullptr);
    copyToDeviceBuffer(&d_xp, (float3*)h_xp, 0, numAtoms, deviceStream, GpuApiCallBehavior::Sync, nullptr);
    if (updateVelocities)
    {
        copyToDeviceBuffer(&d_v, (float3*)h_v, 0, numAtoms, deviceStream, GpuApiCallBehavior::Sync, nullptr);
    }
    settleGpu->apply(d_x, d_xp, updateVelocities, d_v, testData->reciprocalTimeStep_, calcVirial,
                     testData->virial_, pbcAiuc);

    copyFromDeviceBuffer((float3*)h_xp, &d_xp, 0, numAtoms, deviceStream, GpuApiCallBehavior::Sync, nullptr);
    if (updateVelocities)
    {
        copyFromDeviceBuffer((float3*)h_v, &d_v, 0, numAtoms, deviceStream,
                             GpuApiCallBehavior::Sync, nullptr);
    }

    freeDeviceBuffer(&d_x);
    freeDeviceBuffer(&d_xp);
    freeDeviceBuffer(&d_v);
}

} // namespace test
} // namespace gmx
