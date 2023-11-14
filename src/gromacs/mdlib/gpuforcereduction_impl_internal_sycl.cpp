/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2021- The GROMACS Authors
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
 * \brief Implements GPU Force Reduction using SYCL
 *
 * \author Alan Gray <alang@nvidia.com>
 * \author Andrey Alekseenko <al42and@gmail.com>
 *
 * \ingroup module_mdlib
 */

#include "gmxpre.h"

#include <utility>

#include "gromacs/gpu_utils/devicebuffer.h"
#include "gromacs/gpu_utils/gmxsycl.h"
#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/gpu_utils/gpueventsynchronizer.h"
#include "gromacs/utility/template_mp.h"

#include "gpuforcereduction_impl_internal.h"

//! \brief Class name for reduction kernel
template<bool addRvecForce, bool accumulateForce>
class ReduceKernel;

namespace gmx
{

using mode = sycl::access_mode;

//! \brief Function returning the force reduction kernel lambda.
template<bool addRvecForce, bool accumulateForce>
static auto reduceKernel(const Float3* __restrict__ gm_nbnxmForce,
                         const Float3* __restrict__ gm_rvecForceToAdd /* used iff addRvecForce */,
                         Float3* __restrict__ gm_forceTotal,
                         const int* __restrict__ gm_cell,
                         const int atomStart)
{
    return [=](sycl::id<1> itemIdx) {
        // Set to nbnxnm force, then perhaps accumulate further to it
        Float3 temp = gm_nbnxmForce[gm_cell[itemIdx]];

        if constexpr (accumulateForce)
        {
            temp += gm_forceTotal[itemIdx + atomStart];
        }

        if constexpr (addRvecForce)
        {
            temp += gm_rvecForceToAdd[itemIdx + atomStart];
        }

        gm_forceTotal[itemIdx + atomStart] = temp;
    };
}

//! \brief Force reduction SYCL kernel launch code.
template<bool addRvecForce, bool accumulateForce>
static void launchReductionKernel_(const int                   numAtoms,
                                   const int                   atomStart,
                                   const DeviceBuffer<Float3>& d_nbnxmForce,
                                   const DeviceBuffer<Float3>& d_rvecForceToAdd,
                                   DeviceBuffer<Float3>&       d_forceTotal,
                                   const DeviceBuffer<int>&    d_cell,
                                   const DeviceStream&         deviceStream)
{
    const sycl::range<1> rangeNumAtoms(numAtoms);
    sycl::queue          queue = deviceStream.stream();

    queue.submit(GMX_SYCL_DISCARD_EVENT[&](sycl::handler & cgh) {
        auto kernel = reduceKernel<addRvecForce, accumulateForce>(d_nbnxmForce.get_pointer(),
                                                                  d_rvecForceToAdd.get_pointer(),
                                                                  d_forceTotal.get_pointer(),
                                                                  d_cell.get_pointer(),
                                                                  atomStart);
        cgh.parallel_for<ReduceKernel<addRvecForce, accumulateForce>>(rangeNumAtoms, kernel);
    });
}

/*! \brief Select templated Force reduction kernel and launch it. */
void launchForceReductionKernel(int                    numAtoms,
                                int                    atomStart,
                                bool                   addRvecForce,
                                bool                   accumulate,
                                DeviceBuffer<Float3>   d_nbnxmForceToAdd,
                                DeviceBuffer<Float3>   d_rvecForceToAdd,
                                DeviceBuffer<Float3>   d_baseForce,
                                DeviceBuffer<int>      d_cell,
                                const DeviceStream&    deviceStream,
                                DeviceBuffer<uint64_t> d_forcesReadyNvshmemFlags,
                                const uint64_t         forcesReadyNvshmemFlagsCounter)
{
    GMX_UNUSED_VALUE(d_forcesReadyNvshmemFlags);
    GMX_UNUSED_VALUE(forcesReadyNvshmemFlagsCounter);
    dispatchTemplatedFunction(
            [&](auto addRvecForce_, auto accumulateForce_) {
                return launchReductionKernel_<addRvecForce_, accumulateForce_>(
                        numAtoms, atomStart, d_nbnxmForceToAdd, d_rvecForceToAdd, d_baseForce, d_cell, deviceStream);
            },
            addRvecForce,
            accumulate);
}

} // namespace gmx
