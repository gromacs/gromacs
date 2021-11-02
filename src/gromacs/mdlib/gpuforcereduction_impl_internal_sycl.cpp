/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2021, by the GROMACS development team, led by
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
 * \brief Implements GPU Force Reduction using SYCL
 *
 * \author Alan Gray <alang@nvidia.com>
 * \author Andrey Alekseenko <al42and@gmail.com>
 *
 * \ingroup module_mdlib
 */

#include "gmxpre.h"

#include "gpuforcereduction_impl_internal.h"

#include <utility>

#include "gromacs/gpu_utils/gmxsycl.h"
#include "gromacs/gpu_utils/devicebuffer.h"
#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/gpu_utils/gpueventsynchronizer.h"
#include "gromacs/utility/template_mp.h"

//! \brief Class name for reduction kernel
template<bool addRvecForce, bool accumulateForce>
class ReduceKernel;

namespace gmx
{

using cl::sycl::access::mode;

//! \brief Function returning the force reduction kernel lambda.
template<bool addRvecForce, bool accumulateForce>
static auto reduceKernel(cl::sycl::handler&                                 cgh,
                         DeviceAccessor<Float3, mode::read>                 a_nbnxmForce,
                         OptionalAccessor<Float3, mode::read, addRvecForce> a_rvecForceToAdd,
                         DeviceAccessor<Float3, accumulateForce ? mode::read_write : mode::write> a_forceTotal,
                         DeviceAccessor<int, cl::sycl::access::mode::read> a_cell,
                         const int                                         atomStart)
{
    a_nbnxmForce.bind(cgh);
    if constexpr (addRvecForce)
    {
        a_rvecForceToAdd.bind(cgh);
    }
    a_forceTotal.bind(cgh);
    a_cell.bind(cgh);

    return [=](cl::sycl::id<1> itemIdx) {
        // Set to nbnxnm force, then perhaps accumulate further to it
        Float3 temp = a_nbnxmForce[a_cell[itemIdx]];

        if constexpr (accumulateForce)
        {
            temp += a_forceTotal[itemIdx + atomStart];
        }

        if constexpr (addRvecForce)
        {
            temp += a_rvecForceToAdd[itemIdx + atomStart];
        }

        a_forceTotal[itemIdx + atomStart] = temp;
    };
}

//! \brief Force reduction SYCL kernel launch code.
template<bool addRvecForce, bool accumulateForce>
static void launchReductionKernel_(const int                   numAtoms,
                                   const int                   atomStart,
                                   const DeviceBuffer<Float3>& b_nbnxmForce,
                                   const DeviceBuffer<Float3>& b_rvecForceToAdd,
                                   DeviceBuffer<Float3>&       b_forceTotal,
                                   const DeviceBuffer<int>&    b_cell,
                                   const DeviceStream&         deviceStream)
{
    const cl::sycl::range<1> rangeNumAtoms(numAtoms);
    cl::sycl::queue          queue = deviceStream.stream();

    // We only need parts of b_rvecForceToAdd and b_forceTotal, so sub-buffers would be appropriate.
    // But hipSYCL does not support them yet, nor plans to. See Issue #4019.

    queue.submit([&](cl::sycl::handler& cgh) {
        auto kernel = reduceKernel<addRvecForce, accumulateForce>(
                cgh, b_nbnxmForce, b_rvecForceToAdd, b_forceTotal, b_cell, atomStart);
        cgh.parallel_for<ReduceKernel<addRvecForce, accumulateForce>>(rangeNumAtoms, kernel);
    });
}

/*! \brief Select templated Force reduction kernel and launch it. */
void launchForceReductionKernel(int                  numAtoms,
                                int                  atomStart,
                                bool                 addRvecForce,
                                bool                 accumulate,
                                DeviceBuffer<Float3> d_nbnxmForceToAdd,
                                DeviceBuffer<Float3> d_rvecForceToAdd,
                                DeviceBuffer<Float3> d_baseForce,
                                DeviceBuffer<int>    d_cell,
                                const DeviceStream&  deviceStream)
{
    dispatchTemplatedFunction(
            [&](auto addRvecForce_, auto accumulateForce_) {
                return launchReductionKernel_<addRvecForce_, accumulateForce_>(
                        numAtoms, atomStart, d_nbnxmForceToAdd, d_rvecForceToAdd, d_baseForce, d_cell, deviceStream);
            },
            addRvecForce,
            accumulate);
}

} // namespace gmx
