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
 * \brief Implements GPU halo exchange using SYCL.
 *
 *
 * \author Alan Gray <alang@nvidia.com>
 * \author Andrey Alekseenko <al42and@gmail.com>
 *
 * \ingroup module_domdec
 */
#include "gmxpre.h"

#include "gromacs/gpu_utils/devicebuffer_sycl.h"
#include "gromacs/gpu_utils/gmxsycl.h"
#include "gromacs/gpu_utils/gputraits_sycl.h"

#include "domdec_struct.h"
#include "gpuhaloexchange_impl_gpu.h"

template<bool usePbc>
class PackSendBufKernel;

template<bool accumulate>
class UnpackRecvBufKernel;

namespace gmx
{

template<bool usePbc>
auto packSendBufKernel(sycl::handler&                                        cgh,
                       DeviceAccessor<Float3, sycl::access_mode::read_write> a_dataPacked,
                       DeviceAccessor<Float3, sycl::access_mode::read>       a_data,
                       DeviceAccessor<int, sycl::access_mode::read>          a_map,
                       int                                                   mapSize,
                       Float3                                                coordinateShift)
{
    a_dataPacked.bind(cgh);
    a_data.bind(cgh);
    a_map.bind(cgh);

    return [=](sycl::id<1> itemIdx) {
        const int threadIndex = itemIdx;
        if (threadIndex < mapSize)
        {
            if constexpr (usePbc)
            {
                a_dataPacked[itemIdx] = a_data[a_map[itemIdx]] + coordinateShift;
            }
            else
            {
                a_dataPacked[itemIdx] = a_data[a_map[itemIdx]];
            }
        }
    };
}

/*! \brief unpack non-local force data buffer on the GPU using pre-populated "map" containing index
 * information.
 * \param[in]  cgh           SYCL's Command group handler
 * \param[out] a_data        full array of force values
 * \param[in]  a_dataPacked  packed array of force values to be transferred
 * \param[in]  a_map         array of indices defining mapping from full to packed array
 * \param[in]  mapSize       number of elements in map array
 */
template<bool accumulate>
auto unpackRecvBufKernel(sycl::handler&                                        cgh,
                         DeviceAccessor<Float3, sycl::access_mode::read_write> a_data,
                         DeviceAccessor<Float3, sycl::access_mode::read>       a_dataPacked,
                         DeviceAccessor<int, sycl::access_mode::read>          a_map,
                         int                                                   mapSize)
{
    a_dataPacked.bind(cgh);
    a_data.bind(cgh);
    a_map.bind(cgh);

    return [=](sycl::id<1> itemIdx) {
        const int threadIndex = itemIdx;
        if (threadIndex < mapSize)
        {
            if constexpr (accumulate)
            {
                a_data[a_map[itemIdx]] += a_dataPacked[itemIdx];
            }
            else
            {
                a_data[a_map[itemIdx]] = a_dataPacked[itemIdx];
            }
        }
    };
}


template<bool usePbc, class... Args>
static void launchPackSendBufKernel(const DeviceStream& deviceStream, int xSendSize, Args&&... args)
{
    using kernelNameType = PackSendBufKernel<usePbc>;

    const sycl::range<1> range(xSendSize);
    sycl::queue          q = deviceStream.stream();

    q.submit(GMX_SYCL_DISCARD_EVENT[&](sycl::handler & cgh) {
        auto kernel = packSendBufKernel<usePbc>(cgh, std::forward<Args>(args)...);
        cgh.parallel_for<kernelNameType>(range, kernel);
    });
}

template<bool accumulateForces, class... Args>
static void launchUnpackRecvBufKernel(const DeviceStream& deviceStream, int fRecvSize, Args&&... args)
{
    using kernelNameType = UnpackRecvBufKernel<accumulateForces>;

    const sycl::range<1> range(fRecvSize);
    sycl::queue          q = deviceStream.stream();

    q.submit(GMX_SYCL_DISCARD_EVENT[&](sycl::handler & cgh) {
        auto kernel = unpackRecvBufKernel<accumulateForces>(cgh, std::forward<Args>(args)...);
        cgh.parallel_for<kernelNameType>(range, kernel);
    });
}

void GpuHaloExchange::Impl::launchPackXKernel(const matrix box)
{
    const int size = xSendSize_;
    // The coordinateShift changes between steps when we have
    // performed a DD partition, or have updated the box e.g. when
    // performing pressure coupling. So, for simplicity, the box
    // is used every step to pass the shift vector as an argument of
    // the packing kernel, even when PBC is not in use.
    const int    boxDimensionIndex = dd_->dim[dimIndex_];
    const Float3 coordinateShift{ box[boxDimensionIndex][XX],
                                  box[boxDimensionIndex][YY],
                                  box[boxDimensionIndex][ZZ] };

    // Avoid launching kernel when there is no work to do
    if (size > 0)
    {
        if (usePBC_)
        {
            launchPackSendBufKernel<true>(
                    *haloStream_, size, d_sendBuf_, d_x_, d_indexMap_, size, coordinateShift);
        }
        else
        {
            launchPackSendBufKernel<false>(
                    *haloStream_, size, d_sendBuf_, d_x_, d_indexMap_, size, coordinateShift);
        }
    }
}

// The following method should be called after non-local buffer operations,
// and before the local buffer operations.
void GpuHaloExchange::Impl::launchUnpackFKernel(bool accumulateForces)
{
    const int size = fRecvSize_;
    if (size > 0)
    {
        if (accumulateForces)
        {
            launchUnpackRecvBufKernel<true>(*haloStream_, size, d_f_, d_recvBuf_, d_indexMap_, size);
        }
        else
        {
            launchUnpackRecvBufKernel<false>(*haloStream_, size, d_f_, d_recvBuf_, d_indexMap_, size);
        }
    }
}

} // namespace gmx
