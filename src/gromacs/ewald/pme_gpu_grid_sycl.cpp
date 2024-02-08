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
 * \brief Implements PME GPU halo exchange and PME GPU - Host FFT grid conversion
 * functions. These functions are used for PME decomposition in mixed-mode
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 *
 * \ingroup module_ewald
 */

#include "gmxpre.h"

#include "config.h"

#include <cstdlib>

#include "gromacs/fft/parallel_3dfft.h"
#include "gromacs/gpu_utils/gmxsycl.h"
#include "gromacs/gpu_utils/syclutils.h"
#include "gromacs/math/vec.h"
#include "gromacs/timing/wallcycle.h"

#include "pme_gpu_grid.h"
#include "pme_gpu_types.h"
#include "pme_gpu_types_host.h"
#include "pme_gpu_types_host_impl.h"

using mode = sycl::access_mode;

/*! \brief Sub-group size for conversion kernels
 *
 * Chosen to match relevant hardware widths on supported hardware.  In
 * particular, PonteVecchio is 32-wide. */
static constexpr int sc_subGroupSize = 32;

//! Handles a kernel which packs non-contiguous overlap data in all 8 neighboring directions
class PackHaloExternal
{
public:
    /*! \brief
     * Returns the pack kernel
     *
     * \param[in] myGridX,myGridY            Local domain size in X and Y dimension
     * \param[in] pmeSize                    Local PME grid size
     * \param[in] gm_realGrid                PME device grid
     * \param[out] gm_transferGridUp         Device array used to pack data to go up
     * \param[out] gm_transferGridDown       Device array used to pack data to go down
     * \param[out] gm_transferGridLeft       Device array used to pack data to go left
     * \param[out] gm_transferGridRight      Device array used to pack data to go right
     * \param[out] gm_transferGridUpLeft     Device array used to pack data to go up+left
     * \param[out] gm_transferGridDownLeft   Device array used to pack data to go down+left
     * \param[out] gm_transferGridUpRight    Device array used to pack data to go up+right
     * \param[out] gm_transferGridDownRight  Device array used to pack data to go down+right
     * \param[in] overlapSizeUp,overlapSizeDown,overlapSizeLeft,overlapSizeRight
     *                                       Halo size in 4 directions
     * \tparam    subGroupSize               Size of the sub-group.
     */
    template<int subGroupSize>
    static auto kernel(size_t      myGridX,
                       size_t      myGridY,
                       sycl::uint3 pmeSize,
                       const float* __restrict__ gm_realGrid,
                       float* __restrict__ gm_transferGridUp,
                       float* __restrict__ gm_transferGridDown,
                       float* __restrict__ gm_transferGridLeft,
                       float* __restrict__ gm_transferGridRight,
                       float* __restrict__ gm_transferGridUpLeft,
                       float* __restrict__ gm_transferGridDownLeft,
                       float* __restrict__ gm_transferGridUpRight,
                       float* __restrict__ gm_transferGridDownRight,
                       const size_t overlapSizeUp,
                       const size_t overlapSizeDown,
                       const size_t overlapSizeLeft,
                       const size_t overlapSizeRight)
    {
        return [=](sycl::nd_item<3> item_ct1) [[intel::reqd_sub_group_size(subGroupSize)]]
        {
            size_t iz = item_ct1.get_local_id(2) + item_ct1.get_group(2) * item_ct1.get_local_range(2);
            size_t iy = item_ct1.get_local_id(1) + item_ct1.get_group(1) * item_ct1.get_local_range(1);
            size_t ix = item_ct1.get_local_id(0) + item_ct1.get_group(0) * item_ct1.get_local_range(0);

            // we might get iz greater than pmeSize.z when pmeSize.z is not multiple of
            // threadsAlongZDim(see below), same for iy when it's not multiple of threadsAlongYDim
            if (iz >= pmeSize.z() || iy >= myGridY)
            {
                return;
            }

            // up
            if (ix < overlapSizeUp)
            {
                size_t pmeIndex = (ix + pmeSize.x() - overlapSizeUp) * pmeSize.y() * pmeSize.z()
                                  + iy * pmeSize.z() + iz;
                size_t packedIndex             = ix * myGridY * pmeSize.z() + iy * pmeSize.z() + iz;
                gm_transferGridUp[packedIndex] = gm_realGrid[pmeIndex];
            }

            // down
            if (ix >= myGridX - overlapSizeDown)
            {
                size_t pmeIndex =
                        (ix + overlapSizeDown) * pmeSize.y() * pmeSize.z() + iy * pmeSize.z() + iz;
                size_t packedIndex = (ix - (myGridX - overlapSizeDown)) * myGridY * pmeSize.z()
                                     + iy * pmeSize.z() + iz;
                gm_transferGridDown[packedIndex] = gm_realGrid[pmeIndex];
            }

            // left
            if (iy < overlapSizeLeft)
            {
                size_t pmeIndex = ix * pmeSize.y() * pmeSize.z()
                                  + (iy + pmeSize.y() - overlapSizeLeft) * pmeSize.z() + iz;
                size_t packedIndex = ix * overlapSizeLeft * pmeSize.z() + iy * pmeSize.z() + iz;
                gm_transferGridLeft[packedIndex] = gm_realGrid[pmeIndex];
            }

            // right
            if (iy >= myGridY - overlapSizeRight)
            {
                size_t pmeIndex =
                        ix * pmeSize.y() * pmeSize.z() + (iy + overlapSizeRight) * pmeSize.z() + iz;
                size_t packedIndex = ix * overlapSizeRight * pmeSize.z()
                                     + (iy - (myGridY - overlapSizeRight)) * pmeSize.z() + iz;
                gm_transferGridRight[packedIndex] = gm_realGrid[pmeIndex];
            }

            // up left
            if (ix < overlapSizeUp && iy < overlapSizeLeft)
            {
                size_t pmeIndex = (ix + pmeSize.x() - overlapSizeUp) * pmeSize.y() * pmeSize.z()
                                  + (iy + pmeSize.y() - overlapSizeLeft) * pmeSize.z() + iz;
                size_t packedIndex = ix * overlapSizeLeft * pmeSize.z() + iy * pmeSize.z() + iz;
                gm_transferGridUpLeft[packedIndex] = gm_realGrid[pmeIndex];
            }

            // down left
            if (ix >= myGridX - overlapSizeDown && iy < overlapSizeLeft)
            {
                size_t pmeIndex = (ix + overlapSizeDown) * pmeSize.y() * pmeSize.z()
                                  + (iy + pmeSize.y() - overlapSizeLeft) * pmeSize.z() + iz;
                size_t packedIndex = (ix - (myGridX - overlapSizeDown)) * overlapSizeLeft * pmeSize.z()
                                     + iy * pmeSize.z() + iz;
                gm_transferGridDownLeft[packedIndex] = gm_realGrid[pmeIndex];
            }

            // up right
            if (ix < overlapSizeUp && iy >= myGridY - overlapSizeRight)
            {
                size_t pmeIndex = (ix + pmeSize.x() - overlapSizeUp) * pmeSize.y() * pmeSize.z()
                                  + (iy + overlapSizeRight) * pmeSize.z() + iz;
                size_t packedIndex = ix * overlapSizeRight * pmeSize.z()
                                     + (iy - (myGridY - overlapSizeRight)) * pmeSize.z() + iz;
                gm_transferGridUpRight[packedIndex] = gm_realGrid[pmeIndex];
            }

            // down right
            if (ix >= myGridX - overlapSizeDown && iy >= myGridY - overlapSizeRight)
            {
                size_t pmeIndex = (ix + overlapSizeDown) * pmeSize.y() * pmeSize.z()
                                  + (iy + overlapSizeRight) * pmeSize.z() + iz;
                size_t packedIndex = (ix - (myGridX - overlapSizeDown)) * overlapSizeRight * pmeSize.z()
                                     + (iy - (myGridY - overlapSizeRight)) * pmeSize.z() + iz;
                gm_transferGridDownRight[packedIndex] = gm_realGrid[pmeIndex];
            }
        };
    }
};

//! Handles a kernel which gathers data from halo region in all 8 neighboring directions
class UnpackHaloExternal
{
public:
    /*! \brief
     * Returns the unpack kernel
     *
     * \param[in] myGridX,myGridY            Local domain size in X and Y dimension
     * \param[in] pmeSize                    Local PME grid size
     * \param[in] gm_realGrid                PME device grid
     * \param[out] gm_transferGridUp         Device array used to pack data to go up
     * \param[out] gm_transferGridDown       Device array used to pack data to go down
     * \param[out] gm_transferGridLeft       Device array used to pack data to go left
     * \param[out] gm_transferGridRight      Device array used to pack data to go right
     * \param[out] gm_transferGridUpLeft     Device array used to pack data to go up+left
     * \param[out] gm_transferGridDownLeft   Device array used to pack data to go down+left
     * \param[out] gm_transferGridUpRight    Device array used to pack data to go up+right
     * \param[out] gm_transferGridDownRight  Device array used to pack data to go down+right
     * \param[in] overlapSizeUp,overlapSizeDown,overlapSizeLeft,overlapSizeRight
     *                                       Halo size in 4 directions
     * \tparam    subGroupSize               Size of the sub-group.
     */
    template<int subGroupSize>
    static auto kernel(size_t      myGridX,
                       size_t      myGridY,
                       sycl::uint3 pmeSize,
                       float* __restrict__ gm_realGrid,
                       const float* __restrict__ gm_transferGridUp,
                       const float* __restrict__ gm_transferGridDown,
                       const float* __restrict__ gm_transferGridLeft,
                       const float* __restrict__ gm_transferGridRight,
                       const float* __restrict__ gm_transferGridUpLeft,
                       const float* __restrict__ gm_transferGridDownLeft,
                       const float* __restrict__ gm_transferGridUpRight,
                       const float* __restrict__ gm_transferGridDownRight,
                       size_t overlapSizeUp,
                       size_t overlapSizeDown,
                       size_t overlapSizeLeft,
                       size_t overlapSizeRight)
    {
        return [=](sycl::nd_item<3> item_ct1) [[intel::reqd_sub_group_size(subGroupSize)]]
        {
            size_t iz = item_ct1.get_local_id(2) + item_ct1.get_group(2) * item_ct1.get_local_range(2);
            size_t iy = item_ct1.get_local_id(1) + item_ct1.get_group(1) * item_ct1.get_local_range(1);
            size_t ix = item_ct1.get_local_id(0) + item_ct1.get_group(0) * item_ct1.get_local_range(0);

            // we might get iz greater than pmeSize.z when pmeSize.z is not multiple of
            // threadsAlongZDim(see below), same for iy when it's not multiple of threadsAlongYDim
            if (iz >= pmeSize.z() || iy >= myGridY)
            {
                return;
            }

            // up
            if (ix < overlapSizeUp)
            {
                size_t pmeIndex = (ix + pmeSize.x() - overlapSizeUp) * pmeSize.y() * pmeSize.z()
                                  + iy * pmeSize.z() + iz;
                size_t packedIndex    = ix * myGridY * pmeSize.z() + iy * pmeSize.z() + iz;
                gm_realGrid[pmeIndex] = gm_transferGridUp[packedIndex];
            }

            // down
            if (ix >= myGridX - overlapSizeDown)
            {
                size_t pmeIndex =
                        (ix + overlapSizeDown) * pmeSize.y() * pmeSize.z() + iy * pmeSize.z() + iz;
                size_t packedIndex = (ix - (myGridX - overlapSizeDown)) * myGridY * pmeSize.z()
                                     + iy * pmeSize.z() + iz;
                gm_realGrid[pmeIndex] = gm_transferGridDown[packedIndex];
            }

            // left
            if (iy < overlapSizeLeft)
            {
                size_t pmeIndex = ix * pmeSize.y() * pmeSize.z()
                                  + (iy + pmeSize.y() - overlapSizeLeft) * pmeSize.z() + iz;
                size_t packedIndex    = ix * overlapSizeLeft * pmeSize.z() + iy * pmeSize.z() + iz;
                gm_realGrid[pmeIndex] = gm_transferGridLeft[packedIndex];
            }

            // right
            if (iy >= myGridY - overlapSizeRight)
            {
                size_t pmeIndex =
                        ix * pmeSize.y() * pmeSize.z() + (iy + overlapSizeRight) * pmeSize.z() + iz;
                size_t packedIndex = ix * overlapSizeRight * pmeSize.z()
                                     + (iy - (myGridY - overlapSizeRight)) * pmeSize.z() + iz;
                gm_realGrid[pmeIndex] = gm_transferGridRight[packedIndex];
            }

            // up left
            if (ix < overlapSizeUp && iy < overlapSizeLeft)
            {
                size_t pmeIndex = (ix + pmeSize.x() - overlapSizeUp) * pmeSize.y() * pmeSize.z()
                                  + (iy + pmeSize.y() - overlapSizeLeft) * pmeSize.z() + iz;
                size_t packedIndex    = ix * overlapSizeLeft * pmeSize.z() + iy * pmeSize.z() + iz;
                gm_realGrid[pmeIndex] = gm_transferGridUpLeft[packedIndex];
            }

            // down left
            if (ix >= myGridX - overlapSizeDown && iy < overlapSizeLeft)
            {
                size_t pmeIndex = (ix + overlapSizeDown) * pmeSize.y() * pmeSize.z()
                                  + (iy + pmeSize.y() - overlapSizeLeft) * pmeSize.z() + iz;
                size_t packedIndex = (ix - (myGridX - overlapSizeDown)) * overlapSizeLeft * pmeSize.z()
                                     + iy * pmeSize.z() + iz;
                gm_realGrid[pmeIndex] = gm_transferGridDownLeft[packedIndex];
            }

            // up right
            if (ix < overlapSizeUp && iy >= myGridY - overlapSizeRight)
            {
                size_t pmeIndex = (ix + pmeSize.x() - overlapSizeUp) * pmeSize.y() * pmeSize.z()
                                  + (iy + overlapSizeRight) * pmeSize.z() + iz;
                size_t packedIndex = ix * overlapSizeRight * pmeSize.z()
                                     + (iy - (myGridY - overlapSizeRight)) * pmeSize.z() + iz;
                gm_realGrid[pmeIndex] = gm_transferGridUpRight[packedIndex];
            }

            // down right
            if (ix >= myGridX - overlapSizeDown && iy >= myGridY - overlapSizeRight)
            {
                size_t pmeIndex = (ix + overlapSizeDown) * pmeSize.y() * pmeSize.z()
                                  + (iy + overlapSizeRight) * pmeSize.z() + iz;
                size_t packedIndex = (ix - (myGridX - overlapSizeDown)) * overlapSizeRight * pmeSize.z()
                                     + (iy - (myGridY - overlapSizeRight)) * pmeSize.z() + iz;
                gm_realGrid[pmeIndex] = gm_transferGridDownRight[packedIndex];
            }
        };
    }
};

//! Handles a kernel which adds grid overlap data received from neighboring ranks
class UnpackAndAddHaloInternal
{
public:
    /*! \brief
     * Returns the unpack kernel
     *
     * \param[in] myGridX,myGridY            Local domain size in X and Y dimension
     * \param[in] pmeSize                    Local PME grid size
     * \param[in] gm_realGrid                PME device grid
     * \param[out] gm_transferGridUp         Device array used to pack data to go up
     * \param[out] gm_transferGridDown       Device array used to pack data to go down
     * \param[out] gm_transferGridLeft       Device array used to pack data to go left
     * \param[out] gm_transferGridRight      Device array used to pack data to go right
     * \param[out] gm_transferGridUpLeft     Device array used to pack data to go up+left
     * \param[out] gm_transferGridDownLeft   Device array used to pack data to go down+left
     * \param[out] gm_transferGridUpRight    Device array used to pack data to go up+right
     * \param[out] gm_transferGridDownRight  Device array used to pack data to go down+right
     * \param[in] overlapSizeX,overlapSizeY,overlapUp,overlapLeft
     *                                       Halo size in 4 directions
     * \tparam    subGroupSize               Size of the sub-group.
     */
    template<int subGroupSize>
    static auto kernel(size_t      myGridX,
                       size_t      myGridY,
                       sycl::uint3 pmeSize,
                       float* __restrict__ gm_realGrid,
                       const float* __restrict__ gm_transferGridUp,
                       const float* __restrict__ gm_transferGridDown,
                       const float* __restrict__ gm_transferGridLeft,
                       const float* __restrict__ gm_transferGridRight,
                       const float* __restrict__ gm_transferGridUpLeft,
                       const float* __restrict__ gm_transferGridDownLeft,
                       const float* __restrict__ gm_transferGridUpRight,
                       const float* __restrict__ gm_transferGridDownRight,
                       size_t overlapSizeX,
                       size_t overlapSizeY,
                       size_t overlapUp,
                       size_t overlapLeft)
    {
        return [=](sycl::nd_item<3> item_ct1) [[intel::reqd_sub_group_size(subGroupSize)]]
        {
            size_t iz = item_ct1.get_local_id(2) + item_ct1.get_group(2) * item_ct1.get_local_range(2);
            size_t iy = item_ct1.get_local_id(1) + item_ct1.get_group(1) * item_ct1.get_local_range(1);
            size_t ix = item_ct1.get_local_id(0) + item_ct1.get_group(0) * item_ct1.get_local_range(0);

            // we might get iz greater than pmeSize.z when pmeSize.z is not multiple of
            // threadsAlongZDim(see below), same for iy when it's not multiple of threadsAlongYDim
            if (iz >= pmeSize.z() || iy >= myGridY)
            {
                return;
            }

            size_t pmeIndex = ix * pmeSize.y() * pmeSize.z() + iy * pmeSize.z() + iz;

            float val = gm_realGrid[pmeIndex];

            // up rank
            if (ix < overlapSizeX)
            {
                size_t packedIndex = ix * myGridY * pmeSize.z() + iy * pmeSize.z() + iz;
                val += gm_transferGridUp[packedIndex];
            }

            // down rank
            if (ix >= myGridX - overlapSizeX && overlapUp > 0)
            {
                size_t packedIndex = (ix - (myGridX - overlapSizeX)) * myGridY * pmeSize.z()
                                     + iy * pmeSize.z() + iz;
                val += gm_transferGridDown[packedIndex];
            }

            // left rank
            if (iy < overlapSizeY)
            {
                size_t packedIndex = ix * overlapSizeY * pmeSize.z() + iy * pmeSize.z() + iz;
                val += gm_transferGridLeft[packedIndex];
            }

            // right rank
            if (iy >= myGridY - overlapSizeY && overlapLeft > 0)
            {
                size_t packedIndex = ix * overlapSizeY * pmeSize.z()
                                     + (iy - (myGridY - overlapSizeY)) * pmeSize.z() + iz;
                val += gm_transferGridRight[packedIndex];
            }

            // up left rank
            if (ix < overlapSizeX && iy < overlapSizeY)
            {
                size_t packedIndex = ix * overlapSizeY * pmeSize.z() + iy * pmeSize.z() + iz;
                val += gm_transferGridUpLeft[packedIndex];
            }

            // up right rank
            if (ix < overlapSizeX && iy >= myGridY - overlapSizeY && overlapLeft > 0)
            {
                size_t packedIndex = ix * overlapSizeY * pmeSize.z()
                                     + (iy - (myGridY - overlapSizeY)) * pmeSize.z() + iz;
                val += gm_transferGridUpRight[packedIndex];
            }

            // down left rank
            if (ix >= myGridX - overlapSizeX && overlapUp > 0 && iy < overlapSizeY)
            {
                size_t packedIndex = (ix - (myGridX - overlapSizeX)) * overlapSizeY * pmeSize.z()
                                     + iy * pmeSize.z() + iz;
                val += gm_transferGridDownLeft[packedIndex];
            }

            // down right rank
            if (ix >= myGridX - overlapSizeX && overlapUp > 0 && iy >= myGridY - overlapSizeY
                && overlapLeft > 0)
            {
                size_t packedIndex = (ix - (myGridX - overlapSizeX)) * overlapSizeY * pmeSize.z()
                                     + (iy - (myGridY - overlapSizeY)) * pmeSize.z() + iz;
                val += gm_transferGridDownRight[packedIndex];
            }

            gm_realGrid[pmeIndex] = val;
        };
    }
};

//! Handles a kernel which packs non-contiguous overlap data in all 8 neighboring directions
class PackHaloInternal
{
public:
    /*! \brief
     * Returns the pack kernel
     *
     * \param[in] myGridX,myGridY            Local domain size in X and Y dimension
     * \param[in] pmeSize                    Local PME grid size
     * \param[in] gm_realGrid                PME device grid
     * \param[out] gm_transferGridUp         Device array used to pack data to go up
     * \param[out] gm_transferGridDown       Device array used to pack data to go down
     * \param[out] gm_transferGridLeft       Device array used to pack data to go left
     * \param[out] gm_transferGridRight      Device array used to pack data to go right
     * \param[out] gm_transferGridUpLeft     Device array used to pack data to go up+left
     * \param[out] gm_transferGridDownLeft   Device array used to pack data to go down+left
     * \param[out] gm_transferGridUpRight    Device array used to pack data to go up+right
     * \param[out] gm_transferGridDownRight  Device array used to pack data to go down+right
     * \param[in] overlapSizeX,overlapSizeY,overlapUp,overlapLeft
     *                                       Halo size in 4 directions
     * \tparam    subGroupSize               Size of the sub-group.
     */
    template<int subGroupSize>
    static auto kernel(size_t      myGridX,
                       size_t      myGridY,
                       sycl::uint3 pmeSize,
                       const float* __restrict__ gm_realGrid,
                       float* __restrict__ gm_transferGridUp,
                       float* __restrict__ gm_transferGridDown,
                       float* __restrict__ gm_transferGridLeft,
                       float* __restrict__ gm_transferGridRight,
                       float* __restrict__ gm_transferGridUpLeft,
                       float* __restrict__ gm_transferGridDownLeft,
                       float* __restrict__ gm_transferGridUpRight,
                       float* __restrict__ gm_transferGridDownRight,
                       size_t overlapSizeX,
                       size_t overlapSizeY,
                       size_t overlapUp,
                       size_t overlapLeft)
    {
        return [=](sycl::nd_item<3> item_ct1) [[intel::reqd_sub_group_size(subGroupSize)]]
        {
            size_t iz = item_ct1.get_local_id(2) + item_ct1.get_group(2) * item_ct1.get_local_range(2);
            size_t iy = item_ct1.get_local_id(1) + item_ct1.get_group(1) * item_ct1.get_local_range(1);
            size_t ix = item_ct1.get_local_id(0) + item_ct1.get_group(0) * item_ct1.get_local_range(0);

            // we might get iz greater than pmeSize.z when pmeSize.z is not multiple of
            // threadsAlongZDim(see below), same for iy when it's not multiple of threadsAlongYDim
            if (iz >= pmeSize.z() || iy >= myGridY)
            {
                return;
            }

            size_t pmeIndex = ix * pmeSize.y() * pmeSize.z() + iy * pmeSize.z() + iz;

            float val = gm_realGrid[pmeIndex];

            // up rank
            if (ix < overlapSizeX)
            {
                size_t packedIndex             = ix * myGridY * pmeSize.z() + iy * pmeSize.z() + iz;
                gm_transferGridUp[packedIndex] = val;
            }

            // down rank
            if (ix >= myGridX - overlapSizeX && overlapUp > 0)
            {
                size_t packedIndex = (ix - (myGridX - overlapSizeX)) * myGridY * pmeSize.z()
                                     + iy * pmeSize.z() + iz;
                gm_transferGridDown[packedIndex] = val;
            }

            // left rank
            if (iy < overlapSizeY)
            {
                size_t packedIndex = ix * overlapSizeY * pmeSize.z() + iy * pmeSize.z() + iz;
                gm_transferGridLeft[packedIndex] = val;
            }

            // right rank
            if (iy >= myGridY - overlapSizeY && overlapLeft > 0)
            {
                size_t packedIndex = ix * overlapSizeY * pmeSize.z()
                                     + (iy - (myGridY - overlapSizeY)) * pmeSize.z() + iz;
                gm_transferGridRight[packedIndex] = val;
            }

            // up left rank
            if (ix < overlapSizeX && iy < overlapSizeY)
            {
                size_t packedIndex = ix * overlapSizeY * pmeSize.z() + iy * pmeSize.z() + iz;
                gm_transferGridUpLeft[packedIndex] = val;
            }

            // down left rank
            if (ix >= myGridX - overlapSizeX && overlapUp > 0 && iy < overlapSizeY)
            {
                size_t packedIndex = (ix - (myGridX - overlapSizeX)) * overlapSizeY * pmeSize.z()
                                     + iy * pmeSize.z() + iz;
                gm_transferGridDownLeft[packedIndex] = val;
            }

            // up right rank
            if (ix < overlapSizeX && iy >= myGridY - overlapSizeY && overlapLeft > 0)
            {
                size_t packedIndex = ix * overlapSizeY * pmeSize.z()
                                     + (iy - (myGridY - overlapSizeY)) * pmeSize.z() + iz;
                gm_transferGridUpRight[packedIndex] = val;
            }

            // down right rank
            if (ix >= myGridX - overlapSizeX && overlapUp > 0 && iy >= myGridY - overlapSizeY
                && overlapLeft > 0)
            {
                size_t packedIndex = (ix - (myGridX - overlapSizeX)) * overlapSizeY * pmeSize.z()
                                     + (iy - (myGridY - overlapSizeY)) * pmeSize.z() + iz;
                gm_transferGridDownRight[packedIndex] = val;
            }
        };
    }
};

/*! \brief Submits a GPU grid kernel
 *
 * \tparam    Kernel           The class containing a static kernel() method to return
 *                             the kernel to execute
 * \tparam    subGroupSize     Size of the sub-group.
 * \param[in] deviceStream     The device stream upon which to submit
 * \param[in] myGridX,myGridY  Local domain size in X and Y dimension
 * \param[in] pmeSize          Local PME grid size
 * \param[in] args             Parameter pack to pass to the kernel
 */
template<typename Kernel, int subGroupSize, class... Args>
static void
submit(const DeviceStream& deviceStream, size_t myGridX, size_t myGridY, sycl::uint3 pmeSize, Args&&... args)
{
    // Having threadsAlongZDim the same as the sub-group size is good
    // for coalescing. Current value is taken from the CUDA version.
    // It is not ideal for 64-wide hardware, but appears to work fine.
    // For Intel GPUs, we force 32-wide execution with
    // reqd_sub_group_size attribute in the kernel builder function.
    const int threadsAlongZDim = subGroupSize;
    const int threadsAlongYDim = 4;

    const sycl::range<3>    localSize{ 1, threadsAlongYDim, threadsAlongZDim };
    const sycl::range<3>    groupRange{ myGridX,
                                     gmx::divideRoundUp<size_t>(myGridY, threadsAlongYDim),
                                     gmx::divideRoundUp<size_t>(pmeSize[ZZ], threadsAlongZDim) };
    const sycl::nd_range<3> range{ groupRange * localSize, localSize };

    sycl::queue q = deviceStream.stream();
    q.submit(GMX_SYCL_DISCARD_EVENT[&](sycl::handler & cgh) {
        auto kernel = Kernel::template kernel<subGroupSize>(
                myGridX, myGridY, pmeSize, std::forward<Args>(args)...);
        cgh.parallel_for<Kernel>(range, kernel);
    });
}

#if GMX_MPI
/*! \brief
 * Utility function to send and recv halo data from neighboring ranks
 */
static void receiveAndSend(float*       sendBuf,
                           int          sendCount,
                           int          dest,
                           MPI_Request* sendRequest,
                           float*       recvBuf,
                           int          recvCount,
                           int          src,
                           MPI_Request* recvRequest,
                           int          tag,
                           MPI_Comm     comm)
{
    // send data to dest rank and recv from src rank
    MPI_Irecv(recvBuf, recvCount, MPI_FLOAT, src, tag, comm, recvRequest);

    MPI_Isend(sendBuf, sendCount, MPI_FLOAT, dest, tag, comm, sendRequest);
}
#endif

void pmeGpuGridHaloExchange(const PmeGpu* pmeGpu, gmx_wallcycle* wcycle)
{
#if GMX_MPI
    // Note here we are assuming that width of the chunks is not so small that we need to
    // transfer to/from multiple ranks i.e. that the distributed grid contains chunks at least order-1 points wide.

    auto*             kernelParamsPtr = pmeGpu->kernelParams.get();
    const sycl::uint3 localPmeSize{ kernelParamsPtr->grid.realGridSizePadded[XX],
                                    kernelParamsPtr->grid.realGridSizePadded[YY],
                                    kernelParamsPtr->grid.realGridSizePadded[ZZ] };

    size_t overlapX = pmeGpu->haloExchange->haloSizeX[gmx::DirectionX::Center];
    size_t overlapY = pmeGpu->haloExchange->haloSizeY[gmx::DirectionY::Center];

    size_t overlapDown = pmeGpu->haloExchange->haloSizeX[gmx::DirectionX::Down];
    size_t overlapUp   = pmeGpu->haloExchange->haloSizeX[gmx::DirectionX::Up];

    size_t overlapRight = pmeGpu->haloExchange->haloSizeY[gmx::DirectionY::Right];
    size_t overlapLeft  = pmeGpu->haloExchange->haloSizeY[gmx::DirectionY::Left];

    size_t myGridX = pmeGpu->haloExchange->gridSizeX;
    size_t myGridY = pmeGpu->haloExchange->gridSizeY;

    size_t sizeX = pmeGpu->common->nnodesX;
    size_t down  = pmeGpu->haloExchange->ranksX[gmx::DirectionX::Down];
    size_t up    = pmeGpu->haloExchange->ranksX[gmx::DirectionX::Up];

    size_t sizeY = pmeGpu->common->nnodesY;
    size_t right = pmeGpu->haloExchange->ranksY[gmx::DirectionY::Right];
    size_t left  = pmeGpu->haloExchange->ranksY[gmx::DirectionY::Left];

    for (int gridIndex = 0; gridIndex < pmeGpu->common->ngrids; gridIndex++)
    {
        MPI_Request          req[16];
        int                  reqCount = 0;
        DeviceBuffer<float>& realGrid = pmeGpu->kernelParams->grid.d_realGrid[gridIndex];
        float *              sendGridUp, *sendGridDown;

        // no need to pack if slab-decomposition in X-dimension as data is already contiguous
        if (pmeGpu->common->nnodesY == 1)
        {
            int sendOffsetDown = myGridX * localPmeSize[YY] * localPmeSize[ZZ];
            int sendOffsetUp = (localPmeSize[XX] - overlapUp) * localPmeSize[YY] * localPmeSize[ZZ];
            sendGridUp       = &asMpiPointer(realGrid)[sendOffsetUp];
            sendGridDown     = &asMpiPointer(realGrid)[sendOffsetDown];
        }
        else
        {
            wallcycle_start(wcycle, WallCycleCounter::LaunchGpuPme);

            // launch packing kernel
            submit<PackHaloExternal, sc_subGroupSize>(
                    pmeGpu->archSpecific->pmeStream_,
                    myGridX,
                    myGridY,
                    localPmeSize,
                    realGrid.get_pointer(),
                    pmeGpu->haloExchange->d_sendGrids[gmx::DirectionX::Up][gmx::DirectionY::Center].get_pointer(),
                    pmeGpu->haloExchange->d_sendGrids[gmx::DirectionX::Down][gmx::DirectionY::Center]
                            .get_pointer(),
                    pmeGpu->haloExchange->d_sendGrids[gmx::DirectionX::Center][gmx::DirectionY::Left]
                            .get_pointer(),
                    pmeGpu->haloExchange
                            ->d_sendGrids[gmx::DirectionX::Center][gmx::DirectionY::Right]
                            .get_pointer(),
                    pmeGpu->haloExchange->d_sendGrids[gmx::DirectionX::Up][gmx::DirectionY::Left].get_pointer(),
                    pmeGpu->haloExchange->d_sendGrids[gmx::DirectionX::Down][gmx::DirectionY::Left].get_pointer(),
                    pmeGpu->haloExchange->d_sendGrids[gmx::DirectionX::Up][gmx::DirectionY::Right].get_pointer(),
                    pmeGpu->haloExchange->d_sendGrids[gmx::DirectionX::Down][gmx::DirectionY::Right].get_pointer(),
                    overlapUp,
                    overlapDown,
                    overlapLeft,
                    overlapRight);
            sendGridUp = asMpiPointer(
                    pmeGpu->haloExchange->d_sendGrids[gmx::DirectionX::Up][gmx::DirectionY::Center]);
            sendGridDown = asMpiPointer(
                    pmeGpu->haloExchange->d_sendGrids[gmx::DirectionX::Down][gmx::DirectionY::Center]);

            wallcycle_stop(wcycle, WallCycleCounter::LaunchGpuPme);
        }

        wallcycle_start(wcycle, WallCycleCounter::WaitGpuPmeSpread);

        // Make sure data is ready on GPU before MPI communication.
        // Wait for spread to finish in case of slab decomposition along X-dimension and
        // wait for packing to finish otherwise.
        // Todo: Consider using events to create dependcy on spread
        pmeGpu->archSpecific->pmeStream_.synchronize();

        wallcycle_stop(wcycle, WallCycleCounter::WaitGpuPmeSpread);


        wallcycle_start(wcycle, WallCycleCounter::PmeHaloExchangeComm);

        // major dimension
        if (sizeX > 1)
        {
            constexpr int mpiTag = 403; // Arbitrarily chosen

            // send data to down rank and recv from up rank
            receiveAndSend(
                    sendGridDown,
                    overlapDown * myGridY * localPmeSize[ZZ],
                    down,
                    &req[reqCount],
                    asMpiPointer(pmeGpu->haloExchange->d_recvGrids[gmx::DirectionX::Up][gmx::DirectionY::Center]),
                    overlapX * myGridY * localPmeSize[ZZ],
                    up,
                    &req[reqCount + 1],
                    mpiTag,
                    pmeGpu->common->mpiCommX);
            reqCount += 2;

            if (overlapUp > 0)
            {
                // send data to up rank and recv from down rank
                receiveAndSend(
                        sendGridUp,
                        overlapUp * myGridY * localPmeSize[ZZ],
                        up,
                        &req[reqCount],
                        asMpiPointer(pmeGpu->haloExchange->d_recvGrids[gmx::DirectionX::Down][gmx::DirectionY::Center]),
                        overlapX * myGridY * localPmeSize[ZZ],
                        down,
                        &req[reqCount + 1],
                        mpiTag,
                        pmeGpu->common->mpiCommX);
                reqCount += 2;
            }
        }

        // minor dimension
        if (sizeY > 1)
        {
            constexpr int mpiTag = 404; // Arbitrarily chosen

            // recv from left rank and send to right rank
            receiveAndSend(
                    asMpiPointer(pmeGpu->haloExchange->d_sendGrids[gmx::DirectionX::Center][gmx::DirectionY::Right]),
                    overlapRight * myGridX * localPmeSize[ZZ],
                    right,
                    &req[reqCount],
                    asMpiPointer(pmeGpu->haloExchange->d_recvGrids[gmx::DirectionX::Center][gmx::DirectionY::Left]),
                    overlapY * myGridX * localPmeSize[ZZ],
                    left,
                    &req[reqCount + 1],
                    mpiTag,
                    pmeGpu->common->mpiCommY);
            reqCount += 2;

            if (overlapLeft > 0)
            {
                // recv from right rank and send data to left rank
                receiveAndSend(
                        asMpiPointer(pmeGpu->haloExchange->d_sendGrids[gmx::DirectionX::Center][gmx::DirectionY::Left]),
                        overlapLeft * myGridX * localPmeSize[ZZ],
                        left,
                        &req[reqCount],
                        asMpiPointer(pmeGpu->haloExchange->d_recvGrids[gmx::DirectionX::Center][gmx::DirectionY::Right]),
                        overlapY * myGridX * localPmeSize[ZZ],
                        right,
                        &req[reqCount + 1],
                        mpiTag,
                        pmeGpu->common->mpiCommY);
                reqCount += 2;
            }
        }

        if (sizeX > 1 && sizeY > 1)
        {
            int rankUpLeft   = up * sizeY + left;
            int rankDownLeft = down * sizeY + left;

            int rankUpRight   = up * sizeY + right;
            int rankDownRight = down * sizeY + right;

            constexpr int mpiTag = 405; // Arbitrarily chosen

            // send data to down rank and recv from up rank
            receiveAndSend(
                    asMpiPointer(pmeGpu->haloExchange->d_sendGrids[gmx::DirectionX::Down][gmx::DirectionY::Right]),
                    overlapDown * overlapRight * localPmeSize[ZZ],
                    rankDownRight,
                    &req[reqCount],
                    asMpiPointer(pmeGpu->haloExchange->d_recvGrids[gmx::DirectionX::Up][gmx::DirectionY::Left]),
                    overlapX * overlapY * localPmeSize[ZZ],
                    rankUpLeft,
                    &req[reqCount + 1],
                    mpiTag,
                    pmeGpu->common->mpiComm);
            reqCount += 2;

            if (overlapLeft > 0)
            {
                // send data to down left rank and recv from up right rank
                receiveAndSend(
                        asMpiPointer(pmeGpu->haloExchange->d_sendGrids[gmx::DirectionX::Down][gmx::DirectionY::Left]),
                        overlapDown * overlapLeft * localPmeSize[ZZ],
                        rankDownLeft,
                        &req[reqCount],
                        asMpiPointer(pmeGpu->haloExchange->d_recvGrids[gmx::DirectionX::Up][gmx::DirectionY::Right]),
                        overlapX * overlapY * localPmeSize[ZZ],
                        rankUpRight,
                        &req[reqCount + 1],
                        mpiTag,
                        pmeGpu->common->mpiComm);
                reqCount += 2;
            }

            if (overlapUp > 0)
            {
                // send data to up right rank and recv from down left rank
                receiveAndSend(
                        asMpiPointer(pmeGpu->haloExchange->d_sendGrids[gmx::DirectionX::Up][gmx::DirectionY::Right]),
                        overlapUp * overlapRight * localPmeSize[ZZ],
                        rankUpRight,
                        &req[reqCount],
                        asMpiPointer(pmeGpu->haloExchange->d_recvGrids[gmx::DirectionX::Down][gmx::DirectionY::Left]),
                        overlapX * overlapY * localPmeSize[ZZ],
                        rankDownLeft,
                        &req[reqCount + 1],
                        mpiTag,
                        pmeGpu->common->mpiComm);
                reqCount += 2;
            }

            if (overlapUp > 0 && overlapLeft > 0)
            {
                // send data to up left rank and recv from down right rank
                receiveAndSend(
                        asMpiPointer(pmeGpu->haloExchange->d_sendGrids[gmx::DirectionX::Up][gmx::DirectionY::Left]),
                        overlapUp * overlapLeft * localPmeSize[ZZ],
                        rankUpLeft,
                        &req[reqCount],
                        asMpiPointer(pmeGpu->haloExchange->d_recvGrids[gmx::DirectionX::Down][gmx::DirectionY::Right]),
                        overlapX * overlapY * localPmeSize[ZZ],
                        rankDownRight,
                        &req[reqCount + 1],
                        mpiTag,
                        pmeGpu->common->mpiComm);
                reqCount += 2;
            }
        }

        MPI_Waitall(reqCount, req, MPI_STATUSES_IGNORE);

        wallcycle_stop(wcycle, WallCycleCounter::PmeHaloExchangeComm);

        wallcycle_start(wcycle, WallCycleCounter::LaunchGpuPme);

        // reduce halo data
        submit<UnpackAndAddHaloInternal, sc_subGroupSize>(
                pmeGpu->archSpecific->pmeStream_,
                myGridX,
                myGridY,
                localPmeSize,
                realGrid.get_pointer(),
                pmeGpu->haloExchange->d_recvGrids[gmx::DirectionX::Up][gmx::DirectionY::Center].get_pointer(),
                pmeGpu->haloExchange->d_recvGrids[gmx::DirectionX::Down][gmx::DirectionY::Center].get_pointer(),
                pmeGpu->haloExchange->d_recvGrids[gmx::DirectionX::Center][gmx::DirectionY::Left].get_pointer(),
                pmeGpu->haloExchange->d_recvGrids[gmx::DirectionX::Center][gmx::DirectionY::Right].get_pointer(),
                pmeGpu->haloExchange->d_recvGrids[gmx::DirectionX::Up][gmx::DirectionY::Left].get_pointer(),
                pmeGpu->haloExchange->d_recvGrids[gmx::DirectionX::Down][gmx::DirectionY::Left].get_pointer(),
                pmeGpu->haloExchange->d_recvGrids[gmx::DirectionX::Up][gmx::DirectionY::Right].get_pointer(),
                pmeGpu->haloExchange->d_recvGrids[gmx::DirectionX::Down][gmx::DirectionY::Right].get_pointer(),
                overlapX,
                overlapY,
                overlapUp,
                overlapLeft);

        wallcycle_stop(wcycle, WallCycleCounter::LaunchGpuPme);
    }
#else
    GMX_UNUSED_VALUE(pmeGpu);
#endif
}

void pmeGpuGridHaloExchangeReverse(const PmeGpu* pmeGpu, gmx_wallcycle* wcycle)
{
#if GMX_MPI
    // Note here we are assuming that width of the chunks is not so small that we need to
    // transfer to/from multiple ranks i.e. that the distributed grid contains chunks at least order-1 points wide.

    auto*             kernelParamsPtr = pmeGpu->kernelParams.get();
    const sycl::uint3 localPmeSize{ kernelParamsPtr->grid.realGridSizePadded[XX],
                                    kernelParamsPtr->grid.realGridSizePadded[YY],
                                    kernelParamsPtr->grid.realGridSizePadded[ZZ] };

    size_t overlapX = pmeGpu->haloExchange->haloSizeX[gmx::DirectionX::Center];
    size_t overlapY = pmeGpu->haloExchange->haloSizeY[gmx::DirectionY::Center];

    size_t overlapDown = pmeGpu->haloExchange->haloSizeX[gmx::DirectionX::Down];
    size_t overlapUp   = pmeGpu->haloExchange->haloSizeX[gmx::DirectionX::Up];

    size_t overlapRight = pmeGpu->haloExchange->haloSizeY[gmx::DirectionY::Right];
    size_t overlapLeft  = pmeGpu->haloExchange->haloSizeY[gmx::DirectionY::Left];

    size_t myGridX = pmeGpu->haloExchange->gridSizeX;
    size_t myGridY = pmeGpu->haloExchange->gridSizeY;

    size_t sizeX = pmeGpu->common->nnodesX;
    size_t down  = pmeGpu->haloExchange->ranksX[gmx::DirectionX::Down];
    size_t up    = pmeGpu->haloExchange->ranksX[gmx::DirectionX::Up];

    size_t sizeY = pmeGpu->common->nnodesY;
    size_t right = pmeGpu->haloExchange->ranksY[gmx::DirectionY::Right];
    size_t left  = pmeGpu->haloExchange->ranksY[gmx::DirectionY::Left];

    for (int gridIndex = 0; gridIndex < pmeGpu->common->ngrids; gridIndex++)
    {
        MPI_Request req[16];
        int         reqCount = 0;

        DeviceBuffer<float>& realGrid = pmeGpu->kernelParams->grid.d_realGrid[gridIndex];
        float *              sendGridUp, *sendGridDown, *recvGridUp, *recvGridDown;

        // no need to pack if slab-decomposition in X-dimension as data is already contiguous
        if (sizeY == 1)
        {
            int sendOffsetUp   = 0;
            int sendOffsetDown = (myGridX - overlapX) * localPmeSize[YY] * localPmeSize[ZZ];
            int recvOffsetUp = (localPmeSize[XX] - overlapUp) * localPmeSize[YY] * localPmeSize[ZZ];
            int recvOffsetDown = myGridX * localPmeSize[YY] * localPmeSize[ZZ];
            sendGridUp         = &asMpiPointer(realGrid)[sendOffsetUp];
            sendGridDown       = &asMpiPointer(realGrid)[sendOffsetDown];
            recvGridUp         = &asMpiPointer(realGrid)[recvOffsetUp];
            recvGridDown       = &asMpiPointer(realGrid)[recvOffsetDown];
        }
        else
        {
            wallcycle_start(wcycle, WallCycleCounter::LaunchGpuPme);

            // launch packing kernel
            submit<PackHaloInternal, sc_subGroupSize>(
                    pmeGpu->archSpecific->pmeStream_,
                    myGridX,
                    myGridY,
                    localPmeSize,
                    realGrid.get_pointer(),
                    pmeGpu->haloExchange->d_sendGrids[gmx::DirectionX::Up][gmx::DirectionY::Center].get_pointer(),
                    pmeGpu->haloExchange->d_sendGrids[gmx::DirectionX::Down][gmx::DirectionY::Center]
                            .get_pointer(),
                    pmeGpu->haloExchange->d_sendGrids[gmx::DirectionX::Center][gmx::DirectionY::Left]
                            .get_pointer(),
                    pmeGpu->haloExchange
                            ->d_sendGrids[gmx::DirectionX::Center][gmx::DirectionY::Right]
                            .get_pointer(),
                    pmeGpu->haloExchange->d_sendGrids[gmx::DirectionX::Up][gmx::DirectionY::Left].get_pointer(),
                    pmeGpu->haloExchange->d_sendGrids[gmx::DirectionX::Down][gmx::DirectionY::Left].get_pointer(),
                    pmeGpu->haloExchange->d_sendGrids[gmx::DirectionX::Up][gmx::DirectionY::Right].get_pointer(),
                    pmeGpu->haloExchange->d_sendGrids[gmx::DirectionX::Down][gmx::DirectionY::Right].get_pointer(),
                    overlapX,
                    overlapY,
                    overlapUp,
                    overlapLeft);

            sendGridUp = asMpiPointer(
                    pmeGpu->haloExchange->d_sendGrids[gmx::DirectionX::Up][gmx::DirectionY::Center]);
            sendGridDown = asMpiPointer(
                    pmeGpu->haloExchange->d_sendGrids[gmx::DirectionX::Down][gmx::DirectionY::Center]);
            recvGridUp = asMpiPointer(
                    pmeGpu->haloExchange->d_recvGrids[gmx::DirectionX::Up][gmx::DirectionY::Center]);
            recvGridDown = asMpiPointer(
                    pmeGpu->haloExchange->d_recvGrids[gmx::DirectionX::Down][gmx::DirectionY::Center]);
            wallcycle_stop(wcycle, WallCycleCounter::LaunchGpuPme);
        }

        wallcycle_start(wcycle, WallCycleCounter::WaitGpuFftToPmeGrid);

        // Make sure data is ready on GPU before MPI communication.
        // Wait for FFT to PME grid conversion to finish in case of slab decomposition along X-dimension and
        // wait for packing to finish otherwise.
        // Todo: Consider using events to create dependcy on FFT->PME grid operation
        pmeGpu->archSpecific->pmeStream_.synchronize();

        wallcycle_stop(wcycle, WallCycleCounter::WaitGpuFftToPmeGrid);

        wallcycle_start(wcycle, WallCycleCounter::PmeHaloExchangeComm);

        // major dimension
        if (sizeX > 1)
        {
            constexpr int mpiTag = 406; // Arbitrarily chosen

            // send data to up rank and recv from down rank
            receiveAndSend(sendGridUp,
                           overlapX * myGridY * localPmeSize[ZZ],
                           up,
                           &req[reqCount],
                           recvGridDown,
                           overlapDown * myGridY * localPmeSize[ZZ],
                           down,
                           &req[reqCount + 1],
                           mpiTag,
                           pmeGpu->common->mpiCommX);
            reqCount += 2;

            if (overlapUp > 0)
            {
                // send data to down rank and recv from up rank
                receiveAndSend(sendGridDown,
                               overlapX * myGridY * localPmeSize[ZZ],
                               down,
                               &req[reqCount],
                               recvGridUp,
                               overlapUp * myGridY * localPmeSize[ZZ],
                               up,
                               &req[reqCount + 1],
                               mpiTag,
                               pmeGpu->common->mpiCommX);
                reqCount += 2;
            }
        }

        // minor dimension
        if (sizeY > 1)
        {
            constexpr int mpiTag = 407; // Arbitrarily chosen

            // recv from right rank and send data to left rank
            receiveAndSend(
                    asMpiPointer(pmeGpu->haloExchange->d_sendGrids[gmx::DirectionX::Center][gmx::DirectionY::Left]),
                    overlapY * myGridX * localPmeSize[ZZ],
                    left,
                    &req[reqCount],
                    asMpiPointer(pmeGpu->haloExchange->d_recvGrids[gmx::DirectionX::Center][gmx::DirectionY::Right]),
                    overlapRight * myGridX * localPmeSize[ZZ],
                    right,
                    &req[reqCount + 1],
                    mpiTag,
                    pmeGpu->common->mpiCommY);
            reqCount += 2;

            if (overlapLeft > 0)
            {
                // recv from left rank and send data to right rank
                receiveAndSend(
                        asMpiPointer(pmeGpu->haloExchange->d_sendGrids[gmx::DirectionX::Center][gmx::DirectionY::Right]),
                        overlapY * myGridX * localPmeSize[ZZ],
                        right,
                        &req[reqCount],
                        asMpiPointer(pmeGpu->haloExchange->d_recvGrids[gmx::DirectionX::Center][gmx::DirectionY::Left]),
                        overlapLeft * myGridX * localPmeSize[ZZ],
                        left,
                        &req[reqCount + 1],
                        mpiTag,
                        pmeGpu->common->mpiCommY);
                reqCount += 2;
            }
        }

        if (sizeX > 1 && sizeY > 1)
        {
            int rankUpLeft   = up * sizeY + left;
            int rankDownLeft = down * sizeY + left;

            int rankUpRight   = up * sizeY + right;
            int rankDownRight = down * sizeY + right;

            constexpr int mpiTag = 408; // Arbitrarily chosen

            // send data to up left and recv from down right rank
            receiveAndSend(
                    asMpiPointer(pmeGpu->haloExchange->d_sendGrids[gmx::DirectionX::Up][gmx::DirectionY::Left]),
                    overlapX * overlapY * localPmeSize[ZZ],
                    rankUpLeft,
                    &req[reqCount],
                    asMpiPointer(pmeGpu->haloExchange->d_recvGrids[gmx::DirectionX::Down][gmx::DirectionY::Right]),
                    overlapDown * overlapRight * localPmeSize[ZZ],
                    rankDownRight,
                    &req[reqCount + 1],
                    mpiTag,
                    pmeGpu->common->mpiComm);
            reqCount += 2;

            if (overlapLeft > 0)
            {
                // send data to up right rank and recv from down left rank
                receiveAndSend(
                        asMpiPointer(pmeGpu->haloExchange->d_sendGrids[gmx::DirectionX::Up][gmx::DirectionY::Right]),
                        overlapX * overlapY * localPmeSize[ZZ],
                        rankUpRight,
                        &req[reqCount],
                        asMpiPointer(pmeGpu->haloExchange->d_recvGrids[gmx::DirectionX::Down][gmx::DirectionY::Left]),
                        overlapDown * overlapLeft * localPmeSize[ZZ],
                        rankDownLeft,
                        &req[reqCount + 1],
                        mpiTag,
                        pmeGpu->common->mpiComm);
                reqCount += 2;
            }

            if (overlapUp > 0)
            {
                // send data to down left rank and recv from up right rank
                receiveAndSend(
                        asMpiPointer(pmeGpu->haloExchange->d_sendGrids[gmx::DirectionX::Down][gmx::DirectionY::Left]),
                        overlapX * overlapY * localPmeSize[ZZ],
                        rankDownLeft,
                        &req[reqCount],
                        asMpiPointer(pmeGpu->haloExchange->d_recvGrids[gmx::DirectionX::Up][gmx::DirectionY::Right]),
                        overlapUp * overlapRight * localPmeSize[ZZ],
                        rankUpRight,
                        &req[reqCount + 1],
                        mpiTag,
                        pmeGpu->common->mpiComm);
                reqCount += 2;
            }

            if (overlapUp > 0 && overlapLeft > 0)
            {
                // send data to down right rank and recv from up left rank
                receiveAndSend(
                        asMpiPointer(pmeGpu->haloExchange->d_sendGrids[gmx::DirectionX::Down][gmx::DirectionY::Right]),
                        overlapX * overlapY * localPmeSize[ZZ],
                        rankDownRight,
                        &req[reqCount],
                        asMpiPointer(pmeGpu->haloExchange->d_recvGrids[gmx::DirectionX::Up][gmx::DirectionY::Left]),
                        overlapUp * overlapLeft * localPmeSize[ZZ],
                        rankUpLeft,
                        &req[reqCount + 1],
                        mpiTag,
                        pmeGpu->common->mpiComm);
                reqCount += 2;
            }
        }

        MPI_Waitall(reqCount, req, MPI_STATUSES_IGNORE);

        wallcycle_stop(wcycle, WallCycleCounter::PmeHaloExchangeComm);

        // data is written at the right place as part of MPI communication if slab-decomposition is
        // used in X-dimension, but we need to unpack if decomposition happens (also) along Y
        if (sizeY > 1)
        {
            wallcycle_start(wcycle, WallCycleCounter::LaunchGpuPme);

            // assign halo data
            submit<UnpackHaloExternal, sc_subGroupSize>(
                    pmeGpu->archSpecific->pmeStream_,
                    myGridX,
                    myGridY,
                    localPmeSize,
                    realGrid.get_pointer(),
                    pmeGpu->haloExchange->d_recvGrids[gmx::DirectionX::Up][gmx::DirectionY::Center].get_pointer(),
                    pmeGpu->haloExchange->d_recvGrids[gmx::DirectionX::Down][gmx::DirectionY::Center]
                            .get_pointer(),
                    pmeGpu->haloExchange->d_recvGrids[gmx::DirectionX::Center][gmx::DirectionY::Left]
                            .get_pointer(),
                    pmeGpu->haloExchange
                            ->d_recvGrids[gmx::DirectionX::Center][gmx::DirectionY::Right]
                            .get_pointer(),
                    pmeGpu->haloExchange->d_recvGrids[gmx::DirectionX::Up][gmx::DirectionY::Left].get_pointer(),
                    pmeGpu->haloExchange->d_recvGrids[gmx::DirectionX::Down][gmx::DirectionY::Left].get_pointer(),
                    pmeGpu->haloExchange->d_recvGrids[gmx::DirectionX::Up][gmx::DirectionY::Right].get_pointer(),
                    pmeGpu->haloExchange->d_recvGrids[gmx::DirectionX::Down][gmx::DirectionY::Right].get_pointer(),
                    overlapUp,
                    overlapDown,
                    overlapLeft,
                    overlapRight);

            wallcycle_stop(wcycle, WallCycleCounter::LaunchGpuPme);
        }
    }
#else
    GMX_UNUSED_VALUE(pmeGpu);
#endif
}

/*! \brief Builds a kernel to convert between PME and FFT grids
 *
 * \tparam  pmeToFft  A boolean which tells if this is conversion from PME grid to FFT grid or reverse
 */
template<bool pmeToFft>
class GridConverter
{
public:
    /*! \brief
     * Returns the conversion kernel
     *
     * \param[in] fftNData           Local FFT grid size without padding
     * \param[inout] d_realGrid      Local PME grid
     * \param[inout] d_fftGrid       Local FFT grid
     * \param[in] fftSize            Local FFT grid padded size
     * \param[in] pmeSize            Local PME grid padded size
     * \tparam    subGroupSize       Size of the sub-group.
     */
    template<int subGroupSize>
    static auto convertKernel(sycl::uint3 fftNData,
                              float* __restrict__ d_realGrid,
                              float* __restrict__ d_fftGrid,
                              sycl::uint3 fftSize,
                              sycl::uint3 pmeSize)
    {
        return [=](sycl::nd_item<3> item_ct1) [[intel::reqd_sub_group_size(subGroupSize)]]
        {
            size_t iz = item_ct1.get_local_id(2) + item_ct1.get_group(2) * item_ct1.get_local_range(2);
            size_t iy = item_ct1.get_local_id(1) + item_ct1.get_group(1) * item_ct1.get_local_range(1);
            size_t ix = item_ct1.get_local_id(0) + item_ct1.get_group(0) * item_ct1.get_local_range(0);

            if (ix >= fftNData.x() || iy >= fftNData.y() || iz >= fftNData.z())
            {
                return;
            }

            size_t fftidx   = ix * fftSize.y() * fftSize.z() + iy * fftSize.z() + iz;
            size_t pmeIndex = ix * pmeSize.y() * pmeSize.z() + iy * pmeSize.z() + iz;

            if constexpr (pmeToFft)
            {
                d_fftGrid[fftidx] = d_realGrid[pmeIndex];
            }
            else
            {
                d_realGrid[pmeIndex] = d_fftGrid[fftidx];
            }
        };
    }

    /*! \brief Submits kernel
     *
     * \param[in] deviceStream   The device stream upon which to submit
     * \param[in] localFftNData  Local FFT grid size without padding
     * \param[in] args           Parameter pack to pass to the kernel builder
     * \tparam    subGroupSize   Size of the sub-group.
     */
    template<int subGroupSize, class... Args>
    static void submit(const DeviceStream& deviceStream, sycl::uint3 localFftNData, Args&&... args)
    {
        // Having threadsAlongZDim the same as the sub-group size is good
        // for coalescing. Current value is taken from the CUDA version.
        // It is not ideal for 64-wide hardware, but appears to work fine.
        // For Intel GPUs, we force 32-wide execution with
        // reqd_sub_group_size attribute in the kernel builder function.
        const int threadsAlongZDim = subGroupSize;
        const int threadsAlongYDim = 4;

        const sycl::range<3>    localSize{ 1, threadsAlongYDim, threadsAlongZDim };
        const sycl::range<3>    groupRange{ localFftNData[XX],
                                         gmx::divideRoundUp<size_t>(localFftNData[YY], threadsAlongYDim),
                                         gmx::divideRoundUp<size_t>(localFftNData[ZZ], threadsAlongZDim) };
        const sycl::nd_range<3> range{ groupRange * localSize, localSize };

        sycl::queue q = deviceStream.stream();

        q.submit(GMX_SYCL_DISCARD_EVENT[&](sycl::handler & cgh) {
            auto kernel = convertKernel<subGroupSize>(localFftNData, std::forward<Args>(args)...);
            cgh.parallel_for<GridConverter<pmeToFft>>(range, kernel);
        });
    }
};

template<bool pmeToFft>
void convertPmeGridToFftGrid(const PmeGpu* pmeGpu, float* h_fftRealGrid, gmx_parallel_3dfft* fftSetup, const int gridIndex)
{
    ivec localFftNDataAsIvec, localFftOffset, localFftSizeAsIvec;

    gmx_parallel_3dfft_real_limits(fftSetup, localFftNDataAsIvec, localFftOffset, localFftSizeAsIvec);
    const sycl::uint3 localFftNData = { localFftNDataAsIvec[XX],
                                        localFftNDataAsIvec[YY],
                                        localFftNDataAsIvec[ZZ] };
    const sycl::uint3 localFftSize  = { localFftSizeAsIvec[XX],
                                       localFftSizeAsIvec[YY],
                                       localFftSizeAsIvec[ZZ] };
    const sycl::uint3 localPmeSize  = { pmeGpu->kernelParams->grid.realGridSizePadded[XX],
                                       pmeGpu->kernelParams->grid.realGridSizePadded[YY],
                                       pmeGpu->kernelParams->grid.realGridSizePadded[ZZ] };

    // this is true in case of slab decomposition
    if (localPmeSize[ZZ] == localFftSize[ZZ] && localPmeSize[YY] == localFftSize[YY])
    {
        const int fftSize = localFftSize[ZZ] * localFftSize[YY] * localFftNData[XX];
        if (pmeToFft)
        {
            copyFromDeviceBuffer(h_fftRealGrid,
                                 &pmeGpu->kernelParams->grid.d_realGrid[gridIndex],
                                 0,
                                 fftSize,
                                 pmeGpu->archSpecific->pmeStream_,
                                 pmeGpu->settings.transferKind,
                                 nullptr);
        }
        else
        {
            copyToDeviceBuffer(&pmeGpu->kernelParams->grid.d_realGrid[gridIndex],
                               h_fftRealGrid,
                               0,
                               fftSize,
                               pmeGpu->archSpecific->pmeStream_,
                               pmeGpu->settings.transferKind,
                               nullptr);
        }
    }
    else
    {
        // launch copy kernel
        GridConverter<pmeToFft>::template submit<sc_subGroupSize>(
                pmeGpu->archSpecific->pmeStream_,
                localFftNData,
                pmeGpu->kernelParams->grid.d_realGrid[gridIndex].get_pointer(),
                h_fftRealGrid,
                localFftSize,
                localPmeSize);
    }

    if (pmeToFft)
    {
        pmeGpu->archSpecific->syncSpreadGridD2H.markEvent(pmeGpu->archSpecific->pmeStream_);
    }
}

template<bool pmeToFft>
void convertPmeGridToFftGrid(const PmeGpu* pmeGpu, DeviceBuffer<float>* d_fftRealGrid, const int gridIndex)
{
    const sycl::uint3 localPmeSize{ pmeGpu->kernelParams->grid.realGridSizePadded[XX],
                                    pmeGpu->kernelParams->grid.realGridSizePadded[YY],
                                    pmeGpu->kernelParams->grid.realGridSizePadded[ZZ] };
    const sycl::uint3 localFftNData{ pmeGpu->archSpecific->localRealGridSize[XX],
                                     pmeGpu->archSpecific->localRealGridSize[YY],
                                     pmeGpu->archSpecific->localRealGridSize[ZZ] };
    const sycl::uint3 localFftSize{ pmeGpu->archSpecific->localRealGridSizePadded[XX],
                                    pmeGpu->archSpecific->localRealGridSizePadded[YY],
                                    pmeGpu->archSpecific->localRealGridSizePadded[ZZ] };

    // this is true in case of slab decomposition
    if (localPmeSize[ZZ] == localFftSize[ZZ] && localPmeSize[YY] == localFftSize[YY])
    {
        int fftSize = localFftSize[ZZ] * localFftSize[YY] * localFftNData[XX];
        if (pmeToFft)
        {
            copyBetweenDeviceBuffers(d_fftRealGrid,
                                     &pmeGpu->kernelParams->grid.d_realGrid[gridIndex],
                                     fftSize,
                                     pmeGpu->archSpecific->pmeStream_,
                                     pmeGpu->settings.transferKind,
                                     nullptr);
        }
        else
        {
            copyBetweenDeviceBuffers(&pmeGpu->kernelParams->grid.d_realGrid[gridIndex],
                                     d_fftRealGrid,
                                     fftSize,
                                     pmeGpu->archSpecific->pmeStream_,
                                     pmeGpu->settings.transferKind,
                                     nullptr);
        }
    }
    else
    {
        // launch copy kernel
        GridConverter<pmeToFft>::template submit<sc_subGroupSize>(
                pmeGpu->archSpecific->pmeStream_,
                localFftNData,
                pmeGpu->kernelParams->grid.d_realGrid[gridIndex].get_pointer(),
                d_fftRealGrid->get_pointer(),
                localFftSize,
                localPmeSize);
    }
}

template void convertPmeGridToFftGrid<true>(const PmeGpu*       pmeGpu,
                                            float*              h_fftRealGrid,
                                            gmx_parallel_3dfft* fftSetup,
                                            const int           gridIndex);

template void convertPmeGridToFftGrid<false>(const PmeGpu*       pmeGpu,
                                             float*              h_fftRealGrid,
                                             gmx_parallel_3dfft* fftSetup,
                                             const int           gridIndex);

template void convertPmeGridToFftGrid<true>(const PmeGpu*        pmeGpu,
                                            DeviceBuffer<float>* d_fftRealGrid,
                                            const int            gridIndex);

template void convertPmeGridToFftGrid<false>(const PmeGpu*        pmeGpu,
                                             DeviceBuffer<float>* d_fftRealGrid,
                                             const int            gridIndex);
