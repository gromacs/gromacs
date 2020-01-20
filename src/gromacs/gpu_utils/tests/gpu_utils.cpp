/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2020, by the GROMACS development team, led by
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
 * \brief
 * Tests for CUDA float3 type layout.
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 */
#include "gmxpre.h"

#include "config.h"

#include <vector>

#ifndef __CUDA_ARCH__
/*! \brief Dummy definition to avoid compiler error
 *
 * \todo Find a better solution. Probably, move asFloat3(...) function to different header.
 */
#    define __CUDA_ARCH__ -1
#    include <cuda_runtime.h>
#    undef __CUDA_ARCH__
#else
#    include <cuda_runtime.h>
#endif
#include <gtest/gtest.h>

#include "gromacs/gpu_utils/cudautils.cuh"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/real.h"

#if GMX_GPU == GMX_GPU_CUDA

namespace gmx
{

namespace test
{

TEST(GpuDataTypesCompatibilityTest, RVecAndFloat3)
{
    std::vector<RVec> dataRVec;
    dataRVec.emplace_back(1.0, 2.0, 3.0);
    dataRVec.emplace_back(4.0, 5.0, 6.0);
    float3* dataFloat3 = asFloat3(dataRVec.data());
    EXPECT_EQ(dataFloat3[0].x, dataRVec[0][XX]);
    EXPECT_EQ(dataFloat3[0].y, dataRVec[0][YY]);
    EXPECT_EQ(dataFloat3[0].z, dataRVec[0][ZZ]);
    EXPECT_EQ(dataFloat3[1].x, dataRVec[1][XX]);
    EXPECT_EQ(dataFloat3[1].y, dataRVec[1][YY]);
    EXPECT_EQ(dataFloat3[1].z, dataRVec[1][ZZ]);
}

} // namespace test
} // namespace gmx

#endif // GMX_GPU == GMX_GPU_CUDA