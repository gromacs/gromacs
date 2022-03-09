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
 * \brief
 * Tests utilities for fft calculations.
 *
 * \author Gaurav Garg <gaugarg@nvidia.com>
 * \ingroup module_fft
 */
#include "gmxpre.h"

#include "config.h"

#include <algorithm>
#include <random>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/fft/fft.h"
#include "gromacs/fft/gpu_3dfft.h"
#include "gromacs/gpu_utils/clfftinitializer.h"
#if GMX_GPU
#    include "gromacs/gpu_utils/devicebuffer.h"
#endif
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/mpitest.h"
#include "testutils/refdata.h"
#include "testutils/test_hardware_environment.h"
#include "testutils/testasserts.h"
#include "testutils/testmatchers.h"

namespace gmx
{
namespace test
{
using GpuFftTestParams = std::tuple<IVec, // size of grid
                                    int,  // domains in x
                                    int,  // domains in y
                                    FftBackend>;

/*! \brief Check that the real grid after forward and backward
 * 3D transforms matches the input real grid. */
static void checkRealGrid(const IVec           realGridSizeFull,
                          const ivec           realGridSize,
                          const ivec           realGridSizePadded,
                          ArrayRef<const real> inputRealGrid,
                          ArrayRef<real>       outputRealGridValues)
{
    // Normalize the output (as the implementation does not
    // normalize either FFT)
    const real normalizationConstant =
            1.0 / (realGridSizeFull[XX] * realGridSizeFull[YY] * realGridSizeFull[ZZ]);
    std::transform(outputRealGridValues.begin(),
                   outputRealGridValues.end(),
                   outputRealGridValues.begin(),
                   [normalizationConstant](const real r) { return r * normalizationConstant; });
    // Check the real grid, skipping unused data from the padding
    const auto realGridTolerance = relativeToleranceAsFloatingPoint(10, 1e-6);
    for (int i = 0; i < realGridSize[XX] * realGridSize[YY]; i++)
    {
        auto expected =
                arrayRefFromArray(inputRealGrid.data() + i * realGridSizePadded[ZZ], realGridSize[ZZ]);
        auto actual = arrayRefFromArray(outputRealGridValues.data() + i * realGridSizePadded[ZZ],
                                        realGridSize[ZZ]);
        EXPECT_THAT(actual, Pointwise(RealEq(realGridTolerance), expected))
                << formatString("checking backward transform part %d", i);
    }
}

class GpuFftTest3D : public ::testing::Test, public ::testing::WithParamInterface<GpuFftTestParams>
{
public:
    GpuFftTest3D() = default;


    //! The whole logic being tested is contained here
    static void runTest(const GpuFftTestParams& param)
    {
        const auto& deviceList = getTestHardwareEnvironment()->getTestDeviceList();

        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        const auto& testDevice = deviceList[rank % deviceList.size()];

        const DeviceContext& deviceContext = testDevice->deviceContext();
        setActiveDevice(testDevice->deviceInfo());
        const DeviceStream& deviceStream = testDevice->deviceStream();

        FftBackend backend;

        int  numDomainsX;
        int  numDomainsY;
        IVec realGridSizeFull;
        std::tie(realGridSizeFull, numDomainsX, numDomainsY, backend) = param;

        // define local grid sizes - this follows same logic as GROMACS implementation
        std::vector<int> localGridSizesX(numDomainsX);
        for (unsigned int i = 0; i < localGridSizesX.size(); ++i)
        {
            localGridSizesX[i] = ((i + 1) * realGridSizeFull[XX] / numDomainsX)
                                 - (i * realGridSizeFull[XX] / numDomainsX);
            ASSERT_GT(localGridSizesX[i], 0);
        }

        std::vector<int> localGridSizesY(numDomainsY);
        for (unsigned int i = 0; i < localGridSizesY.size(); ++i)
        {
            localGridSizesY[i] = ((i + 1) * realGridSizeFull[YY] / numDomainsY)
                                 - (i * realGridSizeFull[YY] / numDomainsY);
            ASSERT_GT(localGridSizesY[i], 0);
        }

        ivec realGridSize;
        ivec realGridSizePadded;
        ivec complexGridSizePadded;

        // Allocate the device buffers
        DeviceBuffer<float> realGrid, complexGrid;

        const bool     performOutOfPlaceFFT = true;
        const MPI_Comm comm                 = MPI_COMM_WORLD;
        const bool     allocateGrid         = true;
        const int      nz                   = realGridSizeFull[ZZ];
        Gpu3dFft       gpu3dFft(backend,
                          allocateGrid,
                          comm,
                          localGridSizesX,
                          localGridSizesY,
                          nz,
                          performOutOfPlaceFFT,
                          deviceContext,
                          deviceStream,
                          realGridSize,
                          realGridSizePadded,
                          complexGridSizePadded,
                          &realGrid,
                          &complexGrid);

        int sizeInReals = realGridSizePadded[0] * realGridSizePadded[1] * realGridSizePadded[2];

        // initialze random input data
        std::vector<real>                in(sizeInReals);
        std::uniform_real_distribution<> dis(-10.0f, 10.0f);
        std::minstd_rand                 gen(time(NULL) + rank);
        std::generate(in.begin(), in.end(), [&dis, &gen]() {
            // random number between -10 to 10
            return dis(gen);
        });

        // Transfer the real grid input data for the FFT
        copyToDeviceBuffer(
                &realGrid, in.data(), 0, in.size(), deviceStream, GpuApiCallBehavior::Sync, nullptr);

        // Do the forward FFT to compute the complex grid
        CommandEvent* timingEvent = nullptr;
        gpu3dFft.perform3dFft(GMX_FFT_REAL_TO_COMPLEX, timingEvent);

        // clear real grid after the forward FFT, so that we know the
        // final grid is one produced by the complex FFT, not just leftovers
        clearDeviceBufferAsync(&realGrid, 0, sizeInReals, deviceStream);

        // Do the back transform
        gpu3dFft.perform3dFft(GMX_FFT_COMPLEX_TO_REAL, timingEvent);
        deviceStream.synchronize();

        // Transfer the real grid back from the device
        std::vector<float> outputRealGridValues(in.size());
        copyFromDeviceBuffer(outputRealGridValues.data(),
                             &realGrid,
                             0,
                             outputRealGridValues.size(),
                             deviceStream,
                             GpuApiCallBehavior::Sync,
                             nullptr);

        checkRealGrid(realGridSizeFull, realGridSize, realGridSizePadded, in, outputRealGridValues);
    }
};

TEST_P(GpuFftTest3D, GpuFftDecomposition)
{
    GMX_MPI_TEST(RequireRankCount<4>);
    GpuFftTestParams params = GetParam();
    runTest(params);
}

std::vector<GpuFftTestParams> const inputs{
    { IVec{ 5, 6, 9 }, 4, 1, FftBackend::HeFFTe_CUDA }, // slab decomposition
    { IVec{ 5, 6, 9 }, 2, 2, FftBackend::HeFFTe_CUDA }  // pencil decomposition
};

INSTANTIATE_TEST_SUITE_P(GpuFft, GpuFftTest3D, ::testing::ValuesIn(inputs));

} // namespace test
} // namespace gmx
