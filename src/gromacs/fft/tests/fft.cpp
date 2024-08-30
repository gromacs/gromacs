/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2012- The GROMACS Authors
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
 * Current reference data is generated in double precision using the Reference
 * build type, except for the compiler (Apple Clang).
 *
 * \author Roland Schulz <roland@utk.edu>
 * \ingroup module_fft
 */
#include "gmxpre.h"

#include "gromacs/fft/fft.h"

#include "config.h"

#include <cstring>

#include <algorithm>
#include <optional>
#include <string>
#include <tuple>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/fft/gpu_3dfft.h"
#include "gromacs/fft/parallel_3dfft.h"
#include "gromacs/gpu_utils/clfftinitializer.h"
#include "gromacs/math/gmxcomplex.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"
#if GMX_GPU
#    include "gromacs/gpu_utils/devicebuffer.h"
#endif
#include "gromacs/utility/strconvert.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/naming.h"
#include "testutils/refdata.h"
#include "testutils/test_hardware_environment.h"
#include "testutils/testasserts.h"
#include "testutils/testmatchers.h"

namespace gmx
{
namespace test
{
namespace
{

/*! \brief Input data for FFT tests.
 *
 * TODO If we require compilers that all support C++11 user literals,
 * then this array could be of type real, initialized with e.g. -3.5_r
 * that does not suffer from implicit narrowing with brace
 * initializers, and we would not have to do so much useless copying
 * during the unit tests below.
 */
const double inputdata[500] = {
    // print ",\n".join([",".join(["%4s"%(random.randint(-99,99)/10.,) for i in range(25)]) for j in range(20)])
    -3.5, 6.3,  1.2,  0.3,  1.1,  -5.7, 5.8,  -1.9, -6.3, -1.4, 7.4,  2.4,  -9.9, -7.2, 5.4,  6.1,
    -1.9, -7.6, 1.4,  -3.5, 0.7,  5.6,  -4.2, -1.1, -4.4, -6.3, -7.2, 4.6,  -3.0, -0.9, 7.2,  2.5,
    -3.6, 6.1,  -3.2, -2.1, 6.5,  -0.4, -9.0, 2.3,  8.4,  4.0,  -5.2, -9.0, 4.7,  -3.7, -2.0, -9.5,
    -3.9, -3.6, 7.1,  0.8,  -0.6, 5.2,  -9.3, -4.5, 5.9,  2.2,  -5.8, 5.0,  1.2,  -0.1, 2.2,  0.2,
    -7.7, 1.9,  -8.4, 4.4,  2.3,  -2.9, 6.7,  2.7,  5.8,  -3.6, 8.9,  8.9,  4.3,  9.1,  9.3,  -8.7,
    4.1,  9.6,  -6.2, 6.6,  -9.3, 8.2,  4.5,  6.2,  9.4,  -8.0, -6.8, -3.3, 7.2,  1.7,  0.6,  -4.9,
    9.8,  1.3,  3.2,  -0.2, 9.9,  4.4,  -9.9, -7.2, 4.4,  4.7,  7.2,  -0.3, 0.3,  -2.1, 8.4,  -2.1,
    -6.1, 4.1,  -5.9, -2.2, -3.8, 5.2,  -8.2, -7.8, -8.8, 6.7,  -9.5, -4.2, 0.8,  8.3,  5.2,  -9.0,
    8.7,  9.8,  -9.9, -7.8, -8.3, 9.0,  -2.8, -9.2, -9.6, 8.4,  2.5,  6.0,  -0.4, 1.3,  -0.5, 9.1,
    -9.5, -0.8, 1.9,  -6.2, 4.3,  -3.8, 8.6,  -1.9, -2.1, -0.4, -7.1, -3.7, 9.1,  -6.4, -0.6, 2.5,
    8.0,  -5.2, -9.8, -4.3, 4.5,  1.7,  9.3,  9.2,  1.0,  5.3,  -4.5, 6.4,  -6.6, 3.1,  -6.8, 2.1,
    2.0,  7.3,  8.6,  5.0,  5.2,  0.4,  -7.1, 4.5,  -9.2, -9.1, 0.2,  -6.3, -1.1, -9.6, 7.4,  -3.7,
    -5.5, 2.6,  -3.5, -0.7, 9.0,  9.8,  -8.0, 3.6,  3.0,  -2.2, -2.8, 0.8,  9.0,  2.8,  7.7,  -0.7,
    -5.0, -1.8, -2.3, -0.4, -6.2, -9.1, -9.2, 0.5,  5.7,  -3.9, 2.1,  0.6,  0.4,  9.1,  7.4,  7.1,
    -2.5, 7.3,  7.8,  -4.3, 6.3,  -0.8, -3.8, -1.5, 6.6,  2.3,  3.9,  -4.6, 5.8,  -7.4, 5.9,  2.8,
    4.7,  3.9,  -5.4, 9.1,  -1.6, -1.9, -4.2, -2.6, 0.6,  -5.1, 1.8,  5.2,  4.0,  -6.2, 6.5,  -9.1,
    0.5,  2.1,  7.1,  -8.6, 7.6,  -9.7, -4.6, -5.7, 6.1,  -1.8, -7.3, 9.4,  8.0,  -2.6, -1.8, 5.7,
    9.3,  -7.9, 7.4,  6.3,  2.0,  9.6,  -4.5, -6.2, 6.1,  2.3,  0.8,  5.9,  -2.8, -3.5, -1.5, 6.0,
    -4.9, 3.5,  7.7,  -4.2, -9.7, 2.4,  8.1,  5.9,  3.4,  -7.5, 7.5,  2.6,  4.7,  2.7,  2.2,  2.6,
    6.2,  7.5,  0.2,  -6.4, -2.8, -0.5, -0.3, 0.4,  1.2,  3.5,  -4.0, -0.5, 9.3,  -7.2, 8.5,  -5.5,
    -1.7, -5.3, 0.3,  3.9,  -3.6, -3.6, 4.7,  -8.1, 1.4,  4.0,  1.3,  -4.3, -8.8, -7.3, 6.3,  -7.5,
    -9.0, 9.1,  4.5,  -1.9, 1.9,  9.9,  -1.7, -9.1, -5.1, 8.5,  -9.3, 2.1,  -5.8, -3.6, -0.8, -0.9,
    -3.3, -2.7, 7.0,  -7.2, -5.0, 7.4,  -1.4, 0.0,  -4.5, -9.7, 0.7,  -1.0, -9.1, -5.3, 4.3,  3.4,
    -6.6, 9.8,  -1.1, 8.9,  5.0,  2.9,  0.2,  -2.9, 0.8,  6.7,  -0.6, 0.6,  4.1,  5.3,  -1.7, -0.3,
    4.2,  3.7,  -8.3, 4.0,  1.3,  6.3,  0.2,  1.3,  -1.1, -3.5, 2.8,  -7.7, 6.2,  -4.9, -9.9, 9.6,
    3.0,  -9.2, -8.0, -3.9, 7.9,  -6.1, 6.0,  5.9,  9.6,  1.2,  6.2,  3.6,  2.1,  5.8,  9.2,  -8.8,
    8.8,  -3.3, -9.2, 4.6,  1.8,  4.6,  2.9,  -2.7, 4.2,  7.3,  -0.4, 7.7,  -7.0, 2.1,  0.3,  3.7,
    3.3,  -8.6, 9.8,  3.6,  3.1,  6.5,  -2.4, 7.8,  7.5,  8.4,  -2.8, -6.3, -5.1, -2.7, 9.3,  -0.8,
    -9.2, 7.9,  8.9,  3.4,  0.1,  -5.3, -6.8, 4.9,  4.3,  -0.7, -2.2, -3.2, -7.5, -2.3, 0.0,  8.1,
    -9.2, -2.3, -5.7, 2.1,  2.6,  2.0,  0.3,  -8.0, -2.0, -7.9, 6.6,  8.4,  4.0,  -6.2, -6.9, -7.2,
    7.7,  -5.0, 5.3,  1.9,  -5.3, -7.5, 8.8,  8.3,  9.0,  8.1,  3.2,  1.2,  -5.4, -0.2, 2.1,  -5.2,
    9.5,  5.9,  5.6,  -7.8,
};


class BaseFFTTest : public ::testing::Test
{
public:
    BaseFFTTest() : flags_(GMX_FFT_FLAG_CONSERVATIVE) {}
    ~BaseFFTTest() override { gmx_fft_cleanup(); }

    TestReferenceData data_;
    std::vector<real> in_, out_;
    int               flags_;
    // TODO: These tolerances are just something that has been observed
    // to be sufficient to pass the tests.  It would be nicer to
    // actually argue about why they are sufficient (or what is).
    // Should work for both one-way and forward+backward transform.
    FloatingPointTolerance defaultTolerance_ = relativeToleranceAsPrecisionDependentUlp(10.0, 64, 512);
};

class FFTTest : public BaseFFTTest
{
public:
    FFTTest() : fft_(nullptr) { checker_.setDefaultTolerance(defaultTolerance_); }
    ~FFTTest() override
    {
        if (fft_)
        {
            gmx_fft_destroy(fft_);
        }
    }
    TestReferenceChecker checker_ = data_.rootChecker();
    gmx_fft_t            fft_;
};

class ManyFFTTest : public BaseFFTTest
{
public:
    ManyFFTTest() : fft_(nullptr) { checker_.setDefaultTolerance(defaultTolerance_); }
    ~ManyFFTTest() override
    {
        if (fft_)
        {
            gmx_many_fft_destroy(fft_);
        }
    }
    TestReferenceChecker checker_ = data_.rootChecker();
    gmx_fft_t            fft_;
};


// TODO: Add tests for aligned/not-aligned input/output memory

class FFTTest1D : public FFTTest, public ::testing::WithParamInterface<int>
{
};

class FFTTest3D : public BaseFFTTest
{
public:
    FFTTest3D() : fft_(nullptr) {}
    ~FFTTest3D() override
    {
        if (fft_)
        {
            gmx_parallel_3dfft_destroy(fft_);
        }
    }
    gmx_parallel_3dfft_t fft_;
};


TEST_P(FFTTest1D, Complex)
{
    const int nx = GetParam();
    ASSERT_LE(nx * 2, static_cast<int>(sizeof(inputdata) / sizeof(inputdata[0])));

    in_ = std::vector<real>(nx * 2);
    std::copy(inputdata, inputdata + nx * 2, in_.begin());
    out_      = std::vector<real>(nx * 2);
    real* in  = &in_[0];
    real* out = &out_[0];

    gmx_fft_init_1d(&fft_, nx, flags_);

    gmx_fft_1d(fft_, GMX_FFT_FORWARD, in, out);
    checker_.checkSequenceArray(nx * 2, out, "forward");
    gmx_fft_1d(fft_, GMX_FFT_BACKWARD, in, out);
    checker_.checkSequenceArray(nx * 2, out, "backward");
}

TEST_P(FFTTest1D, Real)
{
    const int rx = GetParam();
    const int cx = (rx / 2 + 1);
    ASSERT_LE(cx * 2, static_cast<int>(sizeof(inputdata) / sizeof(inputdata[0])));

    in_ = std::vector<real>(cx * 2);
    std::copy(inputdata, inputdata + cx * 2, in_.begin());
    out_      = std::vector<real>(cx * 2);
    real* in  = &in_[0];
    real* out = &out_[0];

    gmx_fft_init_1d_real(&fft_, rx, flags_);

    gmx_fft_1d_real(fft_, GMX_FFT_REAL_TO_COMPLEX, in, out);
    checker_.checkSequenceArray(cx * 2, out, "forward");
    gmx_fft_1d_real(fft_, GMX_FFT_COMPLEX_TO_REAL, in, out);
    checker_.checkSequenceArray(rx, out, "backward");
}

INSTANTIATE_TEST_SUITE_P(7_8_25_36_60, FFTTest1D, ::testing::Values(7, 8, 25, 36, 60));


TEST_F(ManyFFTTest, Complex1DLength48Multi5Test)
{
    const int nx = 48;
    const int N  = 5;

    in_ = std::vector<real>(nx * 2 * N);
    std::copy(inputdata, inputdata + nx * 2 * N, in_.begin());
    out_      = std::vector<real>(nx * 2 * N);
    real* in  = &in_[0];
    real* out = &out_[0];

    gmx_fft_init_many_1d(&fft_, nx, N, flags_);

    gmx_fft_many_1d(fft_, GMX_FFT_FORWARD, in, out);
    checker_.checkSequenceArray(nx * 2 * N, out, "forward");
    gmx_fft_many_1d(fft_, GMX_FFT_BACKWARD, in, out);
    checker_.checkSequenceArray(nx * 2 * N, out, "backward");
}

TEST_F(ManyFFTTest, Real1DLength48Multi5Test)
{
    const int rx = 48;
    const int cx = (rx / 2 + 1);
    const int N  = 5;

    in_ = std::vector<real>(cx * 2 * N);
    std::copy(inputdata, inputdata + cx * 2 * N, in_.begin());
    out_      = std::vector<real>(cx * 2 * N);
    real* in  = &in_[0];
    real* out = &out_[0];

    gmx_fft_init_many_1d_real(&fft_, rx, N, flags_);

    gmx_fft_many_1d_real(fft_, GMX_FFT_REAL_TO_COMPLEX, in, out);
    checker_.checkSequenceArray(cx * 2 * N, out, "forward");
    gmx_fft_many_1d_real(fft_, GMX_FFT_COMPLEX_TO_REAL, in, out);
    checker_.checkSequenceArray(rx * N, out, "backward");
}

TEST_F(FFTTest, Real2DLength18_15Test)
{
    const int rx = 18;
    const int cx = (rx / 2 + 1);
    const int ny = 15;

    in_ = std::vector<real>(cx * 2 * ny);
    std::copy(inputdata, inputdata + cx * 2 * ny, in_.begin());
    out_      = std::vector<real>(cx * 2 * ny);
    real* in  = &in_[0];
    real* out = &out_[0];

    gmx_fft_init_2d_real(&fft_, rx, ny, flags_);

    gmx_fft_2d_real(fft_, GMX_FFT_REAL_TO_COMPLEX, in, out);
    checker_.checkSequenceArray(cx * 2 * ny, out, "forward");
    //    known to be wrong for gmx_fft_mkl. And not used.
    //    gmx_fft_2d_real(_fft,GMX_FFT_COMPLEX_TO_REAL,in,out);
    //    _checker.checkSequenceArray(rx*ny, out, "backward");
}

using FFTTest3DParameters = std::tuple<int, int, int>;

//! Tuple of formatters to name the parameterized test cases
const test::NameOfTestFromTuple<FFTTest3DParameters> sc_testNamer{
    std::make_tuple(intToString, intToString, intToString)
};

class ParameterizedFFTTest3D : public FFTTest3D, public ::testing::WithParamInterface<FFTTest3DParameters>
{
};

/*! \brief Check that the real grid after forward and backward
 * 3D transforms matches the input real grid. */
void checkRealGrid(const ivec           realGridSize,
                   const ivec           realGridSizePadded,
                   ArrayRef<const real> inputRealGrid,
                   ArrayRef<real>       outputRealGridValues)
{
    SCOPED_TRACE("Checking the real grid");
    // Normalize the output (as the implementation does not
    // normalize either FFT)
    const real normalizationConstant = 1.0 / (realGridSize[XX] * realGridSize[YY] * realGridSize[ZZ]);
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

// TODO: test with threads and more than 1 MPI ranks
TEST_P(ParameterizedFFTTest3D, RunsOnHost)
{
    ivec realGridSize = { std::get<0>(GetParam()), std::get<1>(GetParam()), std::get<2>(GetParam()) };
    MPI_Comm   comm[] = { MPI_COMM_NULL, MPI_COMM_NULL };
    real*      rdata;
    t_complex* cdata;
    ivec       local_ndata, offset, realGridSizePadded, complexGridSizePadded, complex_order;

    gmx_parallel_3dfft_init(&fft_, realGridSize, &rdata, &cdata, comm, TRUE, 1);

    gmx_parallel_3dfft_real_limits(fft_, local_ndata, offset, realGridSizePadded);
    gmx_parallel_3dfft_complex_limits(fft_, complex_order, local_ndata, offset, complexGridSizePadded);
    // Check complex grid against reference data only when it
    // exists. The real grid is always checked by doing the
    // identity transform via forward and inverse transform.
    std::optional<TestReferenceChecker> checker;
    if ((realGridSize[XX] == 5 && realGridSize[YY] == 6 && realGridSize[ZZ] == 9)
        || (realGridSize[XX] == 5 && realGridSize[YY] == 5 && realGridSize[ZZ] == 10))
    {
        checker = TestReferenceChecker(data_.rootChecker()); // Must be inside the loop to avoid warnings
        checker.value().setDefaultTolerance(defaultTolerance_);
        checker.value().checkVector(realGridSizePadded, "realGridSizePadded");
        checker.value().checkVector(complexGridSizePadded, "complexGridSizePadded");
    }
    int size = complexGridSizePadded[0] * complexGridSizePadded[1] * complexGridSizePadded[2];
    int sizeInBytes = size * sizeof(t_complex);
    int sizeInReals = sizeInBytes / sizeof(real);

    // Prepare the real grid
    in_.resize(sizeInReals);
    ASSERT_LT(sizeInReals, sizeof(inputdata) / sizeof(inputdata[0]));
    // Use std::copy to convert from double to real easily
    std::copy(inputdata, inputdata + sizeInReals, in_.begin());
    // Use memcpy to convert to t_complex easily
    memcpy(rdata, in_.data(), sizeInBytes);

    // Do the forward FFT to compute the complex grid
    gmx_parallel_3dfft_execute(fft_, GMX_FFT_REAL_TO_COMPLEX, 0, nullptr);

    // Check the complex grid (NB this data has not been normalized)
    ArrayRef<real> complexGridValues = arrayRefFromArray(reinterpret_cast<real*>(cdata), size * 2);
    if (checker.has_value())
    {
        checker.value().checkSequence(
                complexGridValues.begin(), complexGridValues.end(), "ComplexGridAfterRealToComplex");
    }

    // Do the back transform
    gmx_parallel_3dfft_execute(fft_, GMX_FFT_COMPLEX_TO_REAL, 0, nullptr);

    ArrayRef<real> outputRealGridValues = arrayRefFromArray(
            rdata, realGridSizePadded[XX] * realGridSizePadded[YY] * realGridSizePadded[ZZ]);
    checkRealGrid(realGridSize, realGridSizePadded, in_, outputRealGridValues);
}

#if GMX_GPU

/*! \brief Whether the FFT is in- or out-of-place
 *
 * For real simulations, we use in-place FFT with BBFFT (see pme_gpu_init_internal)
 * and out-of-place FFT otherwise, so that's what we test.
 */
constexpr bool sc_performOutOfPlaceFFT = (GMX_GPU_FFT_BBFFT == 0);

/*! \brief Return the output grid depending on whether in- or out-of
 * place FFT is used
 *
 * Some versions of clang complain of unused code if we would just
 * branch on the value of sc_performOutOfPlaceFFT at run time, because
 * in any single configuration there would indeed be unused code. So
 * the two template specializations are needed so that the compiler
 * only compiles the template that is used. */
template<bool performOutOfPlaceFFT>
DeviceBuffer<float>* actualOutputGrid(DeviceBuffer<float>* realGrid, DeviceBuffer<float>* complexGrid);

#    if GMX_SYCL_DPCPP && (GMX_GPU_FFT_BBFFT)

template<>
DeviceBuffer<float>* actualOutputGrid<false>(DeviceBuffer<float>* realGrid,
                                             DeviceBuffer<float>* /* complexGrid */)
{
    return realGrid;
};

#    else

template<>
DeviceBuffer<float>* actualOutputGrid<true>(DeviceBuffer<float>* /* realGrid */, DeviceBuffer<float>* complexGrid)
{
    return complexGrid;
}

#    endif

// This is not the same test case as RunsOnHost because the
// dimension ordering of the complex grid differs for solve
// when run on the GPU vs CPU, because the latter caters to a 5D
// pencil decomposition when using multiple PME ranks.
// On the CPU, the 3DFFT code handles
// the transpose from XYZ to YXZ and back, however on the GPU that
// is part of spread and gather routines. So the reference data
// is in a different order in these cases that test the 3DFFT output.
TEST_P(ParameterizedFFTTest3D, RunsOnDevices)
{
    // Ensure library resources are managed appropriately
    ClfftInitializer clfftInitializer;
    for (const auto& testDevice : getTestHardwareEnvironment()->getTestDeviceList())
    {
        const DeviceContext& deviceContext = testDevice->deviceContext();
        const DeviceStream&  deviceStream  = testDevice->deviceStream();
        deviceContext.activate();

        ivec realGridSize = { std::get<0>(GetParam()), std::get<1>(GetParam()), std::get<2>(GetParam()) };
        // Note the real-grid padding differs from that on the CPU
        ivec realGridSizePadded = { realGridSize[XX], realGridSize[YY], ((realGridSize[ZZ] / 2) + 1) * 2 };
        ivec complexGridSizePadded = { realGridSize[XX], realGridSize[YY], (realGridSize[ZZ] / 2) + 1 };

        // Check complex grid against reference data only when it
        // exists. The real grid is always checked by doing the
        // identity transform via forward and inverse transform.
        std::optional<TestReferenceChecker> checker;
        if ((realGridSize[XX] == 5 && realGridSize[YY] == 6 && realGridSize[ZZ] == 9)
            || (realGridSize[XX] == 5 && realGridSize[YY] == 5 && realGridSize[ZZ] == 10))
        {
            checker = TestReferenceChecker(data_.rootChecker()); // Must be inside the loop to avoid warnings
            checker.value().setDefaultTolerance(defaultTolerance_);
            checker.value().checkVector(realGridSizePadded, "realGridSizePadded");
            checker.value().checkVector(complexGridSizePadded, "complexGridSizePadded");
        }

        int size = complexGridSizePadded[0] * complexGridSizePadded[1] * complexGridSizePadded[2];
        int sizeInReals = size * 2;
        GMX_RELEASE_ASSERT(sizeof(inputdata) / sizeof(inputdata[0]) >= size_t(sizeInReals),
                           "Size of inputdata is too small");

        // Set up the complex grid. Complex numbers take twice the
        // memory.
        std::vector<float> complexGridValues(sizeInReals);

        // Prepare the real grid
        in_.resize(sizeInReals);
        ASSERT_LT(sizeInReals, sizeof(inputdata) / sizeof(inputdata[0]));
        // Use std::copy to convert from double to real easily
        std::copy(inputdata, inputdata + sizeInReals, in_.begin());

#    if GMX_GPU_CUDA
        const FftBackend backend = FftBackend::Cufft;
#    elif GMX_GPU_HIP
#        if GMX_GPU_FFT_VKFFT
        const FftBackend backend = FftBackend::HipVkfft;
#        else
        const FftBackend backend = FftBackend::Hipfft;
#        endif
#    elif GMX_GPU_OPENCL
#        if GMX_GPU_FFT_VKFFT
        const FftBackend backend = FftBackend::OclVkfft;
#        else
        const FftBackend backend = FftBackend::Ocl;
#        endif
#    elif GMX_GPU_SYCL
#        if GMX_GPU_FFT_MKL
        const FftBackend backend = FftBackend::SyclMkl;
#        elif GMX_GPU_FFT_ONEMKL
        const FftBackend backend = FftBackend::SyclOneMkl;
#        elif GMX_GPU_FFT_BBFFT
        const FftBackend backend = FftBackend::SyclBbfft;
#        elif GMX_GPU_FFT_ROCFFT
        const FftBackend backend = FftBackend::SyclRocfft;
#        elif GMX_GPU_FFT_VKFFT
        const FftBackend backend = FftBackend::SyclVkfft;
#        else
        // Use stub backend so compilation succeeds
        const FftBackend backend = FftBackend::Sycl;
        // Don't complain about unused reference data
        if (checker.has_value())
        {
            checker.value().disableUnusedEntriesCheck();
        }
        // Skip the rest of the test
        GTEST_SKIP() << "No supported GPU FFT backend detected";
#        endif
#    endif

        SCOPED_TRACE("Allocating the device buffers");
        DeviceBuffer<float> realGrid, complexGrid = nullptr;
        allocateDeviceBuffer(&realGrid, in_.size(), deviceContext);

        MPI_Comm           comm                    = MPI_COMM_NULL;
        const bool         allocateGrid            = false;
        std::array<int, 1> gridSizesInXForEachRank = { 0 };
        std::array<int, 1> gridSizesInYForEachRank = { 0 };
        const int          nz                      = realGridSize[ZZ];
        Gpu3dFft           gpu3dFft(backend,
                          allocateGrid,
                          comm,
                          gridSizesInXForEachRank,
                          gridSizesInYForEachRank,
                          nz,
                          sc_performOutOfPlaceFFT,
                          deviceContext,
                          deviceStream,
                          realGridSize,
                          realGridSizePadded,
                          complexGridSizePadded,
                          &realGrid,
                          actualOutputGrid<sc_performOutOfPlaceFFT>(&realGrid, &complexGrid));

        // Transfer the real grid input data for the FFT
        copyToDeviceBuffer(
                &realGrid, in_.data(), 0, in_.size(), deviceStream, GpuApiCallBehavior::Sync, nullptr);

        // Do the forward FFT to compute the complex grid
        CommandEvent* timingEvent = nullptr;
        gpu3dFft.perform3dFft(GMX_FFT_REAL_TO_COMPLEX, timingEvent);
        deviceStream.synchronize();

        SCOPED_TRACE("Checking the complex grid - NB this data has not been normalized");
        copyFromDeviceBuffer(complexGridValues.data(),
                             actualOutputGrid<sc_performOutOfPlaceFFT>(&realGrid, &complexGrid),
                             0,
                             complexGridValues.size(),
                             deviceStream,
                             GpuApiCallBehavior::Sync,
                             nullptr);
        if (checker.has_value())
        {
            checker.value().checkSequence(complexGridValues.begin(),
                                          complexGridValues.end(),
                                          "ComplexGridAfterRealToComplex");
        }

        std::vector<float> outputRealGridValues(in_.size());
        if (sc_performOutOfPlaceFFT)
        {
            // Clear the real grid input data for the FFT so we can
            // compute the back transform into it and observe that it did
            // the work expected.
            copyToDeviceBuffer(&realGrid,
                               outputRealGridValues.data(),
                               0,
                               outputRealGridValues.size(),
                               deviceStream,
                               GpuApiCallBehavior::Sync,
                               nullptr);
        }

        SCOPED_TRACE("Doing the back transform");
        gpu3dFft.perform3dFft(GMX_FFT_COMPLEX_TO_REAL, timingEvent);
        deviceStream.synchronize();

        // Transfer the real grid back from the device
        copyFromDeviceBuffer(outputRealGridValues.data(),
                             &realGrid,
                             0,
                             outputRealGridValues.size(),
                             deviceStream,
                             GpuApiCallBehavior::Sync,
                             nullptr);

        checkRealGrid(realGridSize, realGridSizePadded, in_, outputRealGridValues);

        SCOPED_TRACE("Cleaning up");
        freeDeviceBuffer(&realGrid);
    }
}

#endif // GMX_GPU

/* Uncomment this to run more exhaustive tests
INSTANTIATE_TEST_SUITE_P(ScanWorks,
                         ParameterizedFFTTest3D,
                         ::testing::Combine(::testing::Range(4, 8, 1),
                                            ::testing::Range(4, 8, 1),
                                            ::testing::Range(4, 8, 1)),
                         sc_testNamer);
*/

INSTANTIATE_TEST_SUITE_P(Works,
                         ParameterizedFFTTest3D,
                         ::testing::Values(FFTTest3DParameters{ 5, 6, 9 }, FFTTest3DParameters{ 5, 5, 10 }),
                         sc_testNamer);

} // namespace
} // namespace test
} // namespace gmx
