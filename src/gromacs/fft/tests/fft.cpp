/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014, by the GROMACS development team, led by
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

#include <vector>

#include <gtest/gtest.h>

#include "gromacs/fft/parallel_3dfft.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"

namespace
{

//! Input data for FFT tests.
const real inputdata[] = { //print ",\n".join([",".join(["%4s"%(random.randint(-99,99)/10.,) for i in range(25)]) for j in range(20)])
    -3.5, 6.3, 1.2, 0.3, 1.1, -5.7, 5.8, -1.9, -6.3, -1.4, 7.4, 2.4, -9.9, -7.2, 5.4, 6.1, -1.9, -7.6, 1.4, -3.5, 0.7, 5.6, -4.2, -1.1, -4.4,
    -6.3, -7.2, 4.6, -3.0, -0.9, 7.2, 2.5, -3.6, 6.1, -3.2, -2.1, 6.5, -0.4, -9.0, 2.3, 8.4, 4.0, -5.2, -9.0, 4.7, -3.7, -2.0, -9.5, -3.9, -3.6,
    7.1, 0.8, -0.6, 5.2, -9.3, -4.5, 5.9, 2.2, -5.8, 5.0, 1.2, -0.1, 2.2, 0.2, -7.7, 1.9, -8.4, 4.4, 2.3, -2.9, 6.7, 2.7, 5.8, -3.6, 8.9,
    8.9, 4.3, 9.1, 9.3, -8.7, 4.1, 9.6, -6.2, 6.6, -9.3, 8.2, 4.5, 6.2, 9.4, -8.0, -6.8, -3.3, 7.2, 1.7, 0.6, -4.9, 9.8, 1.3, 3.2, -0.2,
    9.9, 4.4, -9.9, -7.2, 4.4, 4.7, 7.2, -0.3, 0.3, -2.1, 8.4, -2.1, -6.1, 4.1, -5.9, -2.2, -3.8, 5.2, -8.2, -7.8, -8.8, 6.7, -9.5, -4.2, 0.8,
    8.3, 5.2, -9.0, 8.7, 9.8, -9.9, -7.8, -8.3, 9.0, -2.8, -9.2, -9.6, 8.4, 2.5, 6.0, -0.4, 1.3, -0.5, 9.1, -9.5, -0.8, 1.9, -6.2, 4.3, -3.8,
    8.6, -1.9, -2.1, -0.4, -7.1, -3.7, 9.1, -6.4, -0.6, 2.5, 8.0, -5.2, -9.8, -4.3, 4.5, 1.7, 9.3, 9.2, 1.0, 5.3, -4.5, 6.4, -6.6, 3.1, -6.8,
    2.1, 2.0, 7.3, 8.6, 5.0, 5.2, 0.4, -7.1, 4.5, -9.2, -9.1, 0.2, -6.3, -1.1, -9.6, 7.4, -3.7, -5.5, 2.6, -3.5, -0.7, 9.0, 9.8, -8.0, 3.6,
    3.0, -2.2, -2.8, 0.8, 9.0, 2.8, 7.7, -0.7, -5.0, -1.8, -2.3, -0.4, -6.2, -9.1, -9.2, 0.5, 5.7, -3.9, 2.1, 0.6, 0.4, 9.1, 7.4, 7.1, -2.5,
    7.3, 7.8, -4.3, 6.3, -0.8, -3.8, -1.5, 6.6, 2.3, 3.9, -4.6, 5.8, -7.4, 5.9, 2.8, 4.7, 3.9, -5.4, 9.1, -1.6, -1.9, -4.2, -2.6, 0.6, -5.1,
    1.8, 5.2, 4.0, -6.2, 6.5, -9.1, 0.5, 2.1, 7.1, -8.6, 7.6, -9.7, -4.6, -5.7, 6.1, -1.8, -7.3, 9.4, 8.0, -2.6, -1.8, 5.7, 9.3, -7.9, 7.4,
    6.3, 2.0, 9.6, -4.5, -6.2, 6.1, 2.3, 0.8, 5.9, -2.8, -3.5, -1.5, 6.0, -4.9, 3.5, 7.7, -4.2, -9.7, 2.4, 8.1, 5.9, 3.4, -7.5, 7.5, 2.6,
    4.7, 2.7, 2.2, 2.6, 6.2, 7.5, 0.2, -6.4, -2.8, -0.5, -0.3, 0.4, 1.2, 3.5, -4.0, -0.5, 9.3, -7.2, 8.5, -5.5, -1.7, -5.3, 0.3, 3.9, -3.6,
    -3.6, 4.7, -8.1, 1.4, 4.0, 1.3, -4.3, -8.8, -7.3, 6.3, -7.5, -9.0, 9.1, 4.5, -1.9, 1.9, 9.9, -1.7, -9.1, -5.1, 8.5, -9.3, 2.1, -5.8, -3.6,
    -0.8, -0.9, -3.3, -2.7, 7.0, -7.2, -5.0, 7.4, -1.4, 0.0, -4.5, -9.7, 0.7, -1.0, -9.1, -5.3, 4.3, 3.4, -6.6, 9.8, -1.1, 8.9, 5.0, 2.9, 0.2,
    -2.9, 0.8, 6.7, -0.6, 0.6, 4.1, 5.3, -1.7, -0.3, 4.2, 3.7, -8.3, 4.0, 1.3, 6.3, 0.2, 1.3, -1.1, -3.5, 2.8, -7.7, 6.2, -4.9, -9.9, 9.6,
    3.0, -9.2, -8.0, -3.9, 7.9, -6.1, 6.0, 5.9, 9.6, 1.2, 6.2, 3.6, 2.1, 5.8, 9.2, -8.8, 8.8, -3.3, -9.2, 4.6, 1.8, 4.6, 2.9, -2.7, 4.2,
    7.3, -0.4, 7.7, -7.0, 2.1, 0.3, 3.7, 3.3, -8.6, 9.8, 3.6, 3.1, 6.5, -2.4, 7.8, 7.5, 8.4, -2.8, -6.3, -5.1, -2.7, 9.3, -0.8, -9.2, 7.9,
    8.9, 3.4, 0.1, -5.3, -6.8, 4.9, 4.3, -0.7, -2.2, -3.2, -7.5, -2.3, 0.0, 8.1, -9.2, -2.3, -5.7, 2.1, 2.6, 2.0, 0.3, -8.0, -2.0, -7.9, 6.6,
    8.4, 4.0, -6.2, -6.9, -7.2, 7.7, -5.0, 5.3, 1.9, -5.3, -7.5, 8.8, 8.3, 9.0, 8.1, 3.2, 1.2, -5.4, -0.2, 2.1, -5.2, 9.5, 5.9, 5.6, -7.8,
};


class BaseFFTTest : public ::testing::Test
{
    public:
        BaseFFTTest()
            : checker_(data_.rootChecker()), flags_(GMX_FFT_FLAG_CONSERVATIVE)
        {
            // TODO: These tolerances are just something that has been observed
            // to be sufficient to pass the tests.  It would be nicer to
            // actually argue about why they are sufficient (or what is).
            checker_.setDefaultTolerance(
                    gmx::test::relativeToleranceAsPrecisionDependentUlp(10.0, 64, 512));
        }
        ~BaseFFTTest()
        {
            gmx_fft_cleanup();
        }

        gmx::test::TestReferenceData    data_;
        gmx::test::TestReferenceChecker checker_;
        std::vector<real>               in_, out_;
        int                             flags_;
};

class FFTTest : public BaseFFTTest
{
    public:
        FFTTest() : fft_(NULL)
        {
        }
        ~FFTTest()
        {
            if (fft_)
            {
                gmx_fft_destroy(fft_);
            }
        }
        gmx_fft_t fft_;
};

class ManyFFTTest : public BaseFFTTest
{
    public:
        ManyFFTTest() : fft_(NULL)
        {
        }
        ~ManyFFTTest()
        {
            if (fft_)
            {
                gmx_many_fft_destroy(fft_);
            }
        }
        gmx_fft_t fft_;
};


//TODO: Add tests for aligned/not-aligned input/output memory

class FFTTest1D : public FFTTest, public ::testing::WithParamInterface<int>
{

};

class FFFTest3D : public BaseFFTTest
{
    public:
        FFFTest3D() : fft_(NULL)
        {
        }
        ~FFFTest3D()
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
    ASSERT_LE(nx*2, static_cast<int>(sizeof(inputdata)/sizeof(real)));

    in_  = std::vector<real>(inputdata, inputdata+nx*2);
    out_ = std::vector<real>(nx*2);
    real* in  = &in_[0];
    real* out = &out_[0];

    gmx_fft_init_1d(&fft_, nx, flags_);

    gmx_fft_1d(fft_, GMX_FFT_FORWARD, in, out);
    checker_.checkSequenceArray(nx*2, out, "forward");
    gmx_fft_1d(fft_, GMX_FFT_BACKWARD, in, out);
    checker_.checkSequenceArray(nx*2, out, "backward");
}

TEST_P(FFTTest1D, Real)
{
    const int rx = GetParam();
    const int cx = (rx/2+1);
    ASSERT_LE(cx*2, static_cast<int>(sizeof(inputdata)/sizeof(real)));

    in_  = std::vector<real>(inputdata, inputdata+cx*2);
    out_ = std::vector<real>(cx*2);
    real* in  = &in_[0];
    real* out = &out_[0];

    gmx_fft_init_1d_real(&fft_, rx, flags_);

    gmx_fft_1d_real(fft_, GMX_FFT_REAL_TO_COMPLEX, in, out);
    checker_.checkSequenceArray(cx*2, out, "forward");
    gmx_fft_1d_real(fft_, GMX_FFT_COMPLEX_TO_REAL, in, out);
    checker_.checkSequenceArray(rx, out, "backward");
}

INSTANTIATE_TEST_CASE_P(7_8_25_36_60,
                        FFTTest1D, ::testing::Values(7, 8, 25, 36, 60));


TEST_F(ManyFFTTest, Complex1DLength48Multi5Test)
{
    const int nx = 48;
    const int N  = 5;

    in_  = std::vector<real>(inputdata, inputdata+nx*2*N);
    out_ = std::vector<real>(nx*2*N);
    real* in  = &in_[0];
    real* out = &out_[0];

    gmx_fft_init_many_1d(&fft_, nx, N, flags_);

    gmx_fft_many_1d(fft_, GMX_FFT_FORWARD, in, out);
    checker_.checkSequenceArray(nx*2*N, out, "forward");
    gmx_fft_many_1d(fft_, GMX_FFT_BACKWARD, in, out);
    checker_.checkSequenceArray(nx*2*N, out, "backward");
}

TEST_F(ManyFFTTest, Real1DLength48Multi5Test)
{
    const int rx = 48;
    const int cx = (rx/2+1);
    const int N  = 5;

    in_  = std::vector<real>(inputdata, inputdata+cx*2*N);
    out_ = std::vector<real>(cx*2*N);
    real* in  = &in_[0];
    real* out = &out_[0];

    gmx_fft_init_many_1d_real(&fft_, rx, N, flags_);

    gmx_fft_many_1d_real(fft_, GMX_FFT_REAL_TO_COMPLEX, in, out);
    checker_.checkSequenceArray(cx*2*N, out, "forward");
    gmx_fft_many_1d_real(fft_, GMX_FFT_COMPLEX_TO_REAL, in, out);
    checker_.checkSequenceArray(rx*N, out, "backward");
}

TEST_F(FFTTest, Real2DLength18_15Test)
{
    const int rx = 18;
    const int cx = (rx/2+1);
    const int ny = 15;

    in_  = std::vector<real>(inputdata, inputdata+cx*2*ny);
    out_ = std::vector<real>(cx*2*ny);
    real* in  = &in_[0];
    real* out = &out_[0];

    gmx_fft_init_2d_real(&fft_, rx, ny, flags_);

    gmx_fft_2d_real(fft_, GMX_FFT_REAL_TO_COMPLEX, in, out);
    checker_.checkSequenceArray(cx*2*ny, out, "forward");
//    known to be wrong for gmx_fft_mkl. And not used.
//    gmx_fft_2d_real(_fft,GMX_FFT_COMPLEX_TO_REAL,in,out);
//    _checker.checkSequenceArray(rx*ny, out, "backward");
}

//TODO: test with threads and more than 1 MPI ranks
TEST_F(FFFTest3D, Real5_6_9)
{
    int        ndata[] = {5, 6, 9};
    MPI_Comm   comm[]  = {MPI_COMM_NULL, MPI_COMM_NULL};
    real     * rdata;
    t_complex* cdata;
    ivec       local_ndata, offset, rsize, csize, complex_order;

    gmx_parallel_3dfft_init(&fft_, ndata, &rdata, &cdata,
                            comm, TRUE, 1);

    gmx_parallel_3dfft_real_limits(fft_, local_ndata, offset, rsize);
    gmx_parallel_3dfft_complex_limits(fft_, complex_order,
                                      local_ndata, offset, csize);
    checker_.checkVector(rsize, "rsize");
    checker_.checkVector(csize, "csize");
    int size = csize[0]*csize[1]*csize[2];

    memcpy(rdata, inputdata, size*sizeof(t_complex));
    gmx_parallel_3dfft_execute(fft_, GMX_FFT_REAL_TO_COMPLEX, 0, NULL);
    //TODO use std::complex and add checkComplex for it
    checker_.checkSequenceArray(size*2,
                                reinterpret_cast<real*>(cdata), "forward");

    memcpy(cdata, inputdata, size*sizeof(t_complex));
    gmx_parallel_3dfft_execute(fft_, GMX_FFT_COMPLEX_TO_REAL, 0, NULL);
    for (int i = 0; i < ndata[0]*ndata[1]; i++) //check sequence but skip unused data
    {
        checker_.checkSequenceArray(ndata[2], rdata+i*rsize[2],
                                    gmx::formatString("backward %d", i).c_str());
    }
}

} // namespace
