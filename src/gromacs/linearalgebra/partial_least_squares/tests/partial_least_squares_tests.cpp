/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
   \brief
   Tests for partial least squares function

   \author Berenger Bramas <berenger.bramas@mpcdf.mpg.de>
 */
#include "gmxpre.h"

#include "gromacs/linearalgebra/partial_least_squares/partial_least_squares.h"

#include <cmath>

#include <vector>

#include <gtest/gtest.h>

#include "gromacs/math/irreg_array_math.h"

namespace gmx
{

namespace partial_least_squares_tests
{

namespace
{

//! \brief The matrix type to use in the following tests
using Pls2dArray = gmx::IrregArray2D<real, typename gmx::SimdSetup<real>::allocator_type>;
//! \brief The vector type to use in the following tests
using Pls1dArray = std::vector<real, typename gmx::SimdSetup<real>::allocator_type>;

TEST(PartialLeastSquaresTest, SimpleTest){
    using allocator_type = gmx::SimdSetup<real>::allocator_type;
    using index_type     = Pls2dArray::index_type;

    const index_type                               n = 4;
    const index_type                               p = 4;
    const index_type                               k = 4;
    std::vector<real, allocator_type>              y(n, 0);
    gmx::IrregArray2D<real, allocator_type>        X(index_type(0), index_type(n-1), index_type(0), index_type(p-1));

    gmx::IrregArray2D<real, allocator_type>        W(index_type(0), index_type(p-1), index_type(0), index_type(k-1));
    std::vector<real, allocator_type>              q(k);

    // We set a shape like:
    // y = [1 2 3 4]'
    // X = [1 1 0 0 ; 0 1 1 0 ; 0 0 1 1; 0 0 0 1]
    for (index_type idx_n = 0; idx_n < n; ++idx_n)
    {
        if (idx_n < p)
        {
            X(idx_n, idx_n) = 1;
            if (idx_n + 1 < p)
            {
                X(idx_n, idx_n+1) = 1;
            }
        }
        y[idx_n] = idx_n+1;
    }

    PartialLeastSquares<Pls2dArray, Pls1dArray>().pls_denham(X, y, n, p, k, &W, &q);

    std::vector<real, allocator_type>       t(p, 0);
    for (index_type i = 0; i < p; ++i)
    {
        for (index_type ik = 0; ik < k; ++ik)
        {
            t[i] += W(i, ik) * q[ik];
        }
    }

    std::vector<real, allocator_type>       Xt(n, 0);
    for (index_type i = 0; i < n; ++i)
    {
        for (index_type ip = 0; ip < p; ++ip)
        {
            Xt[i] += X(i, ip) * t[ip];
        }
        // We know from the input that Xt and y must be close
        EXPECT_NEAR(Xt[i], y[i], 2e-6);
    }

    for (index_type i = 0; i < n; ++i)
    {
        // The objective of the function is to find such that
        // the difference between y and X t is low
        EXPECT_NEAR(std::abs(y[i] - Xt[i]), 0, 2e-6);
    }
}

TEST(PartialLeastSquaresTest, SimpleTestN){
    using allocator_type = gmx::SimdSetup<real>::allocator_type;
    using index_type     = Pls2dArray::index_type;

    const index_type                               n = 5;
    const index_type                               p = 4;
    const index_type                               k = 4;
    std::vector<real, allocator_type>              y(n, 0);
    gmx::IrregArray2D<real, allocator_type>        X(index_type(0), index_type(n-1), index_type(0), index_type(p-1));

    gmx::IrregArray2D<real, allocator_type>        W(index_type(0), index_type(p-1), index_type(0), index_type(k-1));
    std::vector<real, allocator_type>              q(k);

    // We set a shape like:
    // y = [1 2 3 4]'
    // X = [1 1 0 0 ; 0 1 1 0 ; 0 0 1 1; 0 0 0 1]
    for (index_type idx_n = 0; idx_n < n; ++idx_n)
    {
        if (idx_n < p)
        {
            X(idx_n, idx_n) = 1;
            if (idx_n + 1 < p)
            {
                X(idx_n, idx_n+1) = 1;
            }
        }
        y[idx_n] = idx_n+1;
    }

    PartialLeastSquares<Pls2dArray, Pls1dArray>().pls_denham(X, y, n, p, k, &W, &q);

    std::vector<real, allocator_type>       t(p, 0);
    for (index_type i = 0; i < p; ++i)
    {
        for (index_type ik = 0; ik < k; ++ik)
        {
            t[i] += W(i, ik) * q[ik];
        }
    }

    std::vector<real, allocator_type>       Xt(n, 0);
    for (index_type i = 0; i < n; ++i)
    {
        for (index_type ip = 0; ip < p; ++ip)
        {
            Xt[i] += X(i, ip) * t[ip];
        }
        if (i != n-1)
        {
            // We know from the input that Xt and y must be close
            EXPECT_NEAR(Xt[i], y[i], 2e-6);
        }
        else
        {
            EXPECT_NEAR(Xt[i], 0, 1e-6);
        }
    }

    for (index_type i = 0; i < n-1; ++i)
    {
        // The objective of the function is to find such that
        // the difference between y and X t is low
        EXPECT_NEAR(std::abs(y[i] - Xt[i]), 0, 2e-6);
    }
}

TEST(PartialLeastSquaresTest, SimpleTestK){
    using allocator_type = gmx::SimdSetup<real>::allocator_type;
    using index_type     = Pls2dArray::index_type;

    const index_type                               n = 4;
    const index_type                               p = 3;
    const index_type                               k = 3;
    std::vector<real, allocator_type>              y(n, 0);
    gmx::IrregArray2D<real, allocator_type>        X(index_type(0), index_type(n-1), index_type(0), index_type(p-1));

    gmx::IrregArray2D<real, allocator_type>        W(index_type(0), index_type(p-1), index_type(0), index_type(k-1));
    std::vector<real, allocator_type>              q(k);

    // We set a shape like:
    // y = [1 2 3 4]'
    // X = [1 1 0 0 ; 0 1 1 0 ; 0 0 1 1; 0 0 0 1]
    for (index_type idx_n = 0; idx_n < n; ++idx_n)
    {
        if (idx_n < p)
        {
            X(idx_n, idx_n) = 1;
            if (idx_n + 1 < p)
            {
                X(idx_n, idx_n+1) = 1;
            }
        }
        y[idx_n] = idx_n+1;
    }

    PartialLeastSquares<Pls2dArray, Pls1dArray>().pls_denham(X, y, n, p, k, &W, &q);

    std::vector<real, allocator_type>       t(p, 0);
    for (index_type i = 0; i < p; ++i)
    {
        for (index_type ik = 0; ik < k; ++ik)
        {
            t[i] += W(i, ik) * q[ik];
        }
    }

    std::vector<real, allocator_type>       Xt(n, 0);
    for (index_type i = 0; i < n; ++i)
    {
        for (index_type ip = 0; ip < p; ++ip)
        {
            Xt[i] += X(i, ip) * t[ip];
        }
        if (i != n-1)
        {
            // We know from the input that Xt and y must be close
            EXPECT_NEAR(Xt[i], y[i], 1e-6);
        }
        else
        {
            EXPECT_NEAR(Xt[i], 0, 1e-6);
        }
    }

    for (index_type i = 0; i < n-1; ++i)
    {
        // The objective of the function is to find such that
        // the difference between y and X t is low
        EXPECT_NEAR(std::abs(y[i] - Xt[i]), 0, 1e-6);
    }
}


TEST(PartialLeastSquaresTest, LargerTest){
    using allocator_type = gmx::SimdSetup<real>::allocator_type;
    using index_type     = Pls2dArray::index_type;

    const index_type                               n = 40;
    const index_type                               p = 40;
    const index_type                               k = 40;
    std::vector<real, allocator_type>              y(n, 0);
    gmx::IrregArray2D<real, allocator_type>        X(index_type(0), index_type(n-1), index_type(0), index_type(p-1));

    gmx::IrregArray2D<real, allocator_type>        W(index_type(0), index_type(p-1), index_type(0), index_type(k-1));
    std::vector<real, allocator_type>              q(k);

    // We set a shape like:
    // y = [1 2 3 4]'
    // X = [1 1 0 0 ; 0 1 1 0 ; 0 0 1 1; 0 0 0 1]
    for (index_type idx_n = 0; idx_n < n; ++idx_n)
    {
        if (idx_n < p)
        {
            X(idx_n, idx_n) = 1;
            if (idx_n + 1 < p)
            {
                X(idx_n, idx_n+1) = 1;
            }
        }
        y[idx_n] = idx_n+1;
    }

    PartialLeastSquares<Pls2dArray, Pls1dArray>().pls_denham(X, y, n, p, k, &W, &q);

    std::vector<real, allocator_type>       t(p, 0);
    for (index_type i = 0; i < p; ++i)
    {
        for (index_type ik = 0; ik < k; ++ik)
        {
            t[i] += W(i, ik) * q[ik];
        }
    }

    std::vector<real, allocator_type>       Xt(n, 0);
    for (index_type i = 0; i < n; ++i)
    {
        for (index_type ip = 0; ip < p; ++ip)
        {
            Xt[i] += X(i, ip) * t[ip];
        }
        // We know from the input that Xt and y must be close
        EXPECT_NEAR(Xt[i], y[i], 4e-5);
    }

    for (index_type i = 0; i < n; ++i)
    {
        // The objective of the function is to find such that
        // the difference between y and X t is low
        EXPECT_NEAR(std::abs(y[i] - Xt[i]), 0, 4e-5);
    }
}

TEST(PartialLeastSquaresTest, ComparePrevious){
    using allocator_type = gmx::SimdSetup<real>::allocator_type;
    using index_type     = Pls2dArray::index_type;

    const index_type                               n = 10;
    const index_type                               p = 10;
    std::vector<real, allocator_type>              y(n, 0);
    gmx::IrregArray2D<real, allocator_type>        X(index_type(0), index_type(n-1), index_type(0), index_type(p-1));

    for (int idx_row = 0; idx_row < n; ++idx_row)
    {
        for (int idx_col = 0; idx_col < idx_row; ++idx_col)
        {
            X(idx_row, idx_col) = 0;
        }
        for (int idx_col = idx_row; idx_col < p; ++idx_col)
        {
            X(idx_row, idx_col) = 1;
        }
    }
    for (int idx_row = 0; idx_row < n; ++idx_row)
    {
        y[idx_row] = 1;
    }

    {
        const index_type                               k = 1;
        gmx::IrregArray2D<real, allocator_type>        W(index_type(0), index_type(p-1), index_type(0), index_type(k-1));
        std::vector<real, allocator_type>              q(k);

        PartialLeastSquares<Pls2dArray, Pls1dArray>().pls_denham(X, y, n, p, k, &W, &q);

        for (index_type i = 0; i < p; ++i)
        {
            EXPECT_TRUE(std::abs((W(i, 0)-(i+1))/(i+1)) < 1e-4);
        }


        EXPECT_TRUE(std::abs((q[0]-0.0226)/0.0226) < 1e-2);
    }

    {
        const index_type                               k = 9;
        gmx::IrregArray2D<real, allocator_type>        W(index_type(0), index_type(p-1), index_type(0), index_type(k-1));
        std::vector<real, allocator_type>              q(k);

        PartialLeastSquares<Pls2dArray, Pls1dArray>().pls_denham(X, y, n, p, k, &W, &q);

        const real W_test[n][k] = {{1.0000,   -0.2443,    0.1117,   -0.0553,    0.0243,   -0.0085,    0.0022,   -0.0004,    0.0000},
                                   {2.0000,   -0.4661,    0.1944,   -0.0829,    0.0288,   -0.0068,    0.0007,    0.0001,   -0.0000},
                                   {3.0000,   -0.6425,    0.2243,   -0.0670,    0.0086,    0.0036,   -0.0021,    0.0004,   -0.0000},
                                   {4.0000,   -0.7511,    0.1873,   -0.0118,   -0.0208,   0.0098,   -0.0012,   -0.0002,    0.0001},
                                   {5.0000,   -0.7692,    0.0846,    0.0556,   -0.0320,    0.0021,    0.0022,   -0.0003,   -0.0001},
                                   {6.0000,   -0.6742,   -0.0627,    0.0925,   -0.0086,   -0.0100,    0.0010,    0.0005,    0.0000},
                                   {7.0000,   -0.4434,   -0.2083,    0.0597,    0.0303,   -0.0049,   -0.0031,   -0.0004,   -0.0000},
                                   {8.0000,   -0.0543,   -0.2759,   -0.0435,    0.0288,    0.0135,    0.0022,    0.0001,    0.0000},
                                   {9.0000,    0.5158,   -0.1542,   -0.1261,   -0.0416,   -0.0073,   -0.0007,   -0.0000,   -0.0000},
                                   {10.0000,    1.2896,    0.3084,    0.0701,    0.0119,    0.0013,    0.0001,    0.0000,    0.0000}};

        for (int idx_row = 0; idx_row < n; ++idx_row)
        {
            // We stop at 5 because the accuracy really drops after
            for (int idx_col = 0; idx_col < 5; ++idx_col)
            {
                if (W_test[idx_row][idx_col] != 0)
                {
                    EXPECT_TRUE(std::abs((W(idx_row, idx_col)-W_test[idx_row][idx_col])/W_test[idx_row][idx_col]) < 1e-1);
                }
                else
                {
                    EXPECT_TRUE(std::abs(W(idx_row, idx_col)) < 1e-4);
                }
            }
        }

        const real q_test[k] = {
            0.0260,
            0.2912,
            0.8007,
            1.3350,
            1.7883,
            2.1876,
            2.5712,
            2.9524,
            3.3333
        };
        for (index_type i = 0; i < k; ++i)
        {
            EXPECT_TRUE(std::abs((q[i]-q_test[i])/q_test[i]) < 1e-2);
        }
    }
}


// The test that are assumed to throw should not be compiled with intel
// because it makes the compilation failed.
#if !defined(__INTEL_COMPILER) && !defined(__ICC) && !defined(__ICL)

TEST(PartialLeastSquaresTest, ParametersFailureTest){
    using allocator_type = gmx::SimdSetup<real>::allocator_type;
    using index_type     = Pls2dArray::index_type;

    PartialLeastSquares<Pls2dArray, Pls1dArray> pls;

    {
        const index_type                               n = 4;
        const index_type                               p = 4;
        const index_type                               k = 4;
        std::vector<real, allocator_type>              y(n+1, 0);
        gmx::IrregArray2D<real, allocator_type>        X(index_type(0), index_type(n-1), index_type(0), index_type(p-1));

        gmx::IrregArray2D<real, allocator_type>        W(index_type(0), index_type(p-1), index_type(0), index_type(k-1));
        std::vector<real, allocator_type>              q(k);

        EXPECT_ANY_THROW(pls.pls_denham(X, y, n, p, k, &W, &q));
    }
    {
        const index_type                               n = 4;
        const index_type                               p = 4;
        const index_type                               k = 4;
        std::vector<real, allocator_type>              y(n, 0);
        gmx::IrregArray2D<real, allocator_type>        X(index_type(0), index_type(n-1), index_type(0), index_type(p-1));

        gmx::IrregArray2D<real, allocator_type>        W(index_type(0), index_type(p-1), index_type(0), index_type(k-1));
        std::vector<real, allocator_type>              q(k+1);

        EXPECT_ANY_THROW(pls.pls_denham(X, y, n, p, k, &W, &q));
    }
    {
        const index_type                               n = 4;
        const index_type                               p = 4;
        const index_type                               k = 4+1;
        std::vector<real, allocator_type>              y(n, 0);
        gmx::IrregArray2D<real, allocator_type>        X(index_type(0), index_type(n-1), index_type(0), index_type(p-1));

        gmx::IrregArray2D<real, allocator_type>        W(index_type(0), index_type(p-1), index_type(0), index_type(k-1));
        std::vector<real, allocator_type>              q(k);

        EXPECT_ANY_THROW(pls.pls_denham(X, y, n, p, k, &W, &q));
    }
    {
        const index_type                               n = 4;
        const index_type                               p = 4;
        const index_type                               k = 4;
        std::vector<real, allocator_type>              y(n, 0);
        gmx::IrregArray2D<real, allocator_type>        X(index_type(0), index_type(n), index_type(0), index_type(p-1));

        gmx::IrregArray2D<real, allocator_type>        W(index_type(0), index_type(p-1), index_type(0), index_type(k-1));
        std::vector<real, allocator_type>              q(k);

        EXPECT_ANY_THROW(pls.pls_denham(X, y, n, p, k, &W, &q));
    }
    {
        const index_type                               n = 4;
        const index_type                               p = 4;
        const index_type                               k = 4;
        std::vector<real, allocator_type>              y(n, 0);
        gmx::IrregArray2D<real, allocator_type>        X(index_type(0), index_type(n-1), index_type(0), index_type(p-1));

        gmx::IrregArray2D<real, allocator_type>        W(index_type(0), index_type(p-1), index_type(0), index_type(k-1));
        std::vector<real, allocator_type>              q(k+2);

        EXPECT_ANY_THROW(pls.pls_denham(X, y, n, p, k, &W, &q));
    }
}

#endif

} // namespace

} // namespace partial_least_squares_tests

} // namespace gmx
