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
            EXPECT_TRUE(std::abs((W(i, 0)-real(i+1))/real(i+1)) < 1e-4);
        }


        EXPECT_TRUE(std::abs((q[0]-real(0.0226))/real(0.0226)) < 1e-2);
    }

    {
        const index_type                               k = 10;
        gmx::IrregArray2D<real, allocator_type>        W(index_type(0), index_type(p-1), index_type(0), index_type(k-1));
        std::vector<real, allocator_type>              q(k);

        PartialLeastSquares<Pls2dArray, Pls1dArray>().pls_denham(X, y, n, p, k, &W, &q);

        const double W_test[n][k] = {{1.000000e+00, -2.443439e-01, 1.117034e-01, -5.526516e-02, 2.432813e-02, -8.542338e-03, 2.240709e-03, -4.115256e-04, 4.726970e-05, -2.558759e-06},
                                     {2.000000e+00, -4.660633e-01, 1.944467e-01, -8.289773e-02, 2.883333e-02, -6.802232e-03, 7.469030e-04, 8.382929e-05, -3.851605e-05, 3.838138e-06},
                                     {3.000000e+00, -6.425339e-01, 2.242820e-01, -6.695586e-02, 8.559896e-03, 3.638403e-03, -2.111437e-03, 4.047841e-04, -1.427518e-05, -3.542897e-06},
                                     {4.000000e+00, -7.511312e-01, 1.872862e-01, -1.177746e-02, -2.081590e-02, 9.817555e-03, -1.223537e-03, -2.459188e-04, 5.602335e-05, 2.361931e-06},
                                     {5.000000e+00, -7.692308e-01, 8.457316e-02, 5.563388e-02, -3.198698e-02, 2.146884e-03, 2.217552e-03, -2.986785e-04, -5.790875e-05, -1.180966e-06},
                                     {6.000000e+00, -6.742081e-01, -6.269393e-02, 9.252793e-02, -8.569091e-03, -1.002498e-02, 1.043847e-03, 5.380023e-04, 3.541861e-05, 4.428621e-07},
                                     {7.000000e+00, -4.434389e-01, -2.082902e-01, 5.973322e-02, 3.032281e-02, -4.884565e-03, -3.130808e-03, -3.764990e-04, -1.416824e-05, -1.215700e-07},
                                     {8.000000e+00, -5.429864e-02, -2.759169e-01, -4.350938e-02, 2.875978e-02, 1.352860e-02, 2.164201e-03, 1.452357e-04, 3.675740e-06, 2.315619e-08},
                                     {9.000000e+00, 5.158371e-01, -1.541889e-01, -1.261035e-01, -4.157664e-02, -7.299406e-03, -6.766977e-04, -3.062993e-05, -5.664125e-07, -2.742181e-09},
                                     {1.000000e+01, 1.289593e+00, 3.083778e-01, 7.005748e-02, 1.187904e-02, 1.303465e-03, 8.354292e-05, 2.784539e-06, 3.960927e-08, 1.523425e-10}};

        for (int idx_row = 0; idx_row < n; ++idx_row)
        {
            // We stop at 6 because the accuracy really drops after
            for (int idx_col = 0; idx_col < 6; ++idx_col)
            {
                if (W_test[idx_row][idx_col] != 0)
                {
                    EXPECT_TRUE(std::abs((W(idx_row, idx_col)-real(W_test[idx_row][idx_col]))/real(W_test[idx_row][idx_col])) < 1e-3);
                }
                else
                {
                    EXPECT_TRUE(std::abs(W(idx_row, idx_col)) < 1e-4);
                }
            }
        }

        const double q_test[k] = {
            2.597403e-02, 2.911726e-01, 8.007007e-01,
            1.335042e+00, 1.788281e+00, 2.187625e+00,
            2.571214e+00, 2.952373e+00, 3.333333e+00,
            3.714286e+00
        };
        for (index_type i = 0; i < k-1; ++i)
        {
            EXPECT_TRUE(std::abs((q[i]-real(q_test[i]))/real(q_test[i])) < 1e-2);
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
