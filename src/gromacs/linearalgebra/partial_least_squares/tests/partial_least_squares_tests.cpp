/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017,2018, by the GROMACS development team, led by
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

namespace gmx
{

namespace partial_least_squares_tests
{

namespace
{

TEST(PartialLeastSquares, SimpleTest){
    using allocator_type = gmx::SimdSetup<real>::allocator_type;
    using index_type     = Pls2dArray::index_type;

    const int                               n = 4;
    const int                               p = 4;
    const int                               k = 4;
    std::vector<real, allocator_type>       y(n, 0);
    gmx::IrregArray2D<real, allocator_type> X(index_type(0), index_type(n-1), index_type(0), index_type(p-1));

    gmx::IrregArray2D<real, allocator_type> W(index_type(0), index_type(p-1), index_type(0), index_type(k-1));
    std::vector<real, allocator_type>       q(k);

    // We set a shape like:
    // y = [1 2 3 4]'
    // X = [1 1 0 0 ; 0 1 1 0 ; 0 0 1 1; 0 0 0 1]
    for (int idx_n = 0; idx_n < n; ++idx_n)
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

    const int               retValue = pls_denham(X, y, n, p, k, &W, &q);
    EXPECT_EQ(retValue, 0);

    std::vector<real, allocator_type>       t(p, 0);
    for (int i = 0; i < p; ++i)
    {
        for (int ik = 0; ik < k; ++ik)
        {
            t[i] += W(i)[ik] * q[ik];
        }
    }

    std::vector<real, allocator_type>       Xt(n, 0);
    for (int i = 0; i < n; ++i)
    {
        for (int ip = 0; ip < p; ++ip)
        {
            Xt[i] += X(i)[ip] * t[ip];
        }
        // We know from the input that Xt and y must be close
        EXPECT_NEAR(Xt[i], y[i], 2e-6);
    }

    for (int i = 0; i < n; ++i)
    {
        // The objective of the function is to find such that
        // the difference between y and X t is low
        EXPECT_NEAR(std::abs(y[i] - Xt[i]), 0, 2e-6);
    }
}

TEST(PartialLeastSquares, SimpleTestN){
    using allocator_type = gmx::SimdSetup<real>::allocator_type;
    using index_type     = Pls2dArray::index_type;

    const int                               n = 5;
    const int                               p = 4;
    const int                               k = 4;
    std::vector<real, allocator_type>       y(n, 0);
    gmx::IrregArray2D<real, allocator_type> X(index_type(0), index_type(n-1), index_type(0), index_type(p-1));

    gmx::IrregArray2D<real, allocator_type> W(index_type(0), index_type(p-1), index_type(0), index_type(k-1));
    std::vector<real, allocator_type>       q(k);

    // We set a shape like:
    // y = [1 2 3 4]'
    // X = [1 1 0 0 ; 0 1 1 0 ; 0 0 1 1; 0 0 0 1]
    for (int idx_n = 0; idx_n < n; ++idx_n)
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

    const int               retValue = pls_denham(X, y, n, p, k, &W, &q);
    EXPECT_EQ(retValue, 0);

    std::vector<real, allocator_type>       t(p, 0);
    for (int i = 0; i < p; ++i)
    {
        for (int ik = 0; ik < k; ++ik)
        {
            t[i] += W(i)[ik] * q[ik];
        }
    }

    std::vector<real, allocator_type>       Xt(n, 0);
    for (int i = 0; i < n; ++i)
    {
        for (int ip = 0; ip < p; ++ip)
        {
            Xt[i] += X(i)[ip] * t[ip];
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

    for (int i = 0; i < n-1; ++i)
    {
        // The objective of the function is to find such that
        // the difference between y and X t is low
        EXPECT_NEAR(std::abs(y[i] - Xt[i]), 0, 2e-6);
    }
}

TEST(PartialLeastSquares, SimpleTestK){
    using allocator_type = gmx::SimdSetup<real>::allocator_type;
    using index_type     = Pls2dArray::index_type;

    const int                               n = 4;
    const int                               p = 3;
    const int                               k = 3;
    std::vector<real, allocator_type>       y(n, 0);
    gmx::IrregArray2D<real, allocator_type> X(index_type(0), index_type(n-1), index_type(0), index_type(p-1));

    gmx::IrregArray2D<real, allocator_type> W(index_type(0), index_type(p-1), index_type(0), index_type(k-1));
    std::vector<real, allocator_type>       q(k);

    // We set a shape like:
    // y = [1 2 3 4]'
    // X = [1 1 0 0 ; 0 1 1 0 ; 0 0 1 1; 0 0 0 1]
    for (int idx_n = 0; idx_n < n; ++idx_n)
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

    const int               retValue = pls_denham(X, y, n, p, k, &W, &q);
    EXPECT_EQ(retValue, 0);

    std::vector<real, allocator_type>       t(p, 0);
    for (int i = 0; i < p; ++i)
    {
        for (int ik = 0; ik < k; ++ik)
        {
            t[i] += W(i)[ik] * q[ik];
        }
    }

    std::vector<real, allocator_type>       Xt(n, 0);
    for (int i = 0; i < n; ++i)
    {
        for (int ip = 0; ip < p; ++ip)
        {
            Xt[i] += X(i)[ip] * t[ip];
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

    for (int i = 0; i < n-1; ++i)
    {
        // The objective of the function is to find such that
        // the difference between y and X t is low
        EXPECT_NEAR(std::abs(y[i] - Xt[i]), 0, 1e-6);
    }
}


TEST(PartialLeastSquares, LargerTest){
    using allocator_type = gmx::SimdSetup<real>::allocator_type;
    using index_type     = Pls2dArray::index_type;

    const int                               n = 40;
    const int                               p = 40;
    const int                               k = 40;
    std::vector<real, allocator_type>       y(n, 0);
    gmx::IrregArray2D<real, allocator_type> X(index_type(0), index_type(n-1), index_type(0), index_type(p-1));

    gmx::IrregArray2D<real, allocator_type> W(index_type(0), index_type(p-1), index_type(0), index_type(k-1));
    std::vector<real, allocator_type>       q(k);

    // We set a shape like:
    // y = [1 2 3 4]'
    // X = [1 1 0 0 ; 0 1 1 0 ; 0 0 1 1; 0 0 0 1]
    for (int idx_n = 0; idx_n < n; ++idx_n)
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

    const int               retValue = pls_denham(X, y, n, p, k, &W, &q);
    EXPECT_EQ(retValue, 0);

    std::vector<real, allocator_type>       t(p, 0);
    for (int i = 0; i < p; ++i)
    {
        for (int ik = 0; ik < k; ++ik)
        {
            t[i] += W(i)[ik] * q[ik];
        }
    }

    std::vector<real, allocator_type>       Xt(n, 0);
    for (int i = 0; i < n; ++i)
    {
        for (int ip = 0; ip < p; ++ip)
        {
            Xt[i] += X(i)[ip] * t[ip];
        }
        // We know from the input that Xt and y must be close
        EXPECT_NEAR(Xt[i], y[i], 4e-5);
    }

    for (int i = 0; i < n; ++i)
    {
        // The objective of the function is to find such that
        // the difference between y and X t is low
        EXPECT_NEAR(std::abs(y[i] - Xt[i]), 0, 4e-5);
    }
}


// The test that are assumed to throw should not be compiled with intel
// because it makes the compilation failed.
#if !defined(__INTEL_COMPILER) && !defined(__ICC) && !defined(__ICL)

TEST(PartialLeastSquares, ParametersFailureTest){
    using allocator_type = gmx::SimdSetup<real>::allocator_type;
    using index_type     = Pls2dArray::index_type;

    {
        const int                               n = 4;
        const int                               p = 4;
        const int                               k = 4;
        std::vector<real, allocator_type>       y(n+1, 0);
        gmx::IrregArray2D<real, allocator_type> X(index_type(0), index_type(n-1), index_type(0), index_type(p-1));

        gmx::IrregArray2D<real, allocator_type> W(index_type(0), index_type(p-1), index_type(0), index_type(k-1));
        std::vector<real, allocator_type>       q(k);

        EXPECT_ANY_THROW(pls_denham(X, y, n, p, k, &W, &q));
    }
    {
        const int                               n = 4;
        const int                               p = 4;
        const int                               k = 4;
        std::vector<real, allocator_type>       y(n, 0);
        gmx::IrregArray2D<real, allocator_type> X(index_type(0), index_type(n-1), index_type(0), index_type(p-1));

        gmx::IrregArray2D<real, allocator_type> W(index_type(0), index_type(p-1), index_type(0), index_type(k-1));
        std::vector<real, allocator_type>       q(k);

        EXPECT_ANY_THROW(pls_denham(X, y, n, p, k, &W, &q));
    }
    {
        const int                               n = 4;
        const int                               p = 4;
        const int                               k = 4+1;
        std::vector<real, allocator_type>       y(n, 0);
        gmx::IrregArray2D<real, allocator_type> X(index_type(0), index_type(n-1), index_type(0), index_type(p-1));

        gmx::IrregArray2D<real, allocator_type> W(index_type(0), index_type(p-1), index_type(0), index_type(k-1));
        std::vector<real, allocator_type>       q(k);

        EXPECT_ANY_THROW(pls_denham(X, y, n, p, k, &W, &q));
    }
    {
        const int                               n = 4;
        const int                               p = 4;
        const int                               k = 4;
        std::vector<real, allocator_type>       y(n, 0);
        gmx::IrregArray2D<real, allocator_type> X(index_type(0), index_type(n-1), index_type(0), index_type(p-1));

        gmx::IrregArray2D<real, allocator_type> W(index_type(0), index_type(p-1), index_type(0), index_type(k-1));
        std::vector<real, allocator_type>       q(k);

        EXPECT_ANY_THROW(pls_denham(X, y, n, p, k, &W, &q));
    }
    {
        const int                               n = 4;
        const int                               p = 4;
        const int                               k = 4;
        std::vector<real, allocator_type>       y(n, 0);
        gmx::IrregArray2D<real, allocator_type> X(index_type(0), index_type(n-1), index_type(0), index_type(p-1));

        gmx::IrregArray2D<real, allocator_type> W(index_type(0), index_type(p-1), index_type(0), index_type(k-1));
        std::vector<real, allocator_type>       q(k+2);

        EXPECT_ANY_THROW(pls_denham(X, y, n, p, k, &W, &q));
    }
}

#endif

} // namespace

} // namespace partial_least_squares_tests

} // namespace gmx
