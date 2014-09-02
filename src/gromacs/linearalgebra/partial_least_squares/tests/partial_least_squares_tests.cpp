/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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
    const int               n = 4;
    const int               p = 4;
    const int               k = 3;
    std::vector<real>       y(n, 0);
    gmx::IrregArray2D<real> X(int(0), n-1, int(0), p-1);

    gmx::IrregArray2D<real> W(int(0), p-1, int(0), k-1);
    std::vector<real>       q(k);

    for (int idx_n = 0; idx_n < n; ++idx_n)
    {
        X(idx_n, idx_n) = 1;
        if (idx_n + 1 != p)
        {
            X(idx_n, idx_n+1) = 1;
        }
        y[idx_n] = idx_n+1;
    }

    const int               retValue = pls_denham(X, y, n, p, k, &W, &q);
    EXPECT_EQ(retValue, 0);


    // WIP
    for (int i = 0; i < k; ++i)
    {
        std::cout <<"q[" << i << "]=" << q[i] << std::endl;
    }
    for (int i = 0; i < p; ++i)
    {
        for (int ik = 0; ik < k; ++ik)
        {
            std::cout <<"W[" << i << "][" << ik << "]=" << W(i)[ik] << std::endl;
        }
    }

    std::vector<real>       t(p, 0);
    for (int i = 0; i < p; ++i)
    {
        for (int ik = 0; ik < k; ++ik)
        {
            t[i] += W(i)[ik] * q[ik];
        }
        std::cout <<"t[" << i << "]=" << t[i] << std::endl;
    }


    std::vector<real>       Xt(n, 0);
    for (int i = 0; i < n; ++i)
    {
        for (int ip = 0; ip < p; ++ip)
        {
            Xt[i] += X(i)[ip] * t[ip];
        }
        std::cout <<"Xt[" << i << "]=" << Xt[i] << std::endl;
    }

    for (int i = 0; i < n; ++i)
    {
        std::cout <<"abs(y-X*t)[" << i << "]=" << std::abs(y[i] - Xt[i]) << std::endl;
    }
}

TEST(PartialLeastSquares, ParametersFailureTest){
}

} // namespace

} // namespace partial_least_squares_tests

} // namespace gmx
