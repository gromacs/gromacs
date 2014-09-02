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

#include <vector>

#include <gtest/gtest.h>

namespace gmx
{

namespace partial_least_squares_tests
{

namespace
{

TEST(PartialLeastSquares, SimpleTest){
    const size_t            n = 10;
    const int               k = 10;
    std::vector<real>       y(n, 1);
    gmx::IrregArray2D<real> X(int(0), n-1, int(0), k-1);
    gmx::IrregArray2D<real> W(int(0), k-1, int(0), k-1);
    std::vector<real>       q(k);

    for (int idx_n = 0; idx_n < n; ++idx_n)
    {
        for (int idx_k = 0; idx_k < n; ++idx_k)
        {
            X(idx_n, idx_k) = idx_n*idx_k;
        }
    }

    const int               retValue = pls_denham(X, y, n, k, &W, &q);
    EXPECT_EQ(retValue, 0);
}

TEST(PartialLeastSquares, ParametersFailureTest){
}

} // namespace

} // namespace partial_least_squares_tests

} // namespace gmx
