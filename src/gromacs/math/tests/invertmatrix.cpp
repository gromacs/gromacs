/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015,2016, by the GROMACS development team, led by
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
 * Tests matrix inversion routines
 *
 * \todo Test error conditions when they throw exceptions
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_math
 */
#include "gmxpre.h"

#include "gromacs/math/invertmatrix.h"

#include <vector>

#include <gtest/gtest.h>

#include "gromacs/math/vec.h"

#include "testutils/testasserts.h"

namespace
{

using gmx::invertMatrix;
using gmx::invertBoxMatrix;
using gmx::test::defaultRealTolerance;

TEST(InvertMatrixTest, IdentityIsImpotent)
{
    matrix in = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    matrix out;

    invertMatrix(in, out);

    EXPECT_REAL_EQ_TOL(out[XX][XX], in[XX][XX], defaultRealTolerance());
    EXPECT_REAL_EQ_TOL(out[XX][YY], in[XX][YY], defaultRealTolerance());
    EXPECT_REAL_EQ_TOL(out[XX][ZZ], in[XX][ZZ], defaultRealTolerance());
    EXPECT_REAL_EQ_TOL(out[YY][XX], in[YY][XX], defaultRealTolerance());
    EXPECT_REAL_EQ_TOL(out[YY][YY], in[YY][YY], defaultRealTolerance());
    EXPECT_REAL_EQ_TOL(out[YY][ZZ], in[YY][ZZ], defaultRealTolerance());
    EXPECT_REAL_EQ_TOL(out[ZZ][XX], in[ZZ][XX], defaultRealTolerance());
    EXPECT_REAL_EQ_TOL(out[ZZ][YY], in[ZZ][YY], defaultRealTolerance());
    EXPECT_REAL_EQ_TOL(out[ZZ][ZZ], in[ZZ][ZZ], defaultRealTolerance());
}

TEST(InvertMatrixTest, ComputesInverse)
{
    matrix in = {{1, 2, 3}, {-1, real(2.5), 1}, {10, -2, real(1.2)}};
    matrix out;
    matrix expected = {{real(-0.12019230769230768),
                        real(0.20192307692307693),
                        real(0.13221153846153844)},
                       {real(-0.26923076923076916),
                        real(0.69230769230769229),
                        real(0.096153846153846145)},
                       {real(0.55288461538461531),
                        real(-0.52884615384615374),
                        real(-0.10817307692307691)}};

    invertMatrix(in, out);

    EXPECT_REAL_EQ_TOL(out[XX][XX], expected[XX][XX], defaultRealTolerance());
    EXPECT_REAL_EQ_TOL(out[XX][YY], expected[XX][YY], defaultRealTolerance());
    EXPECT_REAL_EQ_TOL(out[XX][ZZ], expected[XX][ZZ], defaultRealTolerance());
    EXPECT_REAL_EQ_TOL(out[YY][XX], expected[YY][XX], defaultRealTolerance());
    EXPECT_REAL_EQ_TOL(out[YY][YY], expected[YY][YY], defaultRealTolerance());
    EXPECT_REAL_EQ_TOL(out[YY][ZZ], expected[YY][ZZ], defaultRealTolerance());
    EXPECT_REAL_EQ_TOL(out[ZZ][XX], expected[ZZ][XX], defaultRealTolerance());
    EXPECT_REAL_EQ_TOL(out[ZZ][YY], expected[ZZ][YY], defaultRealTolerance());
    EXPECT_REAL_EQ_TOL(out[ZZ][ZZ], expected[ZZ][ZZ], defaultRealTolerance());
}

TEST(InvertBoxMatrixTest, IdentityIsImpotent)
{
    matrix in = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};

    invertBoxMatrix(in, in);

    EXPECT_REAL_EQ_TOL(in[XX][XX], in[XX][XX], defaultRealTolerance());
    EXPECT_REAL_EQ_TOL(in[XX][YY], in[XX][YY], defaultRealTolerance());
    EXPECT_REAL_EQ_TOL(in[XX][ZZ], in[XX][ZZ], defaultRealTolerance());
    EXPECT_REAL_EQ_TOL(in[YY][XX], in[YY][XX], defaultRealTolerance());
    EXPECT_REAL_EQ_TOL(in[YY][YY], in[YY][YY], defaultRealTolerance());
    EXPECT_REAL_EQ_TOL(in[YY][ZZ], in[YY][ZZ], defaultRealTolerance());
    EXPECT_REAL_EQ_TOL(in[ZZ][XX], in[ZZ][XX], defaultRealTolerance());
    EXPECT_REAL_EQ_TOL(in[ZZ][YY], in[ZZ][YY], defaultRealTolerance());
    EXPECT_REAL_EQ_TOL(in[ZZ][ZZ], in[ZZ][ZZ], defaultRealTolerance());
}

TEST(InvertBoxMatrixTest, ComputesInverseInPlace)
{
    matrix in       = {{1, 0, 0}, {-1, real(2.5), 0}, {10, -2, real(1.2)}};
    matrix expected = {{1, 0, 0},
                       {real(0.4), real(0.4), 0},
                       {real(-23.0/3.0), real(2.0/3.0), real(5.0/6.0)}};

    invertBoxMatrix(in, in);

    EXPECT_REAL_EQ_TOL(expected[XX][XX], in[XX][XX], defaultRealTolerance());
    EXPECT_REAL_EQ_TOL(expected[XX][YY], in[XX][YY], defaultRealTolerance());
    EXPECT_REAL_EQ_TOL(expected[XX][ZZ], in[XX][ZZ], defaultRealTolerance());
    EXPECT_REAL_EQ_TOL(expected[YY][XX], in[YY][XX], defaultRealTolerance());
    EXPECT_REAL_EQ_TOL(expected[YY][YY], in[YY][YY], defaultRealTolerance());
    EXPECT_REAL_EQ_TOL(expected[YY][ZZ], in[YY][ZZ], defaultRealTolerance());
    EXPECT_REAL_EQ_TOL(expected[ZZ][XX], in[ZZ][XX], defaultRealTolerance());
    EXPECT_REAL_EQ_TOL(expected[ZZ][YY], in[ZZ][YY], defaultRealTolerance());
    EXPECT_REAL_EQ_TOL(expected[ZZ][ZZ], in[ZZ][ZZ], defaultRealTolerance());
}

} // namespace
