/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2015- The GROMACS Authors
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
 * Tests matrix inversion routines
 *
 * \todo Test error conditions when they throw exceptions
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_math
 */
#include "gmxpre.h"

#include "gromacs/math/invertmatrix.h"

#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/real.h"

#include "testutils/testasserts.h"

namespace gmx
{
namespace test
{
namespace
{

using gmx::invertMatrix;
using gmx::test::defaultRealTolerance;

TEST(InvertMatrixTest, IdentityIsImpotent)
{
    matrix in = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
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
    matrix in = { { 1, 2, 3 }, { -1, real(2.5), 1 }, { 10, -2, real(1.2) } };
    matrix out;
    matrix expected = {
        { real(-0.12019230769230768), real(0.20192307692307693), real(0.13221153846153844) },
        { real(-0.26923076923076916), real(0.69230769230769229), real(0.096153846153846145) },
        { real(0.55288461538461531), real(-0.52884615384615374), real(-0.10817307692307691) }
    };

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

} // namespace
} // namespace test
} // namespace gmx
