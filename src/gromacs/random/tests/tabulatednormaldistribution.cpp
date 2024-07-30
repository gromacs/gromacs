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
 * \brief Tests for GROMACS tabulated normal distribution
 *
 * \author Erik Lindahl <erik.lindahl@gmail.com>
 * \ingroup module_random
 */
#include "gmxpre.h"

#include "gromacs/random/tabulatednormaldistribution.h"

#include <cmath>
#include <cstddef>

#include <array>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/random/seed.h"
#include "gromacs/random/threefry.h"
#include "gromacs/utility/real.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"

namespace gmx
{
namespace test
{
namespace
{

TEST(TabulatedNormalDistributionTest, Output14)
{
    gmx::test::TestReferenceData    data;
    gmx::test::TestReferenceChecker checker(data.rootChecker());

    gmx::ThreeFry2x64<2>               rng(123456, gmx::RandomDomain::Other);
    gmx::TabulatedNormalDistribution<> dist(2.0, 5.0); // Use default 14-bit resolution
    std::vector<float>                 result;

    result.reserve(10);
    for (int i = 0; i < 10; i++)
    {
        result.push_back(dist(rng));
    }
    checker.checkSequence(result.begin(), result.end(), "TabulatedNormalDistribution14");
}

TEST(TabulatedNormalDistributionTest, Output16)
{
    gmx::test::TestReferenceData    data;
    gmx::test::TestReferenceChecker checker(data.rootChecker());

    gmx::ThreeFry2x64<2>                        rng(123456, gmx::RandomDomain::Other);
    gmx::TabulatedNormalDistribution<float, 16> dist(2.0, 5.0); // Use larger 16-bit table
    std::vector<float>                          result;

    result.reserve(10);
    for (int i = 0; i < 10; i++)
    {
        result.push_back(dist(rng));
    }
    checker.checkSequence(result.begin(), result.end(), "TabulatedNormalDistribution16");
}

TEST(TabulatedNormalDistributionTest, OutputDouble14)
{
    gmx::test::TestReferenceData    data;
    gmx::test::TestReferenceChecker checker(data.rootChecker());

    gmx::ThreeFry2x64<2>                     rng(123456, gmx::RandomDomain::Other);
    gmx::TabulatedNormalDistribution<double> dist(2.0, 5.0);
    std::vector<double>                      result;

    result.reserve(10);
    for (int i = 0; i < 10; i++)
    {
        result.push_back(dist(rng));
    }
    checker.checkSequence(result.begin(), result.end(), "TabulatedNormalDistributionDouble14");
}

TEST(TabulatedNormalDistributionTest, Logical)
{
    gmx::ThreeFry2x64<2>               rng(123456, gmx::RandomDomain::Other);
    gmx::TabulatedNormalDistribution<> distA(2.0, 5.0);
    gmx::TabulatedNormalDistribution<> distB(2.0, 5.0);
    gmx::TabulatedNormalDistribution<> distC(3.0, 5.0);
    gmx::TabulatedNormalDistribution<> distD(2.0, 4.0);

    EXPECT_EQ(distA, distB);
    EXPECT_NE(distA, distC);
    EXPECT_NE(distA, distD);
}


TEST(TabulatedNormalDistributionTest, Reset)
{
    gmx::ThreeFry2x64<2>                            rng(123456, gmx::RandomDomain::Other);
    gmx::TabulatedNormalDistribution<>              distA(2.0, 5.0);
    gmx::TabulatedNormalDistribution<>              distB(2.0, 5.0);
    gmx::TabulatedNormalDistribution<>::result_type valA, valB;

    valA = distA(rng);

    distB(rng);
    rng.restart();
    distB.reset();

    valB = distB(rng);

    EXPECT_REAL_EQ_TOL(valA, valB, gmx::test::ulpTolerance(0));
}

TEST(TabulatedNormalDistributionTest, AltParam)
{
    gmx::ThreeFry2x64<2>                           rngA(123456, gmx::RandomDomain::Other);
    gmx::ThreeFry2x64<2>                           rngB(123456, gmx::RandomDomain::Other);
    gmx::TabulatedNormalDistribution<>             distA(2.0, 5.0);
    gmx::TabulatedNormalDistribution<>             distB;
    gmx::TabulatedNormalDistribution<>::param_type paramA(2.0, 5.0);

    EXPECT_NE(distA(rngA), distB(rngB));
    rngA.restart();
    rngB.restart();
    distA.reset();
    distB.reset();
    EXPECT_REAL_EQ_TOL(distA(rngA), distB(rngB, paramA), gmx::test::ulpTolerance(0));
}

TEST(TabulatedNormalDistributionTableTest, HasValidProperties)
{
    auto table = TabulatedNormalDistribution<real>::makeTable();

    EXPECT_EQ(table.size() % 2, 0) << "Table must have even number of entries";

    size_t halfSize     = table.size() / 2;
    double sumOfSquares = 0.0;
    // accept errors of a few ULP since the exact value of the summation
    // below will depend on whether the compiler issues FMA instructions
    const auto elementTolerance = gmx::test::ulpTolerance(10);
    for (size_t i = 0, iFromEnd = table.size() - 1; i < halfSize; ++i, --iFromEnd)
    {
        EXPECT_REAL_EQ_TOL(table.at(i), -table.at(iFromEnd), elementTolerance)
                << "Table is not an odd-valued function for entries " << i << " and " << iFromEnd;
        // Add up the squares of the table values in order of ascending
        // magnitude (to minimize accumulation of round-off error).
        sumOfSquares += table.at(i) * table.at(i) + table.at(iFromEnd) * table.at(iFromEnd);
    }

    /* We calculate the sum of N = table.size() positive values in ascending order.
     * On average, we have N / 2 single-bit differences. Up- and down errors will cancel,
     * reducing the error by a factor on sqrt(N), leading to the tolerance of sqrt(N)/2.
     * See analysis by Erik Lindahl in #4700. */
    const auto   varianceTolerance = gmx::test::ulpTolerance(std::sqrt(table.size()) / 2);
    const double variance          = sumOfSquares / table.size();
    EXPECT_REAL_EQ_TOL(1.0, variance, varianceTolerance) << "Table should have unit variance";
}

} // namespace
} // namespace test
} // namespace gmx
