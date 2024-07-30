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
 * \brief Tests for GROMACS gamma distribution
 *
 * \author Erik Lindahl <erik.lindahl@gmail.com>
 * \ingroup module_random
 */
#include "gmxpre.h"

#include "gromacs/random/gammadistribution.h"

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

TEST(GammaDistributionTest, Output)
{
    gmx::test::TestReferenceData    data;
    gmx::test::TestReferenceChecker checker(data.rootChecker());

    gmx::ThreeFry2x64<8>         rng(123456, gmx::RandomDomain::Other);
    gmx::GammaDistribution<real> distA(0.1, 5.0);
    gmx::GammaDistribution<real> distB(1.0, 5.0);
    gmx::GammaDistribution<real> distC(2.0, 5.0);
    std::vector<real>            result;

    result.reserve(10);
    for (int i = 0; i < 10; i++)
    {
        result.push_back(distA(rng));
    }
    for (int i = 0; i < 10; i++)
    {
        result.push_back(distB(rng));
    }
    for (int i = 0; i < 10; i++)
    {
        result.push_back(distC(rng));
    }

    // Use wider tolerance in double precision
    checker.setDefaultTolerance(test::relativeToleranceAsPrecisionDependentUlp(1.0, 4, 8));
    checker.checkSequence(result.begin(), result.end(), "GammaDistribution");
}

TEST(GammaDistributionTest, Logical)
{
    gmx::ThreeFry2x64<8>         rng(123456, gmx::RandomDomain::Other);
    gmx::GammaDistribution<real> distA(2.0, 5.0);
    gmx::GammaDistribution<real> distB(2.0, 5.0);
    gmx::GammaDistribution<real> distC(3.0, 5.0);
    gmx::GammaDistribution<real> distD(2.0, 4.0);

    EXPECT_EQ(distA, distB);
    EXPECT_NE(distA, distC);
    EXPECT_NE(distA, distD);
}


TEST(GammaDistributionTest, Reset)
{
    gmx::ThreeFry2x64<8>                  rng(123456, gmx::RandomDomain::Other);
    gmx::GammaDistribution<real>          distA(2.0, 5.0);
    gmx::GammaDistribution<real>          distB(2.0, 5.0);
    gmx::GammaDistribution<>::result_type valA, valB;

    valA = distA(rng);

    distB(rng);
    rng.restart();
    distB.reset();

    valB = distB(rng);

    EXPECT_REAL_EQ_TOL(valA, valB, gmx::test::ulpTolerance(0));
}

TEST(GammaDistributionTest, AltParam)
{
    gmx::ThreeFry2x64<8>                     rngA(123456, gmx::RandomDomain::Other);
    gmx::ThreeFry2x64<8>                     rngB(123456, gmx::RandomDomain::Other);
    gmx::GammaDistribution<real>             distA(2.0, 5.0);
    gmx::GammaDistribution<real>             distB; // default parameters
    gmx::GammaDistribution<real>::param_type paramA(2.0, 5.0);

    EXPECT_NE(distA(rngA), distB(rngB));
    rngA.restart();
    rngB.restart();
    distA.reset();
    distB.reset();
    EXPECT_REAL_EQ_TOL(distA(rngA), distB(rngB, paramA), gmx::test::ulpTolerance(0));
}

} // namespace
} // namespace test
} // namespace gmx
