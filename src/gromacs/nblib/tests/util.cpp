/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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
 * This implements basic nblib utility tests
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */
#include "gmxpre.h"

#include "gromacs/nblib/util.h"

#include <vector>

#include "gromacs/math/vectypes.h"

#include "testutils/testasserts.h"

namespace nblib
{
namespace test
{
namespace
{

TEST(NBlibTest, checkNumericValues)
{
    std::vector<gmx::RVec> vec;
    vec.push_back({1., 1., 1.});
    vec.push_back({ 2., 2., 2. });

    bool ret = checkNumericValues(vec);
    EXPECT_EQ(ret, true);
}

TEST(NBlibTest, checkNumericValuesHasNan)
{
    std::vector<gmx::RVec> vec;
    vec.push_back({1., 1., 1.});
    vec.push_back({ 2., 2., 2. });

    vec.push_back({ (real)NAN, (real)NAN, (real)NAN });

    bool ret = checkNumericValues(vec);
    EXPECT_EQ(ret, false);
}

TEST(NBlibTest, checkNumericValuesHasInf)
{
    std::vector<gmx::RVec> vec;
    vec.push_back({1., 1., 1.});
    vec.push_back({ 2., 2., 2. });

    vec.push_back({ (real)INFINITY, (real)INFINITY, (real)INFINITY });

    bool ret = checkNumericValues(vec);
    EXPECT_EQ(ret, false);
}


TEST(NBlibTest, generateVelocity)
{
    constexpr size_t N = 10;
    std::vector<real> masses(N, 1.0);
    auto out = generateVelocity(300.0, 1, masses);

    std::vector<gmx::RVec> expected = {
        { 2.698104, 0.553971 ,-0.669996 },
        { 1.208262, -0.947503 ,-0.393945 },
        { -0.256397, 0.150944 ,-1.902301 },
        { -2.665339, 1.028487, 1.863356 },
        { -0.519059, -1.580847, 0.596605 },
        { -1.535892, -4.004550, 2.329542 },
        { 2.046137, -0.657188 ,-0.847896 },
        { 0.524716, 2.047179, 1.075778 },
        { -0.530676, 1.008563, 1.509182 },
        { 0.710458, -1.426227, 2.217572 }
    };

    for(size_t i = 0; i < N; i++) {
        for (int m = 0; (m < DIM); m++)
        {
            EXPECT_REAL_EQ_TOL(out[XX][XX], expected[XX][XX], gmx::test::defaultRealTolerance());
        }
    }
}
TEST(NBlibTest, generateVelocitySize)
{
    constexpr int N = 10;
    std::vector<real> masses(N, 1.0);
    auto out = generateVelocity(300.0, 1, masses);
    EXPECT_EQ(out.size(), N);
}

TEST(NBlibTest, generateVelocityCheckNumbers)
{
    constexpr int N = 10;
    std::vector<real> masses(N, 1.0);
    auto out = generateVelocity(300.0, 1, masses);
    bool ret = checkNumericValues(out);
    EXPECT_EQ(ret, true);
}

}  // namespace
}  // namespace test
}  // namespace nblib
