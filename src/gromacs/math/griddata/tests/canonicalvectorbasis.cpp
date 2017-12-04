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
 * \brief
 * Tests canonical vector basis
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \ingroup module_math
 */
#include "gmxpre.h"

#include "gromacs/math/griddata/canonicalvectorbasis.h"

#include <vector>
#include <string>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "testutils/testasserts.h"

// #include "testutils/testfilemanager.h"

namespace gmx
{

namespace test
{

namespace internal
{

namespace
{

TEST(CanonicalVectorBasisTest, canConstruct)
{
    CanonicalVectorBasis<1>({1});
    CanonicalVectorBasis<3>({3, -2, 1});
}

TEST(CanonicalVectorBasisTest, throwsWhenZero)
{
    EXPECT_THROW_GMX(CanonicalVectorBasis<1>({0}), gmx::RangeError);
    EXPECT_THROW_GMX(CanonicalVectorBasis<3>({1, 0, 4}), gmx::RangeError);
}

TEST(CanonicalVectorBasisTest, identityMapping)
{
    CanonicalVectorBasis<1>           basis({1});
    CanonicalVectorBasis<1>::NdVector x {
        {
            1
        }
    };
    ASSERT_EQ(x[0], (basis.transformFromBasis(x))[0]);

    CanonicalVectorBasis<3>           basis3d({1, 1, 1});
    CanonicalVectorBasis<3>::NdVector x3d {
        {
            1, 1, 1
        }
    };
    ASSERT_THAT(x3d, testing::ContainerEq(basis3d.transformFromBasis(x3d)));
    ASSERT_THAT(x3d, testing::ContainerEq(basis3d.transformIntoBasis(x3d)));

}

TEST(CanonicalVectorBasisTest, transformIntoIsInverseTransformFrom)
{
    CanonicalVectorBasis<3>           basis3d({-1.5, -1, 7.3});
    CanonicalVectorBasis<3>::NdVector x3d {
        {
            1.1, -5.2, 1.4
        }
    };
    auto transformedIntoAndFrom = basis3d.transformIntoBasis(basis3d.transformFromBasis(x3d));
    auto actual                 = std::begin(transformedIntoAndFrom);
    for (const auto &x : x3d)
    {
        ASSERT_FLOAT_EQ(x, *actual);
        ++actual;
    }
}

TEST(CanonicalVectorBasisTest, transformFromIsInverseTransformInto)
{
    CanonicalVectorBasis<3>           basis3d({-1.5, -1, 7.3});
    CanonicalVectorBasis<3>::NdVector x3d {
        {
            1.1, -5.2, 1.4
        }
    };
    auto transformedIntoAndFrom = basis3d.transformFromBasis(basis3d.transformIntoBasis(x3d));
    auto actual                 = std::begin(transformedIntoAndFrom);
    for (const auto &x : x3d)
    {
        ASSERT_FLOAT_EQ(x, *actual);
        ++actual;
    }
}


} // namespace

} // internal

} // test

} // gmx
