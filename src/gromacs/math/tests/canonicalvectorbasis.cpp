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
 * \brief
 * Tests canonical vector basis
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \ingroup module_math
 */
#include "gmxpre.h"

#include "gromacs/math/canonicalvectorbasis.h"

#include <string>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/inmemoryserializer.h"

#include "testutils/testasserts.h"
#include "testutils/testmatchers.h"

namespace gmx
{

namespace test
{

namespace internal
{

namespace
{

TEST(MdFloatVector, fromRVecIn3D)
{
    FloatingPointTolerance  tolerance(defaultRealTolerance());
    gmx::RVec               r = {-1.5, -1, 7.3};
    MdFloatVector<3>        mdFloat {
        r
    };
    EXPECT_FLOAT_EQ_TOL(r[0], mdFloat[0], tolerance);
    EXPECT_FLOAT_EQ_TOL(r[1], mdFloat[1], tolerance);
    EXPECT_FLOAT_EQ_TOL(r[2], mdFloat[2], tolerance);
}


TEST(CanonicalVectorBasis, canConstruct)
{
    CanonicalVectorBasis<1>({1.7});
    CanonicalVectorBasis<3>({3, -2, 1});
}

TEST(CanonicalVectorBasis, identityMapping)
{
    FloatingPointTolerance  tolerance(defaultFloatTolerance());

    CanonicalVectorBasis<1> basis({1});
    MdFloatVector<1>        x({1});
    EXPECT_EQ(x[0], (basis.transformFromBasis(x))[0]);

    CanonicalVectorBasis<3> basis3d({1, 1, 1});
    MdFloatVector<3>        x3d({ 1, 1, 1 });
    EXPECT_THAT(basis3d.transformFromBasis(x3d), Pointwise(FloatEq(tolerance), x3d));
    EXPECT_THAT(basis3d.transformIntoBasis(x3d), Pointwise(FloatEq(tolerance), x3d));
}

TEST(CanonicalVectorBasis, transformIntoIsInverseTransformFrom)
{
    FloatingPointTolerance  tolerance(defaultFloatTolerance());

    CanonicalVectorBasis<3> basis3d({-1.5, -1, 7.3});
    MdFloatVector<3>        x3d ({ 1.1, -5.2, 1.4 });

    const auto              transformedFromAndInto = basis3d.transformIntoBasis(basis3d.transformFromBasis(x3d));
    EXPECT_THAT(transformedFromAndInto, Pointwise(FloatEq(tolerance), x3d));
}

TEST(CanonicalVectorBasis, transformFromIsInverseTransformInto)
{
    FloatingPointTolerance  tolerance(defaultFloatTolerance());

    CanonicalVectorBasis<3> basis3d({-1.5, -1, 7.3});
    MdFloatVector<3>        x3d({ 1.1, -5.2, 1.4 });

    const auto              transformedIntoAndFrom = basis3d.transformFromBasis(basis3d.transformIntoBasis(x3d));
    EXPECT_THAT(transformedIntoAndFrom, Pointwise(FloatEq(tolerance), x3d));
}

TEST(CanonicalVectorBasis, volumeIsCorrect)
{
    FloatingPointTolerance  tolerance(defaultFloatTolerance());

    CanonicalVectorBasis<4> basis({2, 2, 2, 2});

    EXPECT_FLOAT_EQ_TOL(16, basis.volume(), tolerance);
}

TEST(CanonicalVectorBasis, inverseIsCorrect)
{
    FloatingPointTolerance  tolerance(defaultFloatTolerance());

    CanonicalVectorBasis<4> basis({1, 2, 3, 4});

    EXPECT_FLOAT_EQ_TOL(1, basis.inverse()[0], tolerance);
    EXPECT_FLOAT_EQ_TOL(0.5, basis.inverse()[1], tolerance);
    EXPECT_FLOAT_EQ_TOL(1./3., basis.inverse()[2], tolerance);
    EXPECT_FLOAT_EQ_TOL(0.25, basis.inverse()[3], tolerance);
}

TEST(CanonicalVectorBasis, scaledCopyIsCorrect)
{
    FloatingPointTolerance  tolerance(defaultFloatTolerance());

    CanonicalVectorBasis<4> basis({1, 2, 4, 8});
    const auto              scaledCopy = basis.scaledCopy({{2., 1., 0.5, 0.25}});
    EXPECT_FLOAT_EQ_TOL(2, scaledCopy[0], tolerance);
    EXPECT_FLOAT_EQ_TOL(2, scaledCopy[1], tolerance);
    EXPECT_FLOAT_EQ_TOL(2, scaledCopy[2], tolerance);
    EXPECT_FLOAT_EQ_TOL(2, scaledCopy[3], tolerance);
}

TEST(CanonicalVectorBasis, serializeDeserializeIsIdentity)
{
    FloatingPointTolerance  tolerance(defaultFloatTolerance());

    CanonicalVectorBasis<4> basis({1, 2, 4, 8});

    // serialisation
    gmx::InMemorySerializer serializer;
    basis.serialize(&serializer);

    const auto &buffer = serializer.finishAndGetBuffer();

    // de-serialisation
    gmx::InMemoryDeserializer deserializer(buffer);
    CanonicalVectorBasis<4>   result;

    result.serialize(&deserializer);

    EXPECT_THAT(result.basisVectorLengths(), Pointwise(FloatEq(tolerance), basis.basisVectorLengths()));

}


} // namespace

} // namespace internal

} // namespace test

} // namespace gmx
