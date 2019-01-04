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
 * Tests coordinate transformations.
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \ingroup module_math
 */
#include "gmxpre.h"

#include "gromacs/math/coordinatetransformations.h"

#include <array>

#include <gtest/gtest.h>

#include "testutils/testasserts.h"
#include "testutils/testmatchers.h"

using ::testing::Pointwise;

namespace gmx
{

namespace test
{

namespace
{

class TransformationTest : public ::testing::Test
{
    public:
        const std::array<float, 2> input       = {4, 5};
        const std::array<float, 2> scale       = {2, 3};
        const std::array<float, 2> translation = {7, 6};
        const                      TranslationTransform < std::array < float, 2>> translationTransform {translation};
        const                      ScaleTransform < std::array < float, 2>> scaleTransform {scale};
        FloatingPointTolerance     tolerance {defaultFloatTolerance()};
};

TEST_F(TransformationTest, TranslationAndInverse)
{
    // testing translating (4,5) + (7,6)
    std::array<float, 2> expected = {11, 11};
    EXPECT_THAT(translationTransform(input), Pointwise(FloatEq(tolerance), expected));
    // inverse transformation of translated result shall be input
    EXPECT_THAT(translationTransform.inverse(expected), Pointwise(FloatEq(tolerance), input));
}

TEST_F(TransformationTest, ScaleAndInverse)
{
    // testing scaling (4,5) with (2,3) = (2*4, 3*5)
    std::array<float, 2> expected = {8, 15};
    EXPECT_THAT(scaleTransform(input), Pointwise(FloatEq(tolerance), expected));
    // inverse transformation of scaled result shall be input
    EXPECT_THAT(scaleTransform.inverse(expected), Pointwise(FloatEq(tolerance), input));
}

TEST_F(TransformationTest, CompositeScaleTranslateAndInverse)
{
    // testing scaling (4,5) followed up with translation (7,6)
    auto                 compositeTransform = composeTransforms(scaleTransform, translationTransform);
    std::array<float, 2> expected           = {15, 21};
    EXPECT_THAT(compositeTransform(input), Pointwise(FloatEq(tolerance), expected));
    // inverse composite transformation of result shall be input
    EXPECT_THAT(compositeTransform.inverse(expected), Pointwise(FloatEq(tolerance), input));
}

TEST_F(TransformationTest, CompositeTranslateScaleAndInverse)
{
    // testing translation followed up with scaling
    auto                 compositeTransform = composeTransforms(translationTransform, scaleTransform);
    std::array<float, 2> expected           = {22, 33};
    EXPECT_THAT(compositeTransform(input), Pointwise(FloatEq(tolerance), expected));
    EXPECT_THAT(compositeTransform.inverse(expected), Pointwise(FloatEq(tolerance), input));
}

TEST_F(TransformationTest, MultipleComposites)
{
    // testing translation followed up with scaling
    auto                 compositeTransform   = composeTransforms(translationTransform, scaleTransform);
    auto                 composeWithComposite = composeTransforms(compositeTransform, translationTransform);
    std::array<float, 2> expected             = {29, 39};
    EXPECT_THAT(composeWithComposite(input), Pointwise(FloatEq(tolerance), expected));
    EXPECT_THAT(composeWithComposite.inverse(expected), Pointwise(FloatEq(tolerance), input));
}

} // namespace

} // namespace test

} // namespace gmx
