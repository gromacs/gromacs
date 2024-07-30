/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2016- The GROMACS Authors
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
 * Implements test of some pulling routines
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_pulling
 */
#include "gmxpre.h"

#include "gromacs/pulling/pull.h"

#include "config.h"

#include <cmath>

#include <algorithm>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/math/functions.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/pull_params.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"

#if HAVE_MUPARSER
#    include <muParser.h>
#endif

#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pulling/pull_internal.h"
#include "gromacs/pulling/transformationcoordinate.h"
#include "gromacs/utility/smalloc.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace gmx
{
namespace test
{
namespace
{

using gmx::test::defaultRealTolerance;

class PullTest : public ::testing::Test
{
protected:
    PullTest() {}

    static void test(PbcType pbcType, matrix box)
    {
        t_pbc pbc;

        // PBC stuff
        set_pbc(&pbc, pbcType, box);

        GMX_ASSERT(pbc.ndim_ePBC >= 1 && pbc.ndim_ePBC <= DIM,
                   "Tests only support PBC along at least x and at most x, y, and z");

        real boxSizeZSquared;
        if (pbc.ndim_ePBC > ZZ)
        {
            boxSizeZSquared = gmx::square(box[ZZ][ZZ]);
        }
        else
        {
            boxSizeZSquared = GMX_REAL_MAX;
        }

        {
            // Distance pulling in all 3 dimensions
            t_pull_coord params;
            params.eGeom      = PullGroupGeometry::Distance;
            params.dim[XX]    = 1;
            params.dim[YY]    = 1;
            params.dim[ZZ]    = 1;
            params.coordIndex = 0;
            pull_coord_work_t pcrd(params);
            clear_dvec(pcrd.spatialData.vec);

            real minBoxSize2 = GMX_REAL_MAX;
            for (int d = 0; d < pbc.ndim_ePBC; d++)
            {
                minBoxSize2 = std::min(minBoxSize2, norm2(box[d]));
            }
            EXPECT_REAL_EQ_TOL(0.25 * minBoxSize2, max_pull_distance2(pcrd, pbc), defaultRealTolerance());
        }

        {
            // Distance pulling along Z
            t_pull_coord params;
            params.eGeom      = PullGroupGeometry::Distance;
            params.dim[XX]    = 0;
            params.dim[YY]    = 0;
            params.dim[ZZ]    = 1;
            params.coordIndex = 0;
            pull_coord_work_t pcrd(params);
            clear_dvec(pcrd.spatialData.vec);
            EXPECT_REAL_EQ_TOL(
                    0.25 * boxSizeZSquared, max_pull_distance2(pcrd, pbc), defaultRealTolerance());
        }

        {
            // Directional pulling along Z
            t_pull_coord params;
            params.eGeom      = PullGroupGeometry::Direction;
            params.dim[XX]    = 1;
            params.dim[YY]    = 1;
            params.dim[ZZ]    = 1;
            params.coordIndex = 0;
            pull_coord_work_t pcrd(params);
            clear_dvec(pcrd.spatialData.vec);
            pcrd.spatialData.vec[ZZ] = 1;
            EXPECT_REAL_EQ_TOL(
                    0.25 * boxSizeZSquared, max_pull_distance2(pcrd, pbc), defaultRealTolerance());
        }

        {
            // Directional pulling along X
            t_pull_coord params;
            params.eGeom      = PullGroupGeometry::Direction;
            params.dim[XX]    = 1;
            params.dim[YY]    = 1;
            params.dim[ZZ]    = 1;
            params.coordIndex = 0;
            pull_coord_work_t pcrd(params);
            clear_dvec(pcrd.spatialData.vec);
            pcrd.spatialData.vec[XX] = 1;

            real minDist2 = square(box[XX][XX]);
            for (int d = XX + 1; d < DIM; d++)
            {
                minDist2 -= square(box[d][XX]);
            }
            EXPECT_REAL_EQ_TOL(0.25 * minDist2, max_pull_distance2(pcrd, pbc), defaultRealTolerance());
        }
    }
};

TEST_F(PullTest, MaxPullDistanceXyzScrewBox)
{
    matrix box = { { 10, 0, 0 }, { 0, 10, 0 }, { 0, 0, 10 } };

    test(PbcType::Screw, box);
}

TEST_F(PullTest, MaxPullDistanceXyzCubicBox)
{
    matrix box = { { 10, 0, 0 }, { 0, 10, 0 }, { 0, 0, 10 } };

    test(PbcType::Xyz, box);
}

TEST_F(PullTest, MaxPullDistanceXyzTricBox)
{
    matrix box = { { 10, 0, 0 }, { 3, 10, 0 }, { 3, 4, 10 } };

    test(PbcType::Xyz, box);
}

TEST_F(PullTest, MaxPullDistanceXyzLongBox)
{
    matrix box = { { 10, 0, 0 }, { 0, 10, 0 }, { 0, 0, 30 } };

    test(PbcType::Xyz, box);
}

TEST_F(PullTest, MaxPullDistanceXySkewedBox)
{
    matrix box = { { 10, 0, 0 }, { 5, 8, 0 }, { 0, 0, 0 } };

    test(PbcType::XY, box);
}

#if HAVE_MUPARSER

/*
 * Simple test case with one transformation coordinate and one pull coordinate.
 *
 */
TEST_F(PullTest, TransformationCoordSimple)
{
    //-----_SETUP-------
    pull_t pull;
    // Create standard pull coordinates
    const double dist = 10; // the pull coord value
    t_pull_coord x1;
    x1.coordIndex = 0;
    pull.coord.emplace_back(x1);
    pull.coord[0].spatialData.value = dist;

    // Create transformation pull coordinates
    t_pull_coord x2;
    x2.eGeom      = PullGroupGeometry::Transformation;
    x2.expression = "-x1";
    x2.coordIndex = 1;
    pull.coord.emplace_back(x2);

    const double t = 0.0;

    //-----TESTS -------
    // 1) check transformation pull coordinate values
    pull.coord[1].spatialData.value = getTransformationPullCoordinateValue(
            &pull.coord[1], constArrayRefFromArray(pull.coord.data(), 1), t);
    EXPECT_REAL_EQ_TOL(pull.coord[1].spatialData.value, -dist, defaultRealTolerance());

    // 2) check that a force on the transformation coord propagates to the original pull coordinate
    const double force        = 10; // The force on the transformation coordinate
    pull.coord[1].scalarForce = force;
    distributeTransformationPullCoordForce(
            &pull.coord[1], gmx::ArrayRef<pull_coord_work_t>(pull.coord).subArray(0, 1));
    EXPECT_REAL_EQ_TOL(force, pull.coord[1].scalarForce, defaultRealTolerance());
    EXPECT_REAL_EQ_TOL(-force,
                       pull.coord[0].scalarForce,
                       // Since we do numerical differentiation we need to increase the tolerance
                       test::relativeToleranceAsFloatingPoint(force, 1e-5));
}

/*
 * More complicated test case of transformation pull coordinates.
 * To account for all the different ways transformation coordinates can depend on each other,
 * we need a relatively complicated setup with multiple pull coordinates.
 * We also test the code for multiple numerical values, which adds a for loop.
 */
TEST_F(PullTest, TransformationCoordAdvanced)
{
    pull_t pull;
    // Create standard pull coordinates
    t_pull_coord x1;
    x1.eGeom      = PullGroupGeometry::Distance;
    x1.coordIndex = 0;
    pull.coord.emplace_back(x1);
    t_pull_coord x2;
    x2.eGeom      = PullGroupGeometry::Angle;
    x2.coordIndex = 1;
    pull.coord.emplace_back(x2);

    // Create transformation pull coordinates
    // x3, a pull coordinate that depends on another pull coordinate
    t_pull_coord x3;
    x3.eGeom                = PullGroupGeometry::Transformation;
    std::string expression1 = "x1^2";
    x3.expression           = expression1;
    x3.coordIndex           = 2;
    pull.coord.emplace_back(x3);

    // x4, the last transformation pull coordinate
    t_pull_coord x4;
    x4.eGeom                = PullGroupGeometry::Transformation;
    std::string expression2 = "x1 - 0.5*x2^3 + x3^2 + 3"; // note that x3^2 is equivalent to x1^4
    x4.expression           = expression2;
    x4.coordIndex           = 3;
    x4.dx                   = 1e-4;
    pull.coord.emplace_back(x4);

    const double t = 0.0;

    // below we set x1 and x2 to different values and make sure that
    // 1) the transformation coordinates are correct, i.e. test getTransformationPullCoordinateValue
    // 2) that the force is accurately distributed from the transformation coord to the normal
    // pull coordinates, i.e. test applyTransformationPullCoordForce
    for (double v1 = 0.5; v1 < 11; v1 += 1)
    {
        clear_pull_forces(&pull);

        double v2 = -v1 * 10;
        // transformation pull coord value
        pull.coord[0].spatialData.value = v1;
        pull.coord[1].spatialData.value = v2;
        pull.coord[2].spatialData.value = getTransformationPullCoordinateValue(
                &pull.coord[2], constArrayRefFromArray(pull.coord.data(), 2), t);
        pull.coord[3].spatialData.value = getTransformationPullCoordinateValue(
                &pull.coord[3], constArrayRefFromArray(pull.coord.data(), 3), t);

        // 1) check transformation pull coordinate values
        // Since we perform numerical differentiation and floating point operations
        // we only expect the results below to be approximately equal
        double expectedX3 = v1 * v1;
        EXPECT_REAL_EQ_TOL(pull.coord[2].spatialData.value, expectedX3, defaultRealTolerance());

        double expectedX4 = v1 - 0.5 * v2 * v2 * v2 + expectedX3 * expectedX3 + 3;
        EXPECT_REAL_EQ_TOL(pull.coord[3].spatialData.value, expectedX4, defaultRealTolerance());

        double transformationForcex3 = v2 / 9 + 1;
        pull.coord[2].scalarForce += transformationForcex3;
        double transformationForcex4 = v1 + 4.5;
        pull.coord[3].scalarForce    = transformationForcex4;

        // Distribute the transformation coordinate forces
        pull.comm.bParticipate = true;
        pull_apply_forces(&pull, {}, nullptr, nullptr);

        // 2) check forces on transformation and normal pull coordinates with a direct relationship
        double tolX4 = 1e-2;

        // Note: sum of direct force from pull.coord[2] and indirect force of pull.coord[3]
        double expectedFx1 =
                transformationForcex3 * 2 * v1 + transformationForcex4 * (1 + 4 * v1 * v1 * v1);
        EXPECT_REAL_EQ_TOL(expectedFx1,
                           pull.coord[0].scalarForce,
                           test::relativeToleranceAsFloatingPoint(expectedFx1, tolX4));

        double expectedFx2 = -1.5 * v2 * v2 * transformationForcex4;
        EXPECT_REAL_EQ_TOL(expectedFx2,
                           pull.coord[1].scalarForce,
                           test::relativeToleranceAsFloatingPoint(expectedFx2, tolX4));

        // Note: sum of direct force plus the force of the transformation from pull.coord[3]
        double expectedFx3 = transformationForcex3 + 2 * expectedX3 * transformationForcex4;
        EXPECT_REAL_EQ_TOL(expectedFx3,
                           pull.coord[2].scalarForce,
                           test::relativeToleranceAsFloatingPoint(expectedFx3, tolX4));

        double expectedFx4 = transformationForcex4;
        EXPECT_REAL_EQ_TOL(expectedFx4, pull.coord[3].scalarForce, defaultRealTolerance());
    }
}

/*
 * Tests that time works corretly as a variable in transformation coordinates
 *
 */
TEST_F(PullTest, TransformationCoordTime)
{
    //-----_SETUP-------
    pull_t pull;
    // Create standard pull coordinates
    const double dist = 2; // the pull coord value
    t_pull_coord x1;
    x1.coordIndex = 0;
    pull.coord.emplace_back(x1);
    pull.coord[0].spatialData.value = dist;

    // Create a transformation pull coordinate using time
    t_pull_coord x2;
    x2.eGeom      = PullGroupGeometry::Transformation;
    x2.expression = "-x1 + sin(2 * t)";
    x2.coordIndex = 1;
    pull.coord.emplace_back(x2, true);

    const double t = 0.75;

    const double ref = -dist + std::sin(2 * t);

    //-----TESTS -------
    // 1) check transformation pull coordinate values
    pull.coord[1].spatialData.value = getTransformationPullCoordinateValue(
            &pull.coord[1], constArrayRefFromArray(pull.coord.data(), 1), t);
    EXPECT_REAL_EQ_TOL(pull.coord[1].spatialData.value, ref, defaultRealTolerance());
}

/*
 * Tests that pull transformation throws when time is not allowed as a variable
 *
 */
TEST_F(PullTest, TransformationCoordTimeNotAllowed)
{
    //-----_SETUP-------
    pull_t pull;
    // Create standard pull coordinates
    const double dist = 2; // the pull coord value
    t_pull_coord x1;
    x1.coordIndex = 0;
    pull.coord.emplace_back(x1);
    pull.coord[0].spatialData.value = dist;

    // Create a transformation pull coordinate using time
    t_pull_coord x2;
    x2.eGeom      = PullGroupGeometry::Transformation;
    x2.expression = "-x1 + sin(2 * t)";
    x2.coordIndex = 1;
    // Pass false to not allow time as a variable in transformation coordinate
    pull.coord.emplace_back(x2, false);

    const double t = 0.75;

    EXPECT_THROW_GMX(getTransformationPullCoordinateValue(
                             &pull.coord[1], constArrayRefFromArray(pull.coord.data(), 1), t),
                     gmx::InconsistentInputError);
}

/*
 * Simple test case with one transformation coordinate set to a constant value
 *
 * The purpose of this test is just to make sure that the code works even in this case
 *
 */
TEST_F(PullTest, TransformationCoordDummyExpression)
{
    //-----_SETUP-------
    pull_t       pull;
    t_pull_coord x;
    x.eGeom      = PullGroupGeometry::Transformation;
    x.coordIndex = 0;
    x.expression = "10";
    pull.coord.emplace_back(x);
    //-----TESTS -------
    // check transformation pull coordinate values
    double value = getTransformationPullCoordinateValue(
            &pull.coord[0], ArrayRef<const pull_coord_work_t>{}, 0.0);
    EXPECT_REAL_EQ_TOL(value, 10, defaultRealTolerance());
}
#endif // HAVE_MUPARSER

} // namespace
} // namespace test
} // namespace gmx
