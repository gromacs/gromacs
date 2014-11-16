/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014, by the GROMACS development team, led by
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
 * Tests for the surface area calculation used by the `sasa` analysis module.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_trajectoryanalysis
 */
#include "gmxpre.h"

#include "gromacs/trajectoryanalysis/modules/surfacearea.h"

#include <gtest/gtest.h>

#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/random/random.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"

namespace
{

/********************************************************************
 * SurfaceAreaTest
 */

class SurfaceAreaTest : public ::testing::Test
{
    public:
        SurfaceAreaTest()
            : rng_(NULL), allocated_(0), x_(NULL), area_(0.0), volume_(0.0),
              atomArea_(NULL), dotCount_(0), dots_(NULL)
        {
            // TODO: Handle errors.
            rng_ = gmx_rng_init(12345);
            clear_mat(box_);
        }
        ~SurfaceAreaTest()
        {
            if (rng_ != NULL)
            {
                gmx_rng_destroy(rng_);
            }
            sfree(x_);
            sfree(atomArea_);
            sfree(dots_);
        }

        void reserveSpace(int count)
        {
            GMX_RELEASE_ASSERT(allocated_ == 0 && x_ == NULL,
                               "Cannot allocate data more than once");
            snew(x_, count);
            radius_.reserve(count);
            index_.reserve(count);
            allocated_ = count;
        }

        void addSphere(real x, real y, real z, real radius,
                       bool bAddToIndex = true)
        {
            const int index = radius_.size();
            GMX_RELEASE_ASSERT(index < allocated_,
                               "Not enough space allocated for test data");
            x_[index][XX] = x;
            x_[index][YY] = y;
            x_[index][ZZ] = z;
            if (bAddToIndex)
            {
                index_.push_back(radius_.size());
            }
            radius_.push_back(radius);
        }

        void generateRandomPosition(rvec x, real *radius)
        {
            rvec fx;
            fx[XX]  = gmx_rng_uniform_real(rng_);
            fx[YY]  = gmx_rng_uniform_real(rng_);
            fx[ZZ]  = gmx_rng_uniform_real(rng_);
            mvmul(box_, fx, x);
            *radius = 1.5*gmx_rng_uniform_real(rng_) + 0.5;
        }

        void addDummySpheres(int count)
        {
            for (int i = 0; i < count; ++i)
            {
                rvec x;
                real radius;
                generateRandomPosition(x, &radius);
                addSphere(x[XX], x[YY], x[ZZ], radius, false);
            }
        }

        void generateRandomPositions(int count)
        {
            reserveSpace(count);
            for (int i = 0; i < count; ++i)
            {
                rvec x;
                real radius;
                generateRandomPosition(x, &radius);
                addSphere(x[XX], x[YY], x[ZZ], radius);
            }
        }
        void translatePoints(real x, real y, real z)
        {
            for (size_t i = 0; i < radius_.size(); ++i)
            {
                x_[i][XX] += x;
                x_[i][YY] += y;
                x_[i][ZZ] += z;
            }
        }

        void calculate(int ndots, int flags, bool bPBC)
        {
            volume_   = 0.0;
            sfree(atomArea_);
            atomArea_ = NULL;
            dotCount_ = 0;
            sfree(dots_);
            dots_     = NULL;
            ASSERT_EQ(0, nsc_dclm_pbc(x_, &radius_[0], index_.size(), ndots, flags,
                                      &area_, &atomArea_, &volume_, &dots_, &dotCount_,
                                      &index_[0], epbcXYZ, bPBC ? box_ : NULL));
        }
        real resultArea() const { return area_; }
        real resultVolume() const { return volume_; }
        real atomArea(int index) const { return atomArea_[index]; }

        void checkReference(gmx::test::TestReferenceChecker *checker, const char *id,
                            bool checkDotCoordinates)
        {
            gmx::test::TestReferenceChecker compound(
                    checker->checkCompound("SASA", id));
            compound.checkReal(area_, "Area");
            if (volume_ > 0.0)
            {
                compound.checkReal(volume_, "Volume");
            }
            if (atomArea_ != NULL)
            {
                compound.checkSequenceArray(index_.size(), atomArea_, "AtomArea");
            }
            if (dots_ != NULL)
            {
                if (checkDotCoordinates)
                {
                    compound.checkSequenceArray(3*dotCount_, dots_, "Dots");
                }
                else
                {
                    compound.checkInteger(dotCount_, "DotCount");
                }
            }
        }

        gmx::test::TestReferenceData    data_;
        matrix                          box_;

    private:
        gmx_rng_t          rng_;
        int                allocated_;
        rvec              *x_;
        std::vector<real>  radius_;
        std::vector<int>   index_;

        real               area_;
        real               volume_;
        real              *atomArea_;
        int                dotCount_;
        real              *dots_;
};

TEST_F(SurfaceAreaTest, ComputesSinglePoint)
{
    gmx::test::FloatingPointTolerance tolerance(
            gmx::test::defaultRealTolerance());
    reserveSpace(1);
    addSphere(1, 1, 1, 1);
    ASSERT_NO_FATAL_FAILURE(calculate(24, FLAG_VOLUME | FLAG_ATOM_AREA, false));
    EXPECT_REAL_EQ_TOL(4*M_PI, resultArea(), tolerance);
    EXPECT_REAL_EQ_TOL(4*M_PI, atomArea(0), tolerance);
    EXPECT_REAL_EQ_TOL(4*M_PI/3, resultVolume(), tolerance);
}

TEST_F(SurfaceAreaTest, ComputesTwoPoints)
{
    gmx::test::FloatingPointTolerance tolerance(
            gmx::test::relativeToleranceAsFloatingPoint(1.0, 0.005));
    reserveSpace(2);
    addSphere(1, 1, 1, 1);
    addSphere(2, 1, 1, 1);
    ASSERT_NO_FATAL_FAILURE(calculate(1000, FLAG_ATOM_AREA, false));
    EXPECT_REAL_EQ_TOL(2*2*M_PI*1.5, resultArea(), tolerance);
    EXPECT_REAL_EQ_TOL(2*M_PI*1.5, atomArea(0), tolerance);
    EXPECT_REAL_EQ_TOL(2*M_PI*1.5, atomArea(1), tolerance);
}

TEST_F(SurfaceAreaTest, ComputesTwoPointsOfUnequalRadius)
{
    gmx::test::FloatingPointTolerance tolerance(
            gmx::test::relativeToleranceAsFloatingPoint(1.0, 0.005));
    reserveSpace(2);
    // Spheres of radius 1 and 2 with intersection at 1.5
    const real dist = 0.5 + sqrt(3.25);
    addSphere(1.0, 1.0, 1.0, 1);
    addSphere(1.0 + dist, 1.0, 1.0, 2);
    ASSERT_NO_FATAL_FAILURE(calculate(1000, FLAG_ATOM_AREA, false));
    EXPECT_REAL_EQ_TOL(2*M_PI*(1.5 + (dist - 0.5 + 2)*2), resultArea(), tolerance);
    EXPECT_REAL_EQ_TOL(2*M_PI*1.5, atomArea(0), tolerance);
    EXPECT_REAL_EQ_TOL(2*M_PI*(dist - 0.5 + 2)*2, atomArea(1), tolerance);
}

TEST_F(SurfaceAreaTest, Computes100Points)
{
    gmx::test::TestReferenceChecker checker(data_.rootChecker());
    checker.setDefaultTolerance(gmx::test::absoluteTolerance(0.001));
    box_[XX][XX] = 10.0;
    box_[YY][YY] = 10.0;
    box_[ZZ][ZZ] = 10.0;
    generateRandomPositions(100);
    ASSERT_NO_FATAL_FAILURE(calculate(24, FLAG_VOLUME | FLAG_ATOM_AREA | FLAG_DOTS, false));
    checkReference(&checker, "100Points", false);
}

TEST_F(SurfaceAreaTest, Computes100PointsWithRectangularPBC)
{
    // TODO: It would be nice to check that this produces the same result as
    // without PBC, without duplicating the reference files.
    gmx::test::TestReferenceChecker checker(data_.rootChecker());
    checker.setDefaultTolerance(gmx::test::absoluteTolerance(0.001));
    box_[XX][XX] = 10.0;
    box_[YY][YY] = 10.0;
    box_[ZZ][ZZ] = 10.0;
    generateRandomPositions(100);
    box_[XX][XX] = 20.0;
    box_[YY][YY] = 20.0;
    box_[ZZ][ZZ] = 20.0;
    // TODO: The volume computation with PBC is broken...
    const int flags = FLAG_ATOM_AREA | FLAG_DOTS;
    ASSERT_NO_FATAL_FAILURE(calculate(24, flags, true));
    checkReference(&checker, "100Points", false);

    translatePoints(15.0, 0, 0);
    ASSERT_NO_FATAL_FAILURE(calculate(24, flags, true));
    checkReference(&checker, "100Points", false);

    translatePoints(-15.0, 15.0, 0);
    ASSERT_NO_FATAL_FAILURE(calculate(24, flags, true));
    checkReference(&checker, "100Points", false);

    translatePoints(0, -15.0, 15.0);
    ASSERT_NO_FATAL_FAILURE(calculate(24, flags, true));
    checkReference(&checker, "100Points", false);
}

TEST_F(SurfaceAreaTest, Computes100PointsWithTriclinicPBC)
{
    // TODO: It would be nice to check that this produces the same result as
    // without PBC, without duplicating the reference files.
    gmx::test::TestReferenceChecker checker(data_.rootChecker());
    checker.setDefaultTolerance(gmx::test::absoluteTolerance(0.001));
    box_[XX][XX] = 10.0;
    box_[YY][YY] = 10.0;
    box_[ZZ][ZZ] = 10.0;
    generateRandomPositions(100);
    box_[XX][XX] = 20.0;
    box_[YY][XX] = 10.0;
    box_[YY][YY] = 10.0*sqrt(3.0);
    box_[ZZ][XX] = 10.0;
    box_[ZZ][YY] = 10.0*sqrt(1.0/3.0);
    box_[ZZ][ZZ] = 20.0*sqrt(2.0/3.0);

    // TODO: The volume computation with PBC is broken...
    const int flags = FLAG_ATOM_AREA | FLAG_DOTS;
    ASSERT_NO_FATAL_FAILURE(calculate(24, flags, true));
    checkReference(&checker, "100Points", false);

    translatePoints(15.0, 0, 0);
    ASSERT_NO_FATAL_FAILURE(calculate(24, flags, true));
    checkReference(&checker, "100Points", false);

    translatePoints(-15.0, box_[YY][YY] - 5.0, 0);
    ASSERT_NO_FATAL_FAILURE(calculate(24, flags, true));
    checkReference(&checker, "100Points", false);

    translatePoints(0, -(box_[YY][YY] - 5.0), 15.0);
    ASSERT_NO_FATAL_FAILURE(calculate(24, flags, true));
    checkReference(&checker, "100Points", false);
}

} // namespace
