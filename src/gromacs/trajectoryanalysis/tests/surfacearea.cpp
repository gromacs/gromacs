/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2014- The GROMACS Authors
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
 * Tests for the surface area calculation used by the `sasa` analysis module.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_trajectoryanalysis
 */
#include "gmxpre.h"

#include "gromacs/trajectoryanalysis/modules/surfacearea.h"

#include <cmath>
#include <cstdlib>

#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/random/threefry.h"
#include "gromacs/random/uniformrealdistribution.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"

namespace gmx
{
namespace test
{
namespace
{

/********************************************************************
 * SurfaceAreaTest
 */

class SurfaceAreaTest : public ::testing::Test
{
public:
    SurfaceAreaTest() :
        box_(), rng_(12345), area_(0.0), volume_(0.0), atomArea_(nullptr), dotCount_(0), dots_(nullptr)
    {
    }
    ~SurfaceAreaTest() override
    {
        sfree(atomArea_);
        sfree(dots_);
    }

    void addSphere(real x, real y, real z, real radius, bool bAddToIndex = true)
    {
        if (bAddToIndex)
        {
            index_.push_back(x_.size());
        }
        x_.emplace_back(x, y, z);
        radius_.push_back(radius);
    }

    void generateRandomPosition(rvec x, real* radius)
    {
        rvec                               fx;
        gmx::UniformRealDistribution<real> dist;

        fx[XX] = dist(rng_);
        fx[YY] = dist(rng_);
        fx[ZZ] = dist(rng_);
        mvmul(box_, fx, x);
        *radius = 1.5 * dist(rng_) + 0.5;
    }

    void generateRandomPositions(int count)
    {
        x_.reserve(count);
        radius_.reserve(count);
        index_.reserve(count);
        for (int i = 0; i < count; ++i)
        {
            rvec x;
            real radius = 0;
            generateRandomPosition(x, &radius);
            addSphere(x[XX], x[YY], x[ZZ], radius);
        }
    }
    void translatePoints(real x, real y, real z)
    {
        for (size_t i = 0; i < x_.size(); ++i)
        {
            x_[i][XX] += x;
            x_[i][YY] += y;
            x_[i][ZZ] += z;
        }
    }

    void calculate(int ndots, int flags, bool bPBC)
    {
        volume_ = 0.0;
        sfree(atomArea_);
        atomArea_ = nullptr;
        dotCount_ = 0;
        sfree(dots_);
        dots_ = nullptr;
        t_pbc pbc;
        if (bPBC)
        {
            set_pbc(&pbc, PbcType::Xyz, box_);
        }
        gmx::SurfaceAreaCalculator calculator;
        calculator.setDotCount(ndots);
        calculator.setRadii(radius_);
        calculator.calculate(as_rvec_array(x_.data()),
                             bPBC ? &pbc : nullptr,
                             index_.size(),
                             index_.data(),
                             flags,
                             &area_,
                             &volume_,
                             &atomArea_,
                             &dots_,
                             &dotCount_);
    }
    real resultArea() const { return area_; }
    real resultVolume() const { return volume_; }
    real atomArea(int index) const { return atomArea_[index]; }

    void checkReference(gmx::test::TestReferenceChecker* checker, const char* id, bool checkDotCoordinates)
    {
        gmx::test::TestReferenceChecker compound(checker->checkCompound("SASA", id));
        compound.checkReal(area_, "Area");
        if (volume_ > 0.0)
        {
            compound.checkReal(volume_, "Volume");
        }
        if (atomArea_ != nullptr)
        {
            compound.checkSequenceArray(index_.size(), atomArea_, "AtomArea");
        }
        if (dots_ != nullptr)
        {
            if (checkDotCoordinates)
            {
                // The algorithm may produce the dots in different order in
                // single and double precision due to some internal
                // sorting...
                std::qsort(dots_, dotCount_, sizeof(rvec), &dotComparer);
                compound.checkSequenceArray(3 * dotCount_, dots_, "Dots");
            }
            else
            {
                compound.checkInteger(dotCount_, "DotCount");
            }
        }
    }

    gmx::test::TestReferenceData data_;
    matrix                       box_;

private:
    static int dotComparer(const void* a, const void* b)
    {
        for (int d = DIM - 1; d >= 0; --d)
        {
            const real ad = reinterpret_cast<const real*>(a)[d];
            const real bd = reinterpret_cast<const real*>(b)[d];
            // A fudge factor is needed to get an ordering that is the same
            // in single and double precision, since the points are not
            // exactly on the same Z plane even though in exact arithmetic
            // they probably would be.
            if (ad < bd - 0.001)
            {
                return -1;
            }
            else if (ad > bd + 0.001)
            {
                return 1;
            }
        }
        return 0;
    }

    gmx::DefaultRandomEngine rng_;
    std::vector<gmx::RVec>   x_;
    std::vector<real>        radius_;
    std::vector<int>         index_;

    real  area_;
    real  volume_;
    real* atomArea_;
    int   dotCount_;
    real* dots_;
};

TEST_F(SurfaceAreaTest, ComputesSinglePoint)
{
    gmx::test::FloatingPointTolerance tolerance(gmx::test::defaultRealTolerance());
    addSphere(1, 1, 1, 1);
    ASSERT_NO_FATAL_FAILURE(calculate(24, FLAG_VOLUME | FLAG_ATOM_AREA, false));
    EXPECT_REAL_EQ_TOL(4 * M_PI, resultArea(), tolerance);
    EXPECT_REAL_EQ_TOL(4 * M_PI, atomArea(0), tolerance);
    EXPECT_REAL_EQ_TOL(4 * M_PI / 3, resultVolume(), tolerance);
}

TEST_F(SurfaceAreaTest, ComputesTwoPoints)
{
    gmx::test::FloatingPointTolerance tolerance(gmx::test::relativeToleranceAsFloatingPoint(1.0, 0.005));
    addSphere(1, 1, 1, 1);
    addSphere(2, 1, 1, 1);
    ASSERT_NO_FATAL_FAILURE(calculate(1000, FLAG_ATOM_AREA, false));
    EXPECT_REAL_EQ_TOL(2 * 2 * M_PI * 1.5, resultArea(), tolerance);
    EXPECT_REAL_EQ_TOL(2 * M_PI * 1.5, atomArea(0), tolerance);
    EXPECT_REAL_EQ_TOL(2 * M_PI * 1.5, atomArea(1), tolerance);
}

TEST_F(SurfaceAreaTest, ComputesTwoPointsOfUnequalRadius)
{
    gmx::test::FloatingPointTolerance tolerance(gmx::test::relativeToleranceAsFloatingPoint(1.0, 0.005));
    // Spheres of radius 1 and 2 with intersection at 1.5
    const real dist = 0.5 + std::sqrt(3.25);
    addSphere(1.0, 1.0, 1.0, 1);
    addSphere(1.0 + dist, 1.0, 1.0, 2);
    ASSERT_NO_FATAL_FAILURE(calculate(1000, FLAG_ATOM_AREA, false));
    EXPECT_REAL_EQ_TOL(2 * M_PI * (1.5 + (dist - 0.5 + 2) * 2), resultArea(), tolerance);
    EXPECT_REAL_EQ_TOL(2 * M_PI * 1.5, atomArea(0), tolerance);
    EXPECT_REAL_EQ_TOL(2 * M_PI * (dist - 0.5 + 2) * 2, atomArea(1), tolerance);
}

TEST_F(SurfaceAreaTest, SurfacePoints12)
{
    gmx::test::TestReferenceChecker checker(data_.rootChecker());
    addSphere(0, 0, 0, 1);
    ASSERT_NO_FATAL_FAILURE(calculate(12, FLAG_DOTS, false));
    checkReference(&checker, "Surface", true);
}

TEST_F(SurfaceAreaTest, SurfacePoints32)
{
    gmx::test::TestReferenceChecker checker(data_.rootChecker());
    addSphere(0, 0, 0, 1);
    ASSERT_NO_FATAL_FAILURE(calculate(32, FLAG_DOTS, false));
    checkReference(&checker, "Surface", true);
}

TEST_F(SurfaceAreaTest, SurfacePoints42)
{
    gmx::test::TestReferenceChecker checker(data_.rootChecker());
    addSphere(0, 0, 0, 1);
    ASSERT_NO_FATAL_FAILURE(calculate(42, FLAG_DOTS, false));
    checkReference(&checker, "Surface", true);
}

TEST_F(SurfaceAreaTest, SurfacePoints122)
{
    gmx::test::TestReferenceChecker checker(data_.rootChecker());
    addSphere(0, 0, 0, 1);
    ASSERT_NO_FATAL_FAILURE(calculate(122, FLAG_DOTS, false));
    checkReference(&checker, "Surface", true);
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
    box_[XX][XX]    = 20.0;
    box_[YY][YY]    = 20.0;
    box_[ZZ][ZZ]    = 20.0;
    const int flags = FLAG_ATOM_AREA | FLAG_VOLUME | FLAG_DOTS;
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
    box_[YY][YY] = 10.0 * std::sqrt(3.0);
    box_[ZZ][XX] = 10.0;
    box_[ZZ][YY] = 10.0 * std::sqrt(1.0 / 3.0);
    box_[ZZ][ZZ] = 20.0 * std::sqrt(2.0 / 3.0);

    const int flags = FLAG_ATOM_AREA | FLAG_VOLUME | FLAG_DOTS;
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
} // namespace test
} // namespace gmx
