/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015,2016, by the GROMACS development team, led by
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
 * Tests for functionality of the "angle" trajectory analysis module.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_applied_forces
 */
#include "gmxpre.h"

#include "gromacs/applied-forces/electricfield.h"

#include <gtest/gtest.h>

#include "gromacs/utility/real.h"

#include "testutils/testasserts.h"

namespace
{

/********************************************************************
 * ElectricFieldTest
 */

class ElectricFieldTest : public ::testing::Test
{
    public:
        ElectricFieldTest() {}
        ~ElectricFieldTest() {}

        void testConstant();
        void testOscillating();
        void testPulsed();
};

void ElectricFieldTest::testConstant()
{
    gmx::test::FloatingPointTolerance tolerance(
            gmx::test::relativeToleranceAsFloatingPoint(1.0, 0.005));
    ElectricField                     efield;
    int  dim    = XX;
    real ampl   = 3.14;
    real omega  = 0;
    real t      = 12;

    efield.setFieldTerm(dim, ampl, omega, 0, 0);
    EXPECT_REAL_EQ_TOL(efield.field(dim, t), ampl, tolerance);
}

void ElectricFieldTest::testOscillating()
{
    gmx::test::FloatingPointTolerance tolerance(
            gmx::test::relativeToleranceAsFloatingPoint(1.0, 0.005));
    ElectricField                     efield;
    int  dim    = YY;
    real ampl   = 3.14;
    real omega  = 13;
    real t      = 12;

    efield.setFieldTerm(dim, ampl, omega, 0, 0);

    EXPECT_REAL_EQ_TOL(efield.field(dim, t), 1.480988, tolerance);
}

void ElectricFieldTest::testPulsed()
{
    gmx::test::FloatingPointTolerance tolerance(
            gmx::test::relativeToleranceAsFloatingPoint(1.0, 0.005));
    ElectricField                     efield;
    int  dim    = YY;
    real ampl   = 3.14;
    real omega  = 13;
    real t0     = 16;
    real sigma  = 6;
    real t      = 12;

    efield.setFieldTerm(dim, ampl, omega, t0, sigma);

    EXPECT_REAL_EQ_TOL(efield.field(dim, t), -0.409810245, tolerance);
}

TEST_F(ElectricFieldTest, Constant)
{
    testConstant();
}

TEST_F(ElectricFieldTest, Oscillating)
{
    testOscillating();
}

TEST_F(ElectricFieldTest, Pulsed)
{
    testPulsed();
}

} // namespace
