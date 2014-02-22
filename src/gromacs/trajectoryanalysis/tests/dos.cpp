/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013,2014,2016, by the GROMACS development team, led by
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
 * \ingroup module_trajectoryanalysis
 */
#include "gmxpre.h"

#include "gromacs/trajectoryanalysis/modules/dos.h"

#include <gtest/gtest.h>

#include "gromacs/math/functions.h"
#include "gromacs/trajectoryanalysis/modules/dosutils.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"

#include "moduletest.h"

namespace
{

using gmx::test::CommandLine;

/********************************************************************
 * Tests for gmx::analysismodules::Dos.
 */

//! Test fixture for the angle analysis module.
typedef gmx::test::TrajectoryAnalysisModuleTestFixture<gmx::analysismodules::DosInfo>
    DosModuleTest;

TEST_F(DosModuleTest, ComputesFD)
{
    gmx::test::FloatingPointTolerance tolerance(
            gmx::test::defaultRealTolerance());
    EXPECT_REAL_EQ_TOL(-1.9799999, FD(100, 0.01), tolerance);
    EXPECT_REAL_EQ_TOL(-1.799981, FD(100, 0.1), tolerance);
    EXPECT_REAL_EQ_TOL(-0.19608453, FD(100, 0.9), tolerance);
    EXPECT_REAL_EQ_TOL(-0.01512, FD(100, 0.99), tolerance);
}

TEST_F(DosModuleTest, ComputesYYY)
{
    gmx::test::FloatingPointTolerance tolerance(
            gmx::test::defaultRealTolerance());
    EXPECT_REAL_EQ_TOL(-1.95019975, YYY(0.01, 0.5), tolerance);
    EXPECT_REAL_EQ_TOL(-1.467168, YYY(0.1,  0.6), tolerance);
    EXPECT_REAL_EQ_TOL(1.131694, YYY(0.9, 0.7), tolerance);
    EXPECT_REAL_EQ_TOL(1.17792213, YYY(0.99, 0.8), tolerance);
}

TEST_F(DosModuleTest, ComputesHardSphereEntropy)
{
    gmx::test::FloatingPointTolerance tolerance(
            gmx::test::defaultRealTolerance());
    EXPECT_REAL_EQ_TOL(-0.020246137, calc_Shs(0.5), tolerance);
    EXPECT_REAL_EQ_TOL(-0.041114569, calc_Shs(0.6), tolerance);
    EXPECT_REAL_EQ_TOL(-0.0877367, calc_Shs(0.7), tolerance);
    EXPECT_REAL_EQ_TOL(-0.2204596, calc_Shs(0.8), tolerance);
}

TEST_F(DosModuleTest, ComputesCompressibility)
{
    gmx::test::FloatingPointTolerance tolerance(
            gmx::test::defaultRealTolerance());

    for (int i = 1; i < 10; i++)
    {
        double y        = i * 0.1;
        double compress = calc_compress(y);
        EXPECT_REAL_EQ_TOL((1+y+gmx::square(y)-gmx::power3(y))/(gmx::power3(1-y)),
                           compress, tolerance);
    }
}

TEST_F(DosModuleTest, ComputesFluidicity)
{
    gmx::test::FloatingPointTolerance tolerance(
            gmx::test::defaultRealTolerance());

    double toler = 1e-4;
    EXPECT_REAL_EQ_TOL(0.5411987305, calc_fluidicity(1, toler), tolerance);
    EXPECT_REAL_EQ_TOL(0.9345093, calc_fluidicity(10, toler), tolerance);
    EXPECT_REAL_EQ_TOL(0.9974976, calc_fluidicity(100, toler), tolerance);
    EXPECT_REAL_EQ_TOL(0.18695068, calc_fluidicity(0.1, toler), tolerance);

}

} // namespace
