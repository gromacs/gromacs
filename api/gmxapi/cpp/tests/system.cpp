/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
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
#include "gmxapi/system.h"

#include <gtest/gtest.h>

#include "programs/mdrun/tests/moduletest.h"

#include "gmxapi/compat/tpr.h"

#include "testingconfiguration.h"

namespace gmxapi
{

namespace testing
{

namespace
{

/*!
 * \brief Check gmxapi::System construction.
 */
TEST_F(GmxApiTest, SystemConstruction)
{
    makeTprFile(1);
    EXPECT_NO_THROW(gmxapi::fromTprFile(runner_.tprFileName_));
    // We have nothing to check at this point other than compilation and
    // error-free execution.
}

/*!
 * \brief Test that OutputControl parameters can be accessed via gmxapi.
 *
 * This test verifies that parameters moved to OutputControl struct are still
 * accessible through the gmxapi parameter interface after the refactoring.
 */
TEST_F(GmxApiTest, OutputControlParameterAccess)
{
    makeTprFile(1);

    // Read TPR file and get parameters
    auto tprReadHandle = gmxapicompat::readTprFile(runner_.tprFileName_);
    ASSERT_NE(tprReadHandle, nullptr);

    auto params = gmxapicompat::getMdParams(*tprReadHandle);
    ASSERT_NE(params, nullptr);

    // Test reading integer OutputControl parameters
    int nstxout = 0;
    EXPECT_NO_THROW(nstxout = gmxapicompat::extractParam(*params, "nstxout", int{}));
    // Default value from test TPR should be 0
    EXPECT_EQ(nstxout, 0);

    int nstvout = 0;
    EXPECT_NO_THROW(nstvout = gmxapicompat::extractParam(*params, "nstvout", int{}));
    EXPECT_EQ(nstvout, 0);

    int nstfout = 0;
    EXPECT_NO_THROW(nstfout = gmxapicompat::extractParam(*params, "nstfout", int{}));
    EXPECT_EQ(nstfout, 0);

    int nstlog = 0;
    EXPECT_NO_THROW(nstlog = gmxapicompat::extractParam(*params, "nstlog", int{}));
    EXPECT_GT(nstlog, 0); // nstlog should have a positive default value

    int nstcalcenergy = 0;
    EXPECT_NO_THROW(nstcalcenergy = gmxapicompat::extractParam(*params, "nstcalcenergy", int{}));
    EXPECT_GT(nstcalcenergy, 0); // nstcalcenergy should have a positive default value

    int nstenergy = 0;
    EXPECT_NO_THROW(nstenergy = gmxapicompat::extractParam(*params, "nstenergy", int{}));
    EXPECT_GT(nstenergy, 0); // nstenergy should have a positive default value

    int nstxoutCompressed = 0;
    EXPECT_NO_THROW(nstxoutCompressed = gmxapicompat::extractParam(*params, "nstxout-compressed", int{}));
    EXPECT_EQ(nstxoutCompressed, 0); // Default is typically 0

    // Test reading real OutputControl parameter
    double compressedXPrecision = 0.0;
    EXPECT_NO_THROW(compressedXPrecision =
                            gmxapicompat::extractParam(*params, "compressed-x-precision", double{}));
    EXPECT_GT(compressedXPrecision, 0.0); // Should have a positive default value

    // Test writing OutputControl parameters
    EXPECT_NO_THROW(gmxapicompat::setParam(params.get(), "nstxout", int64_t{ 500 }));
    EXPECT_NO_THROW(gmxapicompat::setParam(params.get(), "nstvout", int64_t{ 1000 }));
    EXPECT_NO_THROW(gmxapicompat::setParam(params.get(), "nstlog", int64_t{ 250 }));
    EXPECT_NO_THROW(gmxapicompat::setParam(params.get(), "compressed-x-precision", 2000.0));

    // Verify the values were set correctly
    int newNstxout = 0;
    EXPECT_NO_THROW(newNstxout = gmxapicompat::extractParam(*params, "nstxout", int{}));
    EXPECT_EQ(newNstxout, 500);

    int newNstvout = 0;
    EXPECT_NO_THROW(newNstvout = gmxapicompat::extractParam(*params, "nstvout", int{}));
    EXPECT_EQ(newNstvout, 1000);

    int newNstlog = 0;
    EXPECT_NO_THROW(newNstlog = gmxapicompat::extractParam(*params, "nstlog", int{}));
    EXPECT_EQ(newNstlog, 250);

    double newCompressedXPrecision = 0.0;
    EXPECT_NO_THROW(newCompressedXPrecision =
                            gmxapicompat::extractParam(*params, "compressed-x-precision", double{}));
    EXPECT_DOUBLE_EQ(newCompressedXPrecision, 2000.0);
}

} // end anonymous namespace

} // end namespace testing

} // end namespace gmxapi
