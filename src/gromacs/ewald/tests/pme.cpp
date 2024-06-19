/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2021- The GROMACS Authors
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
 * Tests for functionality of the pme related classes:
 * class SeparatePmeRanksPermitted
 *
 * \author Dmitry Morozov <dmitry.morozov@jyu.fi>
 * \ingroup module_ewald
 */
#include "gmxpre.h"

#include "gromacs/ewald/pme.h"

#include <string>

#include <gtest/gtest.h>

#include "gromacs/utility/real.h"
#include "gromacs/utility/stringcompare.h"

#include "testutils/testasserts.h"
#include "testutils/testmatchers.h"

#include "pmetestcommon.h"

namespace gmx
{

namespace test
{

class SeparatePmeRanksPermittedTest : public ::testing::Test
{
public:
    void disableFirstReason() { separatePmeRanksPermitted_.disablePmeRanks("First reason"); }

    void disableSecondReason() { separatePmeRanksPermitted_.disablePmeRanks("Second reason"); }

    void disableEmptyReason() { separatePmeRanksPermitted_.disablePmeRanks(""); }

protected:
    SeparatePmeRanksPermitted separatePmeRanksPermitted_;
};

TEST_F(SeparatePmeRanksPermittedTest, ZeroPmeDisableReasons)
{
    // Expect that SeparatePmeRanksPermitted is enabled by default
    EXPECT_TRUE(separatePmeRanksPermitted_.permitSeparatePmeRanks());
}

TEST_F(SeparatePmeRanksPermittedTest, CanBeDisabled)
{
    // Test if disablePmeRanks works
    EXPECT_NO_THROW(disableFirstReason(););
}

TEST_F(SeparatePmeRanksPermittedTest, OneDisableReasonFlag)
{
    disableFirstReason();

    // Expect that SeparatePmeRanksPermitted is disabled now
    EXPECT_FALSE(separatePmeRanksPermitted_.permitSeparatePmeRanks());
}

TEST_F(SeparatePmeRanksPermittedTest, OneDisableReasonText)
{
    disableFirstReason();

    // Expect that reasonsWhyDisabled works with one reason
    EXPECT_TRUE(separatePmeRanksPermitted_.reasonsWhyDisabled() == "First reason");
}

TEST_F(SeparatePmeRanksPermittedTest, TwoDisableReasonText)
{
    disableFirstReason();
    disableSecondReason();

    // Expect that reasonsWhyDisabled works with two reasons
    EXPECT_TRUE(separatePmeRanksPermitted_.reasonsWhyDisabled() == "First reason; Second reason");
}

TEST_F(SeparatePmeRanksPermittedTest, EmptyDisableReasonText)
{
    disableEmptyReason();

    // Expect that reasonsWhyDisabled works with empty reason
    EXPECT_TRUE(separatePmeRanksPermitted_.reasonsWhyDisabled().empty());
}

} // namespace test

} // namespace gmx
