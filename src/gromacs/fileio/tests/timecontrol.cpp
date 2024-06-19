/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2022- The GROMACS Authors
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
 * Tests for time control value setting.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \ingroup module_fileio
 */
#include "gmxpre.h"

#include "gromacs/fileio/timecontrol.h"

#include <iostream>
#include <optional>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "testutils/include/testutils/testasserts.h"

TEST(TimeControlTest, UnSetHasNoValue)
{
    auto value = timeValue(TimeControl::Begin);
    EXPECT_FALSE(value.has_value());
}

TEST(TimeControlTest, CanSetValue)
{
    setTimeValue(TimeControl::Begin, 13.37);
    auto value = timeValue(TimeControl::Begin);
    ASSERT_TRUE(value.has_value());
    EXPECT_FLOAT_EQ(*value, 13.37);
    auto otherValue = timeValue(TimeControl::End);
    EXPECT_FALSE(otherValue.has_value());
}

TEST(TimeControlTest, CanUnsetValueAgain)
{
    setTimeValue(TimeControl::Begin, 13.37);
    setTimeValue(TimeControl::End, 42.23);
    auto value      = timeValue(TimeControl::Begin);
    auto otherValue = timeValue(TimeControl::End);
    EXPECT_TRUE(value.has_value());
    EXPECT_TRUE(otherValue.has_value());
    unsetTimeValue(TimeControl::Begin);
    auto newValue      = timeValue(TimeControl::Begin);
    auto newOtherValue = timeValue(TimeControl::End);
    EXPECT_FALSE(newValue.has_value());
    EXPECT_TRUE(newOtherValue.has_value());
}
