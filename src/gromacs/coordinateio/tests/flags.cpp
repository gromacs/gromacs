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
/*!\file
 * \internal
 * \brief
 * Tests for flag setting method
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \ingroup module_coordinateio
 */

#include "gmxpre.h"

#include "flags.h"

#include "gromacs/coordinateio/outputadapters/outputselector.h"
#include "gromacs/utility/arrayref.h"

namespace gmx
{

namespace test
{

TEST_F(FlagTest, EmptyFlagsWork)
{
    OutputAdapters adapters = getModules();
    EXPECT_EQ(adapters.size(), 0);
    EXPECT_NO_THROW(runTest("test.tng", std::move(adapters)));
}

TEST_F(FlagTest, RunsWithSelection)
{
    Selection sel;
    ASSERT_NO_THROW(addOptionForSelection(&sel, true));
    ASSERT_NO_THROW(setSelectionOptionValues(getOption(), &sel, true));
    OutputAdapters adapters = getModules();
    adapters.emplace_back(compat::make_unique<OutputSelector>(sel));
    EXPECT_EQ(adapters.size(), 1);
    EXPECT_NO_THROW(runTest("test.tng", std::move(adapters)));
}

TEST_F(FlagTest, AddsModuleWhenNeeded)
{
    std::string              option = "atoms";
    std::string              value  = "yes";
    setModuleFlag(option, value, getOption(), false);
    OutputAdapters           adapters = getModules();
    EXPECT_EQ(adapters.size(), 1);
    EXPECT_NO_THROW(runTest("test.tng", std::move(adapters)));
}

TEST_F(FlagTest, DoesntAddModuleWhenNotNeeded)
{
    std::string              option = "atoms";
    std::string              value  = "unchanged";
    setModuleFlag(option, value, getOption(), false);
    OutputAdapters           adapters = getModules();
    EXPECT_EQ(adapters.size(), 0);
    EXPECT_NO_THROW(runTest("test.tng", std::move(adapters)));
}

TEST_F(FlagTest, CanAddNewBox)
{
    std::string              option = "newbox";
    std::string              value  = "3 3 3";
    setModuleFlag(option, value, getOption(), true);
    OutputAdapters           adapters = getModules();
    EXPECT_EQ(adapters.size(), 1);
    EXPECT_NO_THROW(runTest("test.tng", std::move(adapters)));
}

TEST_F(FlagTest, CannotAddBoxWithoutVector)
{
    std::string option = "box";
    std::string value  = "yes";
    setModuleFlag(option, value, getOption(), false);
    EXPECT_THROW(OutputAdapters adapters = getModules(),
                 InconsistentInputError);
}

} // namespace test

} // namespace gmx
