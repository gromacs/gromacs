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

namespace gmx
{

namespace test
{

TEST_F(FlagTest, CanSetSimpleFlag)
{
    std::string              option = "atoms";
    std::string              value  = "yes";
    setModuleFlag(option, value, &options_, TestEnums::efTestString);
    checkRequirementsOptions(&requirements_);
    EXPECT_EQ(requirements_.atoms, ChangeAtomsType::efUserYes);
}

TEST_F(FlagTest, CanAddNewBox)
{
    std::string              option = "newbox";
    std::string              value  = "3 3 3";
    setModuleFlag(option, value, &options_, TestEnums::efTestFloat);
    checkRequirementsOptions(&requirements_);
    EXPECT_EQ(requirements_.box, ChangeFrameInfoType::efUserYes);
}

TEST_F(FlagTest, CannotAddBoxWithoutVector)
{
    std::string option = "box";
    std::string value  = "yes";
    setModuleFlag(option, value, &options_, TestEnums::efTestString);
    EXPECT_THROW(checkRequirementsOptions(&requirements_),
                 InconsistentInputError);
}

TEST_F(FlagTest, SetsImplicitPrecisionChange)
{
    std::string option = "newprec";
    std::string value  = "5";
    setModuleFlag(option, value, &options_, TestEnums::efTestInt);
    checkRequirementsOptions(&requirements_);
    EXPECT_EQ(requirements_.precision, ChangeFrameInfoType::efUserYes);
}

TEST_F(FlagTest, SetsImplicitStartTimeChange)
{
    std::string option = "starttime";
    std::string value  = "20";
    setModuleFlag(option, value, &options_, TestEnums::efTestFloat);
    checkRequirementsOptions(&requirements_);
    EXPECT_EQ(requirements_.frameTime, ChangeFrameTimeType::efStartTime);
}

TEST_F(FlagTest, SetsImplicitTimeStepChange)
{
    std::string option = "timestep";
    std::string value  = "20";
    setModuleFlag(option, value, &options_, TestEnums::efTestFloat);
    checkRequirementsOptions(&requirements_);
    EXPECT_EQ(requirements_.frameTime, ChangeFrameTimeType::efTimeStep);
}

} // namespace test

} // namespace gmx
