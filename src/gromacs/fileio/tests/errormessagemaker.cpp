/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015, by the GROMACS development team, led by
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
 * Tests for helper class for making error messages when reading line-based input
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_fileio
 */
#include "gmxpre.h"

#include "gromacs/fileio/errormessagemaker.h"

#include <string>

#include <gtest/gtest.h>

#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

namespace
{

TEST(ErrorMessageMakerTest, ConstructsUsefulMessages)
{
    gmx::ErrorMessageMaker errorMessage;
    EXPECT_STREQ("Before reading a line", errorMessage.make("Before reading a line").what());

    std::string currentLine("here's the line we read");
    errorMessage.setCurrentLine(&currentLine);
    EXPECT_STREQ("After reading a line\n on line 1 which was\n 'here's the line we read'",
                 errorMessage.make("After reading a line").what());

    currentLine = "here's the line we read\n"; // now on line 2, and including a newline to strip
    errorMessage.setCurrentLine(&currentLine);
    EXPECT_STREQ("After reading a line\n on line 2 which was\n 'here's the line we read'",
                 errorMessage.make("After reading a line").what());
}

} // namespace
