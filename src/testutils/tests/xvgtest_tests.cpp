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
 * Tests utilities for testing xvg files
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_testutils
 */
#include "gmxpre.h"

#include "testutils/xvgtest.h"

#include <vector>

#include <gtest/gtest.h>
#include <gtest/gtest-spi.h>

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/stringstream.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"

namespace
{

using gmx::test::checkXvgFile;
using gmx::test::XvgMatchSettings;

//! Input testing data - an inline xvg file.
const char * const input[] = {
    "0     2905.86    -410.199",
    "0.2     6656.67    -430.437",
    "0.4     5262.44    -409.399",
    "0.6     5994.69    -405.763",
    "0.8     5941.37    -408.337",
    "1     5869.87    -411.124"
};

TEST(XvgTests, CreateFile)
{
    {
        // Create new data
        gmx::test::TestReferenceData    data(gmx::test::erefdataUpdateAll);
        gmx::test::TestReferenceChecker checker(data.rootChecker());
        // Convert char array to a stream and add it to the checker
        gmx::StringInputStream          sis(input);
        checkXvgFile(&sis, &checker, XvgMatchSettings());
    }
    {
        // Now read it back
        gmx::test::TestReferenceData    data(gmx::test::erefdataCompare);
        gmx::test::TestReferenceChecker checker(data.rootChecker());
        // Convert char array to a stream and add it to the checker
        gmx::StringInputStream          sis(input);
        checkXvgFile(&sis, &checker, XvgMatchSettings());
    }
}

TEST(XvgTests, CheckMissing)
{
    {
        // Create new data
        gmx::test::TestReferenceData    data(gmx::test::erefdataUpdateAll);
        gmx::test::TestReferenceChecker checker(data.rootChecker());
        // Convert char array to a stream and add it to the checker
        gmx::StringInputStream          sis(input);
        checkXvgFile(&sis, &checker, XvgMatchSettings());
    }
    {
        const char * const              input[] = {
            "0     2905.86    -410.199",
            "0.2     6656.67    -430.437",
            "0.4     5262.44    -409.399"
        };
        // Now check with missing data
        gmx::test::TestReferenceData    data(gmx::test::erefdataCompare);
        gmx::test::TestReferenceChecker checker(data.rootChecker());
        gmx::StringInputStream          sis(input);
        EXPECT_NONFATAL_FAILURE(checkXvgFile(&sis, &checker, XvgMatchSettings()), "not used in test");
    }
}

TEST(XvgTests, CheckExtra)
{
    {
        // Create new data
        gmx::test::TestReferenceData    data(gmx::test::erefdataUpdateAll);
        gmx::test::TestReferenceChecker checker(data.rootChecker());
        // Convert char array to a stream and add it to the checker
        gmx::StringInputStream          sis(input);
        checkXvgFile(&sis, &checker, XvgMatchSettings());
    }
    {
        const char * const              input[] = {
            "0     2905.86    -410.199",
            "0.2     6656.67    -430.437",
            "0.4     5262.44    -409.399",
            "0.6     5994.69    -405.763",
            "0.8     5941.37    -408.337",
            "1     5869.87    -411.124",
            "1.2     5889.87    -413.124"
        };
        // Now check with missing data
        gmx::test::TestReferenceData    data(gmx::test::erefdataCompare);
        gmx::test::TestReferenceChecker checker(data.rootChecker());
        gmx::StringInputStream          sis(input);
        EXPECT_NONFATAL_FAILURE(checkXvgFile(&sis, &checker, XvgMatchSettings()), "Row6");
    }
}

TEST(XvgTests, ReadIncorrect)
{
    {
        // Create new data
        gmx::test::TestReferenceData    data(gmx::test::erefdataUpdateAll);
        gmx::test::TestReferenceChecker checker(data.rootChecker());
        // Convert char array to a stream and add it to the checker
        gmx::StringInputStream          sis(input);
        checkXvgFile(&sis, &checker, XvgMatchSettings());
    }
    {
        const char * const              input[] = {
            "0     2905.86    -410.199",
            "0.2     6656.67    -430.437",
            "0.4     5262.44    -409.399",
            "0.6     5994.69    -405.763",
            "0.8     5941.37    -408.337",
            "1     5869.87    -421.124"
        };
        // Now check with incorrect data
        gmx::test::TestReferenceData    data(gmx::test::erefdataCompare);
        gmx::test::TestReferenceChecker checker(data.rootChecker());
        gmx::StringInputStream          sis(input);
        EXPECT_NONFATAL_FAILURE(checkXvgFile(&sis, &checker, XvgMatchSettings()), "-411");
    }
}

} // namespace
