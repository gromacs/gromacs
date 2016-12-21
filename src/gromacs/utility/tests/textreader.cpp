/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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
 * Tests for gmx::TextReader.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_utility
 */
#include "gmxpre.h"

#include "gromacs/utility/textreader.h"

#include <functional>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/stringstream.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/refdata.h"
#include "testutils/stringtest.h"

namespace gmx
{
namespace test
{
namespace
{

//! Test fixture.
class TextReaderTest : public StringTestBase
{
    public:
        //! Convenience type for callbacks.
        using testCallbackFunc = void(*)(TextReader &);
        /*! \brief Runs the test.
         *
         * Reads lines of \c input via a string-based TextReader,
         * configures the behaviour of the reader with the \c
         * callback, then reads the input, and checks each resulting
         * line against reference data.
         */
        void runTest(const std::vector<std::string>   &input,
                     testCallbackFunc                  callback,
                     const char                       *id)
        {
            StringInputStream   stream(input);
            TextReader          reader(&stream);
            // Configure the reader with the intended behaviour
            callback(reader);
            // Read the input and store the output
            std::vector<std::string> output;
            std::string              line;
            while (reader.readLine(&line))
            {
                output.push_back(line);
            }
            // Check the output
            TestReferenceChecker compound(checker().checkCompound("Strings", id));
            int                  i = 0;
            for (auto &s : output)
            {
                StringTestBase::checkText(&compound, s, formatString("%d", i).c_str());
                ++i;
            }
        }
};

TEST_F(TextReaderTest, CanTrimWhiteSpace)
{
    std::vector<std::string> input =
    {
        "",
        " \t ",
        "expected text",
        " expected text ",
        "expected text \t",
        " \t expected text",
        " \t expected text \t",
    };
    runTest(input, [](TextReader &r)
            {
                GMX_UNUSED_VALUE(r);
            }, "NothingTrimmed");
    runTest(input, [](TextReader &r)
            {
                r.setTrimLeadingWhiteSpace(true);
            }, "LeadingWhiteSpaceTrimmed");
    runTest(input, [](TextReader &r)
            {
                r.setTrimTrailingWhiteSpace(true);
            }, "TrailingWhiteSpaceTrimmed");
    runTest(input, [](TextReader &r)
            {
                r.setTrimTrailingWhiteSpace(true);
                r.setTrimLeadingWhiteSpace(true);
            }, "LeadingAndTrailingWhiteSpaceTrimmed");
}

TEST_F(TextReaderTest, CanTrimComments)
{
    std::vector<std::string> input =
    {
        "",
        " \t ",
        "expected text",
        "expected text \t ",
        "expected text#",
        "expected text\t #",
        "expected text   # not expected",
        "#",
        "\t #",
        "   # not expected",
    };
    runTest(input, [](TextReader &r)
            {
                r.setCommentChar('#');
                r.setTrimTrailingComment(true);
            }, "TrailingCommentTrimmed");
    runTest(input, [](TextReader &r)
            {
                r.setCommentChar('#');
                r.setTrimTrailingComment(true);
                r.setTrimTrailingWhiteSpace(true);
            }, "TrailingCommentAndWhiteSpaceTrimmed");
}

} // namespace
} // namespace
} // namespace
