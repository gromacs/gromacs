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

#include <gmock/gmock.h>
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

//! Convenience name.
using Container = std::vector<std::string>;
//! Convenience type for callbacks.
using TestCallbackFunc = void(*)(TextReader &);

//! Helper struct.
struct TextReaderTestParams
{
    //! Input data.
    const Container  input;
    //! Callback to configure the reader with the behaviour being tested.
    TestCallbackFunc callback;
    //! Output to expect from the configured reader acting on the \c input.
    const Container  expectedOutput;
};

//! Test fixture.
class TextReaderTest : public ::testing::TestWithParam<TextReaderTestParams>
{
};

TEST_P(TextReaderTest, UsingDifferentConfigurations)
{
    const auto &params = GetParam();

    // Prepare the reader with the input lines.
    StringInputStream stream(params.input);
    TextReader        reader(&stream);
    // Configure the intended reading behaviour.
    params.callback(reader);
    // Read the input and store what is read.
    Container   readLines;
    std::string line;
    while (reader.readLine(&line))
    {
        readLines.push_back(line);
    }
    // Check the results
    EXPECT_THAT(readLines, ::testing::ElementsAreArray(params.expectedOutput));
}

//! Test input data. Some configurations will remove comments delimited by '#'.
const Container g_inputs =
{
    "",
    " \t ",
    "expected text",
    " expected text ",
    "expected text \t",
    " \t expected text",
    " \t expected text \t",
    "expected text#",
    "expected text\t #",
    "expected text# ",
    "expected text   # not expected ",
    "#",
    "\t #",
    "   # not expected ",
};

/*! \brief A case of expected output data that is produced from two
 * different configurations.
 *
 * Note that the implementation of StringInputStream joins the input
 * container with "\n", so the inputs are always changed before being
 * read. The name of this variable reflects that TextReader does not
 * change them during reading. */
const Container g_unchangedOutputs =
{
    "\n",
    " \t \n",
    "expected text\n",
    " expected text \n",
    "expected text \t\n",
    " \t expected text\n",
    " \t expected text \t\n",
    "expected text#\n",
    "expected text\t #\n",
    "expected text# \n",
    "expected text   # not expected \n",
    "#\n",
    "\t #\n",
    "   # not expected \n",
};
INSTANTIATE_TEST_CASE_P(ParsesLinesDifferently, TextReaderTest,
                            ::testing::Values(TextReaderTestParams {
                                                  g_inputs,
                                                  [](TextReader &r)
                                                  {
                                                      GMX_UNUSED_VALUE(r);
                                                  },
                                                  g_unchangedOutputs
                                              },
                                              TextReaderTestParams {
                                                  g_inputs,
                                                  [](TextReader &r)
                                                  {
                                                      r.setTrimLeadingWhiteSpace(true);
                                                  },
                                                  { "",
                                                    "",
                                                    "expected text\n",
                                                    "expected text \n",
                                                    "expected text \t\n",
                                                    "expected text\n",
                                                    "expected text \t\n",
                                                    "expected text#\n",
                                                    "expected text\t #\n",
                                                    "expected text# \n",
                                                    "expected text   # not expected \n",
                                                    "#\n",
                                                    "#\n",
                                                    "# not expected \n", }
                                              },
                                              TextReaderTestParams {
                                                  g_inputs,
                                                  [](TextReader &r)
                                                  {
                                                      r.setTrimTrailingWhiteSpace(true);
                                                  },
                                                  { "",
                                                    "",
                                                    "expected text",
                                                    " expected text",
                                                    "expected text",
                                                    " \t expected text",
                                                    " \t expected text",
                                                    "expected text#",
                                                    "expected text\t #",
                                                    "expected text#",
                                                    "expected text   # not expected",
                                                    "#",
                                                    "\t #",
                                                    "   # not expected", }
                                              },
                                              TextReaderTestParams {
                                                  g_inputs,
                                                  [](TextReader &r)
                                                  {
                                                      r.setTrimTrailingWhiteSpace(true);
                                                      r.setTrimLeadingWhiteSpace(true);
                                                  },
                                                  { "",
                                                    "",
                                                    "expected text",
                                                    "expected text",
                                                    "expected text",
                                                    "expected text",
                                                    "expected text",
                                                    "expected text#",
                                                    "expected text\t #",
                                                    "expected text#",
                                                    "expected text   # not expected",
                                                    "#",
                                                    "#",
                                                    "# not expected", }
                                              },
                                              TextReaderTestParams {
                                                  g_inputs,
                                                  [](TextReader &r)
                                                  {
                                                      r.setTrimTrailingComment(true, '#');
                                                  },
                                                  { "\n",
                                                    " \t \n",
                                                    "expected text\n",
                                                    " expected text \n",
                                                    "expected text \t\n",
                                                    " \t expected text\n",
                                                    " \t expected text \t\n",
                                                    "expected text",
                                                    "expected text\t ",
                                                    "expected text",
                                                    "expected text   ",
                                                    "",
                                                    "\t ",
                                                    "   ", }
                                              },
                                              TextReaderTestParams {
                                                  g_inputs,
                                                  [](TextReader &r)
                                                  {
                                                      r.setTrimTrailingComment(true, '#');
                                                      r.setTrimTrailingComment(false, 0);
                                                  },
                                                  g_unchangedOutputs
                                              },
                                              TextReaderTestParams {
                                                  g_inputs,
                                                  [](TextReader &r)
                                                  {
                                                      r.setTrimTrailingComment(true, '#');
                                                      r.setTrimTrailingWhiteSpace(true);
                                                  },
                                                  { "",
                                                    "",
                                                    "expected text",
                                                    " expected text",
                                                    "expected text",
                                                    " \t expected text",
                                                    " \t expected text",
                                                    "expected text",
                                                    "expected text",
                                                    "expected text",
                                                    "expected text",
                                                    "",
                                                    "",
                                                    "",   }
                                              }
                                              ));

} // namespace
} // namespace
} // namespace
