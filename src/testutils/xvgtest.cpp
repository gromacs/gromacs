/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2015- The GROMACS Authors
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
 * Implements routine to check the content of xvg files.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_testutils
 */
#include "gmxpre.h"

#include "testutils/xvgtest.h"

#include <cerrno>
#include <cstdlib>

#include <string>
#include <vector>

#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textstream.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/textblockmatchers.h"

namespace gmx
{
namespace test
{

namespace
{

class XvgMatcher : public ITextBlockMatcher
{
public:
    explicit XvgMatcher(const XvgMatchSettings& settings) : settings_(settings) {}

    void checkStream(TextInputStream* stream, TestReferenceChecker* checker) override
    {
        checkXvgFile(stream, checker, settings_);
    }

private:
    XvgMatchSettings settings_;
};

//! Helper function to identify which @ lines in xvg files should be tested.
bool isRelevantXvgCommand(const std::string& line)
{
    return contains(line, " title ") || contains(line, " subtitle ") || contains(line, " label ")
           || contains(line, "@TYPE ") || contains(line, " legend \"");
}

//! Helper function to check a single xvg value in a sequence.
void checkXvgDataPoint(TestReferenceChecker* checker, const std::string& value)
{
    checker->checkRealFromString(value, nullptr);
}

} // namespace

void checkXvgFile(TextInputStream* input, TestReferenceChecker* checker, const XvgMatchSettings& settings)
{
    TestReferenceChecker legendChecker(checker->checkCompound("XvgLegend", "Legend"));
    TestReferenceChecker dataChecker(checker->checkCompound("XvgData", "Data"));
    dataChecker.setDefaultTolerance(settings.tolerance);

    std::string legendText;
    int         dataRowCount = 0;
    std::string line;
    while (input->readLine(&line))
    {
        // Ignore comments, as they contain dynamic content, and very little of
        // that would be useful to test (and in particular, not with every
        // output file).
        if (startsWith(line, "#"))
        {
            continue;
        }
        // Ignore ampersand dataset separators (for now)
        // Later, when we need to test code that writes multiple
        // datasets, we might want to introduce that new concept
        // to this testing code.
        if (startsWith(line, "&"))
        {
            continue;
        }
        if (startsWith(line, "@"))
        {
            if (isRelevantXvgCommand(line))
            {
                legendText.append(stripString(line.substr(1)));
                legendText.append("\n");
            }
            continue;
        }
        if (!settings.testData)
        {
            break;
        }
        const std::vector<std::string> columns = splitString(line);
        const std::string              id      = formatString("Row%d", dataRowCount);
        dataChecker.checkSequence(columns.begin(), columns.end(), id.c_str(), &checkXvgDataPoint);
        ++dataRowCount;
    }
    dataChecker.checkUnusedEntries();
    legendChecker.checkTextBlock(legendText, "XvgLegend");
}

TextBlockMatcherPointer XvgMatch::createMatcher() const
{
    return TextBlockMatcherPointer(new XvgMatcher(settings_));
}

} // namespace test
} // namespace gmx
