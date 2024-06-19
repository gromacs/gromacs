/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2017- The GROMACS Authors
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
 * Implements routine to check the content of conf files.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_testutils
 */
#include "gmxpre.h"

#include "testutils/conftest.h"

#include <cstdio>
#include <cstdlib>

#include <string>

#include <gtest/gtest.h>

#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textreader.h"
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

class ConfMatcher : public ITextBlockMatcher
{
public:
    explicit ConfMatcher(const ConfMatchSettings& settings) : settings_(settings) {}

    void checkStream(TextInputStream* stream, TestReferenceChecker* checker) override
    {
        checkConfFile(stream, checker, settings_);
    }

private:
    ConfMatchSettings settings_;
};

} // namespace

void checkConfFile(TextInputStream* input, TestReferenceChecker* checker, const ConfMatchSettings& settings)
{

    TestReferenceChecker groChecker(checker->checkCompound("GroFile", "Header"));
    // Check the first two lines of the output file
    std::string line;
    EXPECT_TRUE(input->readLine(&line));
    line = stripSuffixIfPresent(line, "\n");
    groChecker.checkString(line, "Title");
    EXPECT_TRUE(input->readLine(&line));
    line = stripSuffixIfPresent(line, "\n");
    groChecker.checkInteger(std::atoi(line.c_str()), "Number of atoms");
    // Check the full configuration only if required
    if (settings.matchFullConfiguration)
    {
        TextReader reader(input);
        checker->checkTextBlock(reader.readAll(), "Configuration");
    }
}

TextBlockMatcherPointer ConfMatch::createMatcher() const
{
    return TextBlockMatcherPointer(new ConfMatcher(settings_));
}

} // namespace test
} // namespace gmx
