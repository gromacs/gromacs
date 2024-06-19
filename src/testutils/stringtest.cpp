/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2012- The GROMACS Authors
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
 * Implements gmx::test::StringTestBase.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_testutils
 */
#include "gmxpre.h"

#include "testutils/stringtest.h"

#include <cstdio>

#include <filesystem>
#include <memory>
#include <string>

#include <gtest/gtest.h>

#include "gromacs/options/basicoptions.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/utility/textreader.h"

#include "testutils/refdata.h"
#include "testutils/testoptions.h"

namespace gmx
{
namespace test
{

namespace
{
//! Stores the -stdout flag value to print out values instead of checking them.
bool g_bWriteToStdOut = false;
} // namespace

// TODO: Only add this option to those test binaries that actually need it
// (depending on the linker, it may or may not appear right now),
// or replace by a generic mechanism in TestReferenceData.
//! \cond
GMX_TEST_OPTIONS(StringTestOptions, options)
{
    options->addOption(
            BooleanOption("stdout").store(&g_bWriteToStdOut).description("Print the test string to stdout instead of checking against reference data"));
}
//! \endcond

/********************************************************************
 * StringTestBase::Impl
 */

class StringTestBase::Impl
{
public:
    TestReferenceData    data_;
    TestReferenceChecker checker_;
};

/********************************************************************
 * StringTestBase
 */

// static
void StringTestBase::checkText(TestReferenceChecker* checker, const std::string& text, const char* id)
{
    if (g_bWriteToStdOut)
    {
        printf("%s:\n", id);
        printf("%s[END]\n", text.c_str());
    }
    else
    {
        checker->checkTextBlock(text, id);
    }
}

StringTestBase::StringTestBase() : impl_(new Impl) {}

StringTestBase::~StringTestBase() {}

TestReferenceChecker& StringTestBase::checker()
{
    if (!impl_->checker_)
    {
        impl_->checker_ = impl_->data_.rootChecker();
    }
    return impl_->checker_;
}

void StringTestBase::checkText(const std::string& text, const char* id)
{
    checkText(&checker(), text, id);
}

void StringTestBase::checkFileContents(const std::filesystem::path& filename, const char* id)
{
    const std::string text = TextReader::readFileToString(filename);
    checkText(text, id);
}

void StringTestBase::testFilesEqual(const std::filesystem::path& refFilename,
                                    const std::filesystem::path& testFilename)
{
    const std::string expectedContents = TextReader::readFileToString(refFilename);
    const std::string contents         = TextReader::readFileToString(testFilename);
    if (g_bWriteToStdOut)
    {
        printf("%s[END]\n", contents.c_str());
    }
    EXPECT_EQ(expectedContents, contents);
}

} // namespace test
} // namespace gmx
