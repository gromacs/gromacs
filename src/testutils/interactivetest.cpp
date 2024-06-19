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
 * Implements classes from interactivetest.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_testutils
 */
#include "gmxpre.h"

#include "testutils/interactivetest.h"

#include <memory>
#include <string>
#include <utility>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textstream.h"

#include "testutils/refdata.h"
#include "testutils/stringtest.h"

namespace gmx
{
namespace test
{

// These two classes cannot be in an unnamed namespace (easily), since
// then their use as members below would trigger warnings.
// But if anyone needs these outside this file, they can easily be moved to a
// separate header.

class MockTextInputStream : public TextInputStream
{
public:
    MOCK_METHOD1(readLine, bool(std::string*));
    MOCK_METHOD0(close, void());
};

class MockTextOutputStream : public TextOutputStream
{
public:
    MOCK_METHOD1(write, void(const char*));
    MOCK_METHOD0(close, void());
};

class InteractiveTestHelper::Impl
{
public:
    explicit Impl(TestReferenceChecker checker) :
        checker_(std::move(checker)), bLastNewline_(true), currentLine_(0), bHasOutput_(false)
    {
        using ::testing::_;
        using ::testing::Invoke;
        EXPECT_CALL(inputStream_, readLine(_)).WillRepeatedly(Invoke(this, &Impl::readInputLine));
        EXPECT_CALL(inputStream_, close()).Times(0);
        EXPECT_CALL(outputStream_, write(_)).WillRepeatedly(Invoke(this, &Impl::addOutput));
        EXPECT_CALL(outputStream_, close()).Times(0);
    }

    bool readInputLine(std::string* line)
    {
        checkOutput();
        line->clear();
        const bool bPresent = (currentLine_ < Index(inputLines_.size()));
        if (bPresent)
        {
            line->assign(inputLines_[currentLine_]);
            if (bLastNewline_ || currentLine_ + 1 < Index(inputLines_.size()))
            {
                line->append("\n");
            }
        }
        ++currentLine_;
        const std::string id = formatString("Input%d", static_cast<int>(currentLine_));
        StringTestBase::checkText(&checker_, *line, id.c_str());
        return bPresent;
    }
    void addOutput(const char* str)
    {
        bHasOutput_ = true;
        currentOutput_.append(str);
    }

    void checkOutput()
    {
        if (bHasOutput_)
        {
            const std::string id = formatString("Output%d", static_cast<int>(currentLine_));
            StringTestBase::checkText(&checker_, currentOutput_, id.c_str());
            bHasOutput_ = false;
        }
        currentOutput_.clear();
    }

    TestReferenceChecker        checker_;
    ArrayRef<const char* const> inputLines_;
    bool                        bLastNewline_;
    Index                       currentLine_;
    bool                        bHasOutput_;
    std::string                 currentOutput_;
    MockTextInputStream         inputStream_;
    MockTextOutputStream        outputStream_;
};

InteractiveTestHelper::InteractiveTestHelper(TestReferenceChecker checker) :
    impl_(new Impl(checker.checkCompound("InteractiveSession", "Interactive")))
{
}

InteractiveTestHelper::~InteractiveTestHelper() {}

void InteractiveTestHelper::setLastNewline(bool bInclude)
{
    impl_->bLastNewline_ = bInclude;
}

void InteractiveTestHelper::setInputLines(const ArrayRef<const char* const>& inputLines)
{
    impl_->inputLines_  = inputLines;
    impl_->currentLine_ = 0;
}

TextInputStream& InteractiveTestHelper::inputStream()
{
    return impl_->inputStream_;
}

TextOutputStream& InteractiveTestHelper::outputStream()
{
    return impl_->outputStream_;
}

void InteractiveTestHelper::checkSession()
{
    impl_->checkOutput();
    impl_->checker_.checkUnusedEntries();
}

} // namespace test
} // namespace gmx
