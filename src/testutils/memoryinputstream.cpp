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
 * Implements classes from memoryinputstream.h.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_testutils
 */
#include "gmxpre.h"

#include "memoryinputstream.h"

#include <string>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textstream.h"

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
        MOCK_METHOD1(readLine, bool(std::string *));
        MOCK_METHOD0(close, void());
};

class MockTextOutputStream : public TextOutputStream
{
    public:
        MOCK_METHOD1(write, void(const char *));
        MOCK_METHOD0(close, void());
};

class MemoryInputStream::Impl
{
    public:
        Impl() : currentLine_(0)
        {
            using ::testing::_;
            using ::testing::Invoke;
            EXPECT_CALL(inputStream_, readLine(_))
                .WillRepeatedly(Invoke(this, &Impl::readInputLine));
            EXPECT_CALL(inputStream_, close()).Times(0);
        }

        bool readInputLine(std::string *line)
        {
            line->clear();
            const bool bPresent = (currentLine_ < inputLines_.size());
            if (bPresent)
            {
                line->assign(inputLines_[currentLine_]);
                line->append("\n");
            }
            ++currentLine_;
            const std::string id = formatString("Input%d", static_cast<int>(currentLine_));
            return bPresent;
        }

        ConstArrayRef<const char *> inputLines_;
        size_t                      currentLine_;
        MockTextInputStream         inputStream_;
};

MemoryInputStream::MemoryInputStream() : impl_(new Impl())
{
}

MemoryInputStream::~MemoryInputStream()
{
}

void MemoryInputStream::setInputLines(const ConstArrayRef<const char *> &inputLines)
{
    impl_->inputLines_  = inputLines;
    impl_->currentLine_ = 0;
}

TextInputStream &MemoryInputStream::inputStream()
{
    return impl_->inputStream_;
}

} // namespace test
} // namespace gmx
