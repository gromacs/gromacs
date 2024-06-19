/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2016- The GROMACS Authors
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
 * Implements gmx::test::LoggerTestHelper.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_testutils
 */
#include "gmxpre.h"

#include "testutils/loggertest.h"

#include <memory>
#include <string>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/logger.h"

namespace gmx
{
namespace test
{

using ::testing::NiceMock;

namespace
{
class MockLogTarget : public ILogTarget
{
public:
    MOCK_METHOD1(writeEntry, void(const LogEntry&));
};
} // namespace

/********************************************************************
 * LoggerTestHelper::Impl
 */

class LoggerTestHelper::Impl
{
public:
    Impl()
    {
        // TODO: Add support for -stdout for echoing the log to stdout.
        logger_.warning = LogLevelHelper(&getTarget(MDLogger::LogLevel::Warning));
        logger_.info    = LogLevelHelper(&getTarget(MDLogger::LogLevel::Info));
    }

    NiceMock<MockLogTarget>& getTarget(MDLogger::LogLevel level)
    {
        return targets_[static_cast<int>(level)];
    }

    NiceMock<MockLogTarget> targets_[MDLogger::LogLevelCount];
    MDLogger                logger_;
};

/********************************************************************
 * LoggerTestHelper
 */

LoggerTestHelper::LoggerTestHelper() : impl_(new Impl) {}

LoggerTestHelper::~LoggerTestHelper() {}

const MDLogger& LoggerTestHelper::logger()
{
    return impl_->logger_;
}

void LoggerTestHelper::expectEntryMatchingRegex(gmx::MDLogger::LogLevel level, const char* re)
{
    using ::testing::ContainsRegex;
    using ::testing::Field;
    auto& target = impl_->getTarget(level);
    EXPECT_CALL(target, writeEntry(Field(&LogEntry::text, ContainsRegex(re))));
}

void LoggerTestHelper::expectNoEntries(gmx::MDLogger::LogLevel level)
{
    auto& target = impl_->getTarget(level);
    EXPECT_CALL(target, writeEntry(testing::_)).Times(0);
}

} // namespace test
} // namespace gmx
