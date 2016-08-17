/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
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
#include "gmxpre.h"

#include "gromacs/utility/logger.h"

#include <gtest/gtest.h>

#include "gromacs/utility/loggerbuilder.h"
#include "gromacs/utility/stringstream.h"

#include "testutils/stringtest.h"
#include "testutils/testfilemanager.h"

namespace
{

//! Test fixture for logging tests.
typedef gmx::test::StringTestBase LoggerTest;

TEST_F(LoggerTest, EmptyLoggerWorks)
{
    gmx::MDLogger logger;
    GMX_LOG(logger.info).appendText("foobar");
    GMX_LOG(logger.warning).appendText("foobar").asParagraph();
}

TEST_F(LoggerTest, LogsToStream)
{
    gmx::StringOutputStream stream;
    gmx::LoggerBuilder      builder;
    builder.addTargetStream(gmx::MDLogger::LogLevel::Info, &stream);
    gmx::LoggerOwner        owner  = builder.build();
    const gmx::MDLogger    &logger = owner.logger();
    GMX_LOG(logger.info).appendText("line");
    GMX_LOG(logger.warning).appendText("par").asParagraph();
    GMX_LOG(logger.info).appendText("line2");
    checkText(stream.toString(), "Output");
}

TEST_F(LoggerTest, LogsToFile)
{
    gmx::test::TestFileManager files;
    std::string                filename(files.getTemporaryFilePath("log.txt"));
    FILE                      *fp = fopen(filename.c_str(), "w");
    {
        gmx::LoggerBuilder      builder;
        builder.addTargetFile(gmx::MDLogger::LogLevel::Info, fp);
        gmx::LoggerOwner        owner  = builder.build();
        const gmx::MDLogger    &logger = owner.logger();
        GMX_LOG(logger.info).appendText("line");
        GMX_LOG(logger.warning).appendText("par").asParagraph();
        GMX_LOG(logger.info).appendText("line2");
    }
    fclose(fp);
    checkFileContents(filename, "Output");
}

TEST_F(LoggerTest, LevelFilteringWorks)
{
    gmx::StringOutputStream stream;
    gmx::LoggerBuilder      builder;
    builder.addTargetStream(gmx::MDLogger::LogLevel::Warning, &stream);
    gmx::LoggerOwner        owner  = builder.build();
    const gmx::MDLogger    &logger = owner.logger();
    GMX_LOG(logger.info).appendText("line");
    GMX_LOG(logger.warning).appendText("par").asParagraph();
    GMX_LOG(logger.info).appendText("line2");
    checkText(stream.toString(), "Output");
}

TEST_F(LoggerTest, LogsToMultipleStreams)
{
    gmx::StringOutputStream stream1;
    gmx::StringOutputStream stream2;
    gmx::LoggerBuilder      builder;
    builder.addTargetStream(gmx::MDLogger::LogLevel::Info, &stream1);
    builder.addTargetStream(gmx::MDLogger::LogLevel::Warning, &stream2);
    gmx::LoggerOwner        owner  = builder.build();
    const gmx::MDLogger    &logger = owner.logger();
    GMX_LOG(logger.info).appendText("line");
    GMX_LOG(logger.warning).appendText("par").asParagraph();
    GMX_LOG(logger.info).appendText("line2");
    checkText(stream1.toString(), "Output1");
    checkText(stream2.toString(), "Output2");
}

} // namespace
