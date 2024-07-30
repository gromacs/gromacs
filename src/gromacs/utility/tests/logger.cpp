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
#include "gmxpre.h"

#include "gromacs/utility/logger.h"

#include <cstdio>

#include <filesystem>
#include <string>

#include <gtest/gtest.h>

#include "gromacs/utility/loggerbuilder.h"
#include "gromacs/utility/stringstream.h"

#include "testutils/stringtest.h"
#include "testutils/testfilemanager.h"

namespace gmx
{
namespace test
{
namespace
{

//! Test fixture for logging tests.
typedef gmx::test::StringTestBase LoggerTest;

TEST_F(LoggerTest, EmptyLoggerWorks)
{
    gmx::MDLogger logger;
    GMX_LOG(logger.info).appendText("foobar");
    GMX_LOG(logger.warning).appendText("foobar").asParagraph();
    GMX_LOG(logger.debug).appendText("foobaz");
    GMX_LOG(logger.error).appendText("baz");
    GMX_LOG(logger.verboseDebug).appendText("verbose");
}

TEST_F(LoggerTest, LogsToStream)
{
    gmx::StringOutputStream stream;
    gmx::LoggerBuilder      builder;
    builder.addTargetStream(gmx::MDLogger::LogLevel::VerboseDebug, &stream);
    gmx::LoggerOwner     owner  = builder.build();
    const gmx::MDLogger& logger = owner.logger();
    GMX_LOG(logger.info).appendText("line");
    GMX_LOG(logger.warning).appendText("par").asParagraph();
    GMX_LOG(logger.info).appendText("line2");
    GMX_LOG(logger.error).appendTextFormatted("%s", "formatted");
    GMX_LOG(logger.debug).appendText("debugline");
    GMX_LOG(logger.verboseDebug).appendText("verbose");
    checkText(stream.toString(), "Output");
}

TEST_F(LoggerTest, LogsToFile)
{
    gmx::test::TestFileManager files;
    std::string                filename(files.getTemporaryFilePath("log.txt").string());
    FILE*                      fp = fopen(filename.c_str(), "w");
    {
        gmx::LoggerBuilder builder;
        builder.addTargetFile(gmx::MDLogger::LogLevel::VerboseDebug, fp);
        gmx::LoggerOwner     owner  = builder.build();
        const gmx::MDLogger& logger = owner.logger();
        GMX_LOG(logger.info).appendText("line");
        GMX_LOG(logger.warning).appendText("par").asParagraph();
        GMX_LOG(logger.info).appendText("line2");
        GMX_LOG(logger.error).appendTextFormatted("%s", "formatted");
        GMX_LOG(logger.debug).appendText("debugline");
        GMX_LOG(logger.verboseDebug).appendText("verbose");
    }
    fclose(fp);
    checkFileContents(filename, "Output");
}

TEST_F(LoggerTest, LevelFilteringWorks)
{
    gmx::StringOutputStream stream;
    gmx::LoggerBuilder      builder;
    builder.addTargetStream(gmx::MDLogger::LogLevel::Warning, &stream);
    gmx::LoggerOwner     owner  = builder.build();
    const gmx::MDLogger& logger = owner.logger();
    GMX_LOG(logger.info).appendText("line");
    GMX_LOG(logger.warning).appendText("par").asParagraph();
    GMX_LOG(logger.info).appendText("line2");
    GMX_LOG(logger.error).appendTextFormatted("%s", "formatted");
    GMX_LOG(logger.debug).appendText("debugline");
    GMX_LOG(logger.verboseDebug).appendText("verbose");
    checkText(stream.toString(), "Output");
}

TEST_F(LoggerTest, LogsToMultipleStreams)
{
    gmx::StringOutputStream stream1;
    gmx::StringOutputStream stream2;
    gmx::StringOutputStream stream3;
    gmx::StringOutputStream stream4;
    gmx::StringOutputStream stream5;
    gmx::LoggerBuilder      builder;
    builder.addTargetStream(gmx::MDLogger::LogLevel::Info, &stream1);
    builder.addTargetStream(gmx::MDLogger::LogLevel::Warning, &stream2);
    builder.addTargetStream(gmx::MDLogger::LogLevel::Error, &stream3);
    builder.addTargetStream(gmx::MDLogger::LogLevel::Debug, &stream4);
    builder.addTargetStream(gmx::MDLogger::LogLevel::VerboseDebug, &stream5);
    gmx::LoggerOwner     owner  = builder.build();
    const gmx::MDLogger& logger = owner.logger();
    GMX_LOG(logger.info).appendText("line");
    GMX_LOG(logger.warning).appendText("par").asParagraph();
    GMX_LOG(logger.info).appendText("line2");
    GMX_LOG(logger.error).appendTextFormatted("%s", "formatted");
    GMX_LOG(logger.debug).appendText("debugline");
    GMX_LOG(logger.verboseDebug).appendText("verbose");

    checkText(stream1.toString(), "Output1");
    checkText(stream2.toString(), "Output2");
    checkText(stream3.toString(), "Output3");
    checkText(stream4.toString(), "Output4");
    checkText(stream5.toString(), "Output5");
}

TEST_F(LoggerTest, LogsToMultipleFiles)
{
    gmx::test::TestFileManager files;
    std::string                filename1(files.getTemporaryFilePath("log.txt").string());
    std::string                filename2(files.getTemporaryFilePath("warn.txt").string());
    std::string                filename3(files.getTemporaryFilePath("error.txt").string());
    std::string                filename4(files.getTemporaryFilePath("debug.txt").string());
    std::string                filename5(files.getTemporaryFilePath("verboseDebug.txt").string());
    FILE*                      fp1 = fopen(filename1.c_str(), "w");
    FILE*                      fp2 = fopen(filename2.c_str(), "w");
    FILE*                      fp3 = fopen(filename3.c_str(), "w");
    FILE*                      fp4 = fopen(filename4.c_str(), "w");
    FILE*                      fp5 = fopen(filename5.c_str(), "w");
    {
        gmx::LoggerBuilder builder;
        builder.addTargetFile(gmx::MDLogger::LogLevel::Info, fp1);
        builder.addTargetFile(gmx::MDLogger::LogLevel::Warning, fp2);
        builder.addTargetFile(gmx::MDLogger::LogLevel::Error, fp3);
        builder.addTargetFile(gmx::MDLogger::LogLevel::Debug, fp4);
        builder.addTargetFile(gmx::MDLogger::LogLevel::VerboseDebug, fp5);
        gmx::LoggerOwner     owner  = builder.build();
        const gmx::MDLogger& logger = owner.logger();
        GMX_LOG(logger.info).appendText("line");
        GMX_LOG(logger.warning).appendText("par").asParagraph();
        GMX_LOG(logger.info).appendText("line2");
        GMX_LOG(logger.error).appendTextFormatted("%s", "formatted");
        GMX_LOG(logger.debug).appendText("debugline");
        GMX_LOG(logger.verboseDebug).appendText("verbose");
    }
    fclose(fp1);
    fclose(fp2);
    fclose(fp3);
    fclose(fp4);
    fclose(fp5);
    checkFileContents(filename1, "Output1");
    checkFileContents(filename2, "Output2");
    checkFileContents(filename3, "Output3");
    checkFileContents(filename4, "Output4");
    checkFileContents(filename5, "Output5");
}

TEST_F(LoggerTest, LogsToStreamAndFile)
{
    gmx::test::TestFileManager files;
    gmx::StringOutputStream    stream;
    std::string                filename(files.getTemporaryFilePath("verboseDebug.txt").string());
    FILE*                      fp = fopen(filename.c_str(), "w");
    {
        gmx::LoggerBuilder builder;
        builder.addTargetFile(gmx::MDLogger::LogLevel::VerboseDebug, fp);
        builder.addTargetStream(gmx::MDLogger::LogLevel::VerboseDebug, &stream);
        gmx::LoggerOwner     owner  = builder.build();
        const gmx::MDLogger& logger = owner.logger();
        GMX_LOG(logger.info).appendText("line");
        GMX_LOG(logger.warning).appendText("par").asParagraph();
        GMX_LOG(logger.info).appendText("line2");
        GMX_LOG(logger.error).appendTextFormatted("%s", "formatted");
        GMX_LOG(logger.debug).appendText("debugline");
        GMX_LOG(logger.verboseDebug).appendText("verbose");
    }
    fclose(fp);
    checkText(stream.toString(), "OutputStream");
    checkFileContents(filename, "OutputFile");
}

} // namespace
} // namespace test
} // namespace gmx
