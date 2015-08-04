/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015, by the GROMACS development team, led by
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
 * Tests file name option implementation dependent on gmx::FileNameOptionManager.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_options
 */
#include "gmxpre.h"

#include "gromacs/options/filenameoptionmanager.h"

#include <gtest/gtest.h>

#include "gromacs/fileio/filenm.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/options.h"
#include "gromacs/options/optionsassigner.h"
#include "gromacs/utility/exceptions.h"

#include "testutils/testasserts.h"
#include "testutils/testfileredirector.h"

namespace
{

using gmx::FileNameOption;

class FileNameOptionManagerTest : public ::testing::Test
{
    public:
        FileNameOptionManagerTest()
            : options_(NULL, NULL)
        {
            manager_.setInputRedirector(&redirector_);
            options_.addManager(&manager_);
        }

        void addExistingFile(const char *filename)
        {
            redirector_.addExistingFile(filename);
        }

        gmx::test::TestFileInputRedirector redirector_;
        gmx::FileNameOptionManager         manager_;
        gmx::Options                       options_;
};

/********************************************************************
 * Actual tests
 */

TEST_F(FileNameOptionManagerTest, AddsMissingExtension)
{
    std::string value;
    ASSERT_NO_THROW_GMX(options_.addOption(
                                FileNameOption("f").store(&value)
                                    .filetype(gmx::eftTrajectory).outputFile()));

    gmx::OptionsAssigner assigner(&options_);
    EXPECT_NO_THROW_GMX(assigner.start());
    EXPECT_NO_THROW_GMX(assigner.startOption("f"));
    EXPECT_NO_THROW_GMX(assigner.appendValue("testfile"));
    EXPECT_NO_THROW_GMX(assigner.finishOption());
    EXPECT_NO_THROW_GMX(assigner.finish());
    EXPECT_NO_THROW_GMX(options_.finish());

    EXPECT_EQ("testfile.xtc", value);
}

TEST_F(FileNameOptionManagerTest, AddsMissingCustomDefaultExtension)
{
    std::string value;
    ASSERT_NO_THROW_GMX(options_.addOption(
                                FileNameOption("f").store(&value)
                                    .filetype(gmx::eftTrajectory).outputFile()
                                    .defaultType(efPDB)));

    gmx::OptionsAssigner assigner(&options_);
    EXPECT_NO_THROW_GMX(assigner.start());
    EXPECT_NO_THROW_GMX(assigner.startOption("f"));
    EXPECT_NO_THROW_GMX(assigner.appendValue("testfile"));
    EXPECT_NO_THROW_GMX(assigner.finishOption());
    EXPECT_NO_THROW_GMX(assigner.finish());
    EXPECT_NO_THROW_GMX(options_.finish());

    EXPECT_EQ("testfile.pdb", value);
}

TEST_F(FileNameOptionManagerTest, GivesErrorOnMissingInputFile)
{
    std::string value;
    ASSERT_NO_THROW_GMX(options_.addOption(
                                FileNameOption("f").store(&value)
                                    .filetype(gmx::eftIndex).inputFile()));
    EXPECT_TRUE(value.empty());

    gmx::OptionsAssigner assigner(&options_);
    EXPECT_NO_THROW_GMX(assigner.start());
    EXPECT_NO_THROW_GMX(assigner.startOption("f"));
    EXPECT_THROW_GMX(assigner.appendValue("missing.ndx"), gmx::InvalidInputError);
    EXPECT_NO_THROW_GMX(assigner.finishOption());
    EXPECT_NO_THROW_GMX(assigner.finish());
    EXPECT_NO_THROW_GMX(options_.finish());

    EXPECT_TRUE(value.empty());
}

TEST_F(FileNameOptionManagerTest, GivesErrorOnMissingGenericInputFile)
{
    std::string value;
    ASSERT_NO_THROW_GMX(options_.addOption(
                                FileNameOption("f").store(&value)
                                    .filetype(gmx::eftTrajectory).inputFile()));
    EXPECT_TRUE(value.empty());

    gmx::OptionsAssigner assigner(&options_);
    EXPECT_NO_THROW_GMX(assigner.start());
    EXPECT_NO_THROW_GMX(assigner.startOption("f"));
    EXPECT_THROW_GMX(assigner.appendValue("missing.trr"), gmx::InvalidInputError);
    EXPECT_NO_THROW_GMX(assigner.finishOption());
    EXPECT_NO_THROW_GMX(assigner.finish());
    EXPECT_NO_THROW_GMX(options_.finish());

    EXPECT_TRUE(value.empty());
}

TEST_F(FileNameOptionManagerTest, GivesErrorOnMissingDefaultInputFile)
{
    std::string value;
    ASSERT_NO_THROW_GMX(options_.addOption(
                                FileNameOption("f").store(&value)
                                    .filetype(gmx::eftIndex).inputFile()
                                    .defaultBasename("missing")));

    gmx::OptionsAssigner assigner(&options_);
    EXPECT_NO_THROW_GMX(assigner.start());
    EXPECT_NO_THROW_GMX(assigner.startOption("f"));
    EXPECT_NO_THROW_GMX(assigner.finishOption());
    EXPECT_NO_THROW_GMX(assigner.finish());
    EXPECT_THROW_GMX(options_.finish(), gmx::InvalidInputError);
}

TEST_F(FileNameOptionManagerTest, GivesErrorOnMissingRequiredInputFile)
{
    std::string value;
    ASSERT_NO_THROW_GMX(options_.addOption(
                                FileNameOption("f").store(&value).required()
                                    .filetype(gmx::eftIndex).inputFile()
                                    .defaultBasename("missing")));
    EXPECT_EQ("missing.ndx", value);

    gmx::OptionsAssigner assigner(&options_);
    EXPECT_NO_THROW_GMX(assigner.start());
    EXPECT_NO_THROW_GMX(assigner.finish());
    EXPECT_THROW_GMX(options_.finish(), gmx::InvalidInputError);
}

TEST_F(FileNameOptionManagerTest, AcceptsMissingInputFileIfSpecified)
{
    std::string value;
    ASSERT_NO_THROW_GMX(options_.addOption(
                                FileNameOption("f").store(&value)
                                    .filetype(gmx::eftIndex).inputFile()
                                    .allowMissing()));
    EXPECT_TRUE(value.empty());

    gmx::OptionsAssigner assigner(&options_);
    EXPECT_NO_THROW_GMX(assigner.start());
    EXPECT_NO_THROW_GMX(assigner.startOption("f"));
    EXPECT_NO_THROW_GMX(assigner.appendValue("missing.ndx"));
    EXPECT_NO_THROW_GMX(assigner.finishOption());
    EXPECT_NO_THROW_GMX(assigner.finish());
    EXPECT_NO_THROW_GMX(options_.finish());

    EXPECT_EQ("missing.ndx", value);
}

TEST_F(FileNameOptionManagerTest, AcceptsMissingDefaultInputFileIfSpecified)
{
    std::string value;
    ASSERT_NO_THROW_GMX(options_.addOption(
                                FileNameOption("f").store(&value)
                                    .filetype(gmx::eftIndex).inputFile()
                                    .defaultBasename("missing")
                                    .allowMissing()));

    gmx::OptionsAssigner assigner(&options_);
    EXPECT_NO_THROW_GMX(assigner.start());
    EXPECT_NO_THROW_GMX(assigner.startOption("f"));
    EXPECT_NO_THROW_GMX(assigner.finishOption());
    EXPECT_NO_THROW_GMX(assigner.finish());
    EXPECT_NO_THROW_GMX(options_.finish());

    EXPECT_EQ("missing.ndx", value);
}

TEST_F(FileNameOptionManagerTest, AcceptsMissingRequiredInputFileIfSpecified)
{
    std::string value;
    ASSERT_NO_THROW_GMX(options_.addOption(
                                FileNameOption("f").store(&value).required()
                                    .filetype(gmx::eftIndex).inputFile()
                                    .defaultBasename("missing")
                                    .allowMissing()));
    EXPECT_EQ("missing.ndx", value);

    gmx::OptionsAssigner assigner(&options_);
    EXPECT_NO_THROW_GMX(assigner.start());
    EXPECT_NO_THROW_GMX(assigner.finish());
    EXPECT_NO_THROW_GMX(options_.finish());

    EXPECT_EQ("missing.ndx", value);
}

TEST_F(FileNameOptionManagerTest, AddsMissingExtensionBasedOnExistingFile)
{
    addExistingFile("testfile.trr");

    std::string value;
    ASSERT_NO_THROW_GMX(options_.addOption(
                                FileNameOption("f").store(&value)
                                    .filetype(gmx::eftTrajectory).inputFile()));

    gmx::OptionsAssigner assigner(&options_);
    EXPECT_NO_THROW_GMX(assigner.start());
    EXPECT_NO_THROW_GMX(assigner.startOption("f"));
    EXPECT_NO_THROW_GMX(assigner.appendValue("testfile"));
    EXPECT_NO_THROW_GMX(assigner.finishOption());
    EXPECT_NO_THROW_GMX(assigner.finish());
    EXPECT_NO_THROW_GMX(options_.finish());

    EXPECT_EQ("testfile.trr", value);
}

TEST_F(FileNameOptionManagerTest,
       AddsMissingExtensionForRequiredDefaultNameBasedOnExistingFile)
{
    addExistingFile("testfile.trr");

    std::string value;
    ASSERT_NO_THROW_GMX(options_.addOption(
                                FileNameOption("f").store(&value).required()
                                    .filetype(gmx::eftTrajectory).inputFile()
                                    .defaultBasename("testfile")));
    EXPECT_EQ("testfile.xtc", value);

    gmx::OptionsAssigner assigner(&options_);
    EXPECT_NO_THROW_GMX(assigner.start());
    EXPECT_NO_THROW_GMX(assigner.startOption("f"));
    EXPECT_NO_THROW_GMX(assigner.finishOption());
    EXPECT_NO_THROW_GMX(assigner.finish());
    EXPECT_NO_THROW_GMX(options_.finish());

    EXPECT_EQ("testfile.trr", value);
}

TEST_F(FileNameOptionManagerTest,
       AddsMissingExtensionForOptionalDefaultNameBasedOnExistingFile)
{
    addExistingFile("testfile.trr");

    std::string value;
    ASSERT_NO_THROW_GMX(options_.addOption(
                                FileNameOption("f").store(&value)
                                    .filetype(gmx::eftTrajectory).inputFile()
                                    .defaultBasename("testfile")));

    gmx::OptionsAssigner assigner(&options_);
    EXPECT_NO_THROW_GMX(assigner.start());
    EXPECT_NO_THROW_GMX(assigner.startOption("f"));
    EXPECT_NO_THROW_GMX(assigner.finishOption());
    EXPECT_NO_THROW_GMX(assigner.finish());
    EXPECT_NO_THROW_GMX(options_.finish());

    EXPECT_EQ("testfile.trr", value);
}

TEST_F(FileNameOptionManagerTest,
       AddsMissingExtensionForRequiredFromDefaultNameOptionBasedOnExistingFile)
{
    addExistingFile("testfile.trr");

    std::string value;
    ASSERT_NO_THROW_GMX(options_.addOption(
                                FileNameOption("f").store(&value).required()
                                    .filetype(gmx::eftTrajectory).inputFile()
                                    .defaultBasename("foo")));
    ASSERT_NO_THROW_GMX(manager_.addDefaultFileNameOption(&options_, "deffnm"));
    EXPECT_EQ("foo.xtc", value);

    gmx::OptionsAssigner assigner(&options_);
    EXPECT_NO_THROW_GMX(assigner.start());
    EXPECT_NO_THROW_GMX(assigner.startOption("deffnm"));
    EXPECT_NO_THROW_GMX(assigner.appendValue("testfile"));
    EXPECT_NO_THROW_GMX(assigner.finishOption());
    EXPECT_NO_THROW_GMX(assigner.finish());
    EXPECT_NO_THROW_GMX(options_.finish());

    EXPECT_EQ("testfile.trr", value);
}

TEST_F(FileNameOptionManagerTest,
       AddsMissingExtensionForOptionalFromDefaultNameOptionBasedOnExistingFile)
{
    addExistingFile("testfile.trr");

    std::string value;
    ASSERT_NO_THROW_GMX(options_.addOption(
                                FileNameOption("f").store(&value)
                                    .filetype(gmx::eftTrajectory).inputFile()
                                    .defaultBasename("foo")));
    ASSERT_NO_THROW_GMX(manager_.addDefaultFileNameOption(&options_, "deffnm"));

    gmx::OptionsAssigner assigner(&options_);
    EXPECT_NO_THROW_GMX(assigner.start());
    EXPECT_NO_THROW_GMX(assigner.startOption("deffnm"));
    EXPECT_NO_THROW_GMX(assigner.appendValue("testfile"));
    EXPECT_NO_THROW_GMX(assigner.finishOption());
    EXPECT_NO_THROW_GMX(assigner.startOption("f"));
    EXPECT_NO_THROW_GMX(assigner.finishOption());
    EXPECT_NO_THROW_GMX(assigner.finish());
    EXPECT_NO_THROW_GMX(options_.finish());

    EXPECT_EQ("testfile.trr", value);
}

TEST_F(FileNameOptionManagerTest, DefaultNameOptionWorksWithoutInputChecking)
{
    std::string value;
    ASSERT_NO_THROW_GMX(manager_.disableInputOptionChecking(true));
    ASSERT_NO_THROW_GMX(options_.addOption(
                                FileNameOption("f").store(&value).required()
                                    .filetype(gmx::eftIndex).inputFile()
                                    .defaultBasename("default")
                                    .allowMissing()));
    ASSERT_NO_THROW_GMX(manager_.addDefaultFileNameOption(&options_, "deffnm"));
    EXPECT_EQ("default.ndx", value);

    gmx::OptionsAssigner assigner(&options_);
    EXPECT_NO_THROW_GMX(assigner.start());
    EXPECT_NO_THROW_GMX(assigner.startOption("deffnm"));
    EXPECT_NO_THROW_GMX(assigner.appendValue("missing"));
    EXPECT_NO_THROW_GMX(assigner.finishOption());
    EXPECT_NO_THROW_GMX(assigner.finish());
    EXPECT_NO_THROW_GMX(options_.finish());

    EXPECT_EQ("missing.ndx", value);
}

} // namespace
