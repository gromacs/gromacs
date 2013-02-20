/*
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 *
 * For more info, check our website at http://www.gromacs.org
 */
/*! \internal \file
 * \brief
 * Tests file name option implementation.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_options
 */
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/options/filenameoption.h"
#include "gromacs/options/options.h"
#include "gromacs/options/optionsassigner.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/file.h"

#include "testutils/testfilemanager.h"

namespace
{

using gmx::FileNameOption;
using gmx::test::TestFileManager;

TEST(FileNameOptionTest, AddsMissingExtension)
{
    gmx::Options           options(NULL, NULL);
    std::string            value;
    ASSERT_NO_THROW(options.addOption(
                            FileNameOption("f").store(&value)
                                .filetype(gmx::eftTrajectory).outputFile()));

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW(assigner.start());
    EXPECT_NO_THROW(assigner.startOption("f"));
    EXPECT_NO_THROW(assigner.appendValue("testfile"));
    EXPECT_NO_THROW(assigner.finishOption());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(options.finish());

    EXPECT_EQ("testfile.xtc", value);
}

TEST(FileNameOptionTest, HandlesRequiredDefaultValueWithoutExtension)
{
    gmx::Options           options(NULL, NULL);
    std::string            value;
    ASSERT_NO_THROW(options.addOption(
                            FileNameOption("f").store(&value).required()
                                .filetype(gmx::eftGenericData).outputFile()
                                .defaultBasename("testfile")));
    EXPECT_EQ("testfile.dat", value);

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW(assigner.start());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(options.finish());

    EXPECT_EQ("testfile.dat", value);
}

TEST(FileNameOptionTest, HandlesOptionalDefaultValueWithoutExtension)
{
    gmx::Options           options(NULL, NULL);
    std::string            value;
    ASSERT_NO_THROW(options.addOption(
                            FileNameOption("f").store(&value)
                                .filetype(gmx::eftIndex).outputFile()
                                .defaultBasename("testfile")));
    EXPECT_TRUE(value.empty());

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW(assigner.start());
    EXPECT_NO_THROW(assigner.startOption("f"));
    EXPECT_NO_THROW(assigner.finishOption());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(options.finish());

    EXPECT_EQ("testfile.ndx", value);
}

TEST(FileNameOptionTest, AddsMissingExtensionBasedOnExistingFile)
{
    gmx::Options           options(NULL, NULL);
    std::string            value;
    ASSERT_NO_THROW(options.addOption(
                            FileNameOption("f").store(&value)
                                .filetype(gmx::eftTrajectory).inputFile()));
    TestFileManager      tempFiles;
    std::string          filename(tempFiles.getTemporaryFilePath(".trr"));
    gmx::File::writeFileFromString(filename, "Dummy trajectory file");
    std::string          inputValue(filename.substr(0, filename.length() - 4));

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW(assigner.start());
    EXPECT_NO_THROW(assigner.startOption("f"));
    EXPECT_NO_THROW(assigner.appendValue(inputValue));
    EXPECT_NO_THROW(assigner.finishOption());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(options.finish());

    EXPECT_EQ(filename, value);
}

} // namespace
