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
 * Tests basic file name option implementation.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_options
 */
#include "gmxpre.h"

#include "gromacs/options/filenameoption.h"

#include <string>

#include <gtest/gtest.h>

#include "gromacs/fileio/filetypes.h"
#include "gromacs/options/optionfiletype.h"
#include "gromacs/options/options.h"
#include "gromacs/options/optionsassigner.h"
#include "gromacs/utility/exceptions.h"

#include "testutils/testasserts.h"

namespace gmx
{
namespace test
{
namespace
{

using gmx::FileNameOption;

TEST(FileNameOptionTest, HandlesRequiredDefaultValueWithoutExtension)
{
    gmx::Options options;
    std::string  value;
    ASSERT_NO_THROW_GMX(options.addOption(FileNameOption("f")
                                                  .store(&value)
                                                  .required()
                                                  .filetype(gmx::OptionFileType::GenericData)
                                                  .outputFile()
                                                  .defaultBasename("testfile")));
    EXPECT_EQ("testfile.dat", value);

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW_GMX(assigner.start());
    EXPECT_NO_THROW_GMX(assigner.finish());
    EXPECT_NO_THROW_GMX(options.finish());

    EXPECT_EQ("testfile.dat", value);
}

TEST(FileNameOptionTest, HandlesRequiredOptionWithoutValue)
{
    gmx::Options options;
    std::string  value;
    ASSERT_NO_THROW_GMX(options.addOption(FileNameOption("f")
                                                  .store(&value)
                                                  .required()
                                                  .filetype(gmx::OptionFileType::GenericData)
                                                  .outputFile()
                                                  .defaultBasename("testfile")));
    EXPECT_EQ("testfile.dat", value);

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW_GMX(assigner.start());
    EXPECT_NO_THROW_GMX(assigner.startOption("f"));
    EXPECT_NO_THROW_GMX(assigner.finishOption());
    EXPECT_NO_THROW_GMX(assigner.finish());
    EXPECT_NO_THROW_GMX(options.finish());

    EXPECT_EQ("testfile.dat", value);
}

TEST(FileNameOptionTest, HandlesOptionalUnsetOption)
{
    gmx::Options options;
    std::string  value;
    ASSERT_NO_THROW_GMX(options.addOption(
            FileNameOption("f").store(&value).filetype(gmx::OptionFileType::Trajectory).outputFile().defaultBasename("testfile")));
    EXPECT_TRUE(value.empty());

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW_GMX(assigner.start());
    EXPECT_NO_THROW_GMX(assigner.finish());
    EXPECT_NO_THROW_GMX(options.finish());

    EXPECT_TRUE(value.empty());
}

TEST(FileNameOptionTest, HandlesOptionalDefaultValueWithoutExtension)
{
    gmx::Options options;
    std::string  value;
    ASSERT_NO_THROW_GMX(options.addOption(
            FileNameOption("f").store(&value).filetype(gmx::OptionFileType::AtomIndex).outputFile().defaultBasename("testfile")));
    EXPECT_TRUE(value.empty());

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW_GMX(assigner.start());
    EXPECT_NO_THROW_GMX(assigner.startOption("f"));
    EXPECT_NO_THROW_GMX(assigner.finishOption());
    EXPECT_NO_THROW_GMX(assigner.finish());
    EXPECT_NO_THROW_GMX(options.finish());

    EXPECT_EQ("testfile.ndx", value);
}

TEST(FileNameOptionTest, HandlesRequiredCustomDefaultExtension)
{
    gmx::Options options;
    std::string  value;
    ASSERT_NO_THROW_GMX(options.addOption(FileNameOption("f")
                                                  .store(&value)
                                                  .required()
                                                  .filetype(gmx::OptionFileType::Trajectory)
                                                  .outputFile()
                                                  .defaultBasename("testfile")
                                                  .defaultType(efPDB)));
    EXPECT_EQ("testfile.pdb", value);

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW_GMX(assigner.start());
    EXPECT_NO_THROW_GMX(assigner.finish());
    EXPECT_NO_THROW_GMX(options.finish());

    EXPECT_EQ("testfile.pdb", value);
}

TEST(FileNameOptionTest, HandlesOptionalCustomDefaultExtension)
{
    gmx::Options options;
    std::string  value;
    ASSERT_NO_THROW_GMX(options.addOption(FileNameOption("f")
                                                  .store(&value)
                                                  .filetype(gmx::OptionFileType::Trajectory)
                                                  .outputFile()
                                                  .defaultBasename("testfile")
                                                  .defaultType(efPDB)));
    EXPECT_TRUE(value.empty());

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW_GMX(assigner.start());
    EXPECT_NO_THROW_GMX(assigner.startOption("f"));
    EXPECT_NO_THROW_GMX(assigner.finishOption());
    EXPECT_NO_THROW_GMX(assigner.finish());
    EXPECT_NO_THROW_GMX(options.finish());

    EXPECT_EQ("testfile.pdb", value);
}

TEST(FileNameOptionTest, GivesErrorOnUnknownFileSuffix)
{
    gmx::Options options;
    std::string  value;
    ASSERT_NO_THROW_GMX(options.addOption(
            FileNameOption("f").store(&value).filetype(gmx::OptionFileType::AtomIndex).outputFile()));
    EXPECT_TRUE(value.empty());

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW_GMX(assigner.start());
    EXPECT_NO_THROW_GMX(assigner.startOption("f"));
    EXPECT_THROW_GMX(assigner.appendValue("testfile.foo"), gmx::InvalidInputError);
    EXPECT_NO_THROW_GMX(assigner.finishOption());
    EXPECT_NO_THROW_GMX(assigner.finish());
    EXPECT_NO_THROW_GMX(options.finish());

    EXPECT_TRUE(value.empty());
}

TEST(FileNameOptionTest, GivesErrorOnInvalidFileSuffix)
{
    gmx::Options options;
    std::string  value;
    ASSERT_NO_THROW_GMX(options.addOption(
            FileNameOption("f").store(&value).filetype(gmx::OptionFileType::Trajectory).outputFile()));
    EXPECT_TRUE(value.empty());

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW_GMX(assigner.start());
    EXPECT_NO_THROW_GMX(assigner.startOption("f"));
    EXPECT_THROW_GMX(assigner.appendValue("testfile.dat"), gmx::InvalidInputError);
    EXPECT_NO_THROW_GMX(assigner.finishOption());
    EXPECT_NO_THROW_GMX(assigner.finish());
    EXPECT_NO_THROW_GMX(options.finish());

    EXPECT_TRUE(value.empty());
}

TEST(FileNameOptionTest, HandlesRequiredCsvValueWithoutExtension)
{
    gmx::Options options;
    std::string  value;
    ASSERT_NO_THROW_GMX(options.addOption(FileNameOption("f")
                                                  .store(&value)
                                                  .required()
                                                  .filetype(gmx::OptionFileType::Csv)
                                                  .outputFile()
                                                  .defaultBasename("testfile")));
    EXPECT_EQ("testfile.csv", value);

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW_GMX(assigner.start());
    EXPECT_NO_THROW_GMX(assigner.finish());
    EXPECT_NO_THROW_GMX(options.finish());

    EXPECT_EQ("testfile.csv", value);
}

TEST(FileNameOptionTest, HandlesRequiredCsvOptionWithoutValue)
{
    gmx::Options options;
    std::string  value;
    ASSERT_NO_THROW_GMX(options.addOption(FileNameOption("f")
                                                  .store(&value)
                                                  .required()
                                                  .filetype(gmx::OptionFileType::Csv)
                                                  .outputFile()
                                                  .defaultBasename("testfile")));
    EXPECT_EQ("testfile.csv", value);

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW_GMX(assigner.start());
    EXPECT_NO_THROW_GMX(assigner.startOption("f"));
    EXPECT_NO_THROW_GMX(assigner.finishOption());
    EXPECT_NO_THROW_GMX(assigner.finish());
    EXPECT_NO_THROW_GMX(options.finish());

    EXPECT_EQ("testfile.csv", value);
}

} // namespace
} // namespace test
} // namespace gmx
