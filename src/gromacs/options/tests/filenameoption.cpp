/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014, by the GROMACS development team, led by
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
 * Tests basic file name option implementation.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_options
 */
#include "gmxpre.h"

#include "gromacs/options/filenameoption.h"

#include <gtest/gtest.h>

#include "gromacs/fileio/filenm.h"
#include "gromacs/options/options.h"
#include "gromacs/options/optionsassigner.h"
#include "gromacs/utility/exceptions.h"

#include "testutils/testasserts.h"

namespace
{

using gmx::FileNameOption;

TEST(FileNameOptionTest, HandlesRequiredDefaultValueWithoutExtension)
{
    gmx::Options           options(NULL, NULL);
    std::string            value;
    ASSERT_NO_THROW_GMX(options.addOption(
                                FileNameOption("f").store(&value).required()
                                    .filetype(gmx::eftGenericData).outputFile()
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
    gmx::Options           options(NULL, NULL);
    std::string            value;
    ASSERT_NO_THROW_GMX(options.addOption(
                                FileNameOption("f").store(&value).required()
                                    .filetype(gmx::eftGenericData).outputFile()
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
    gmx::Options           options(NULL, NULL);
    std::string            value;
    ASSERT_NO_THROW_GMX(options.addOption(
                                FileNameOption("f").store(&value)
                                    .filetype(gmx::eftTrajectory).outputFile()
                                    .defaultBasename("testfile")));
    EXPECT_TRUE(value.empty());

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW_GMX(assigner.start());
    EXPECT_NO_THROW_GMX(assigner.finish());
    EXPECT_NO_THROW_GMX(options.finish());

    EXPECT_TRUE(value.empty());
}

TEST(FileNameOptionTest, HandlesOptionalDefaultValueWithoutExtension)
{
    gmx::Options           options(NULL, NULL);
    std::string            value;
    ASSERT_NO_THROW_GMX(options.addOption(
                                FileNameOption("f").store(&value)
                                    .filetype(gmx::eftIndex).outputFile()
                                    .defaultBasename("testfile")));
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
    gmx::Options           options(NULL, NULL);
    std::string            value;
    ASSERT_NO_THROW_GMX(options.addOption(
                                FileNameOption("f").store(&value).required()
                                    .filetype(gmx::eftTrajectory).outputFile()
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
    gmx::Options           options(NULL, NULL);
    std::string            value;
    ASSERT_NO_THROW_GMX(options.addOption(
                                FileNameOption("f").store(&value)
                                    .filetype(gmx::eftTrajectory).outputFile()
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
    gmx::Options           options(NULL, NULL);
    std::string            value;
    ASSERT_NO_THROW_GMX(options.addOption(
                                FileNameOption("f").store(&value)
                                    .filetype(gmx::eftIndex).outputFile()));
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
    gmx::Options           options(NULL, NULL);
    std::string            value;
    ASSERT_NO_THROW_GMX(options.addOption(
                                FileNameOption("f").store(&value)
                                    .filetype(gmx::eftTrajectory).outputFile()));
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

} // namespace
