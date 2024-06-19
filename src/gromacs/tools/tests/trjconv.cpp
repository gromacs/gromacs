/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2013- The GROMACS Authors
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
 * Tests for gmx trjconv.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 */
#include "gmxpre.h"

#include "gromacs/tools/trjconv.h"

#include "config.h"

#include <cstring>

#include <string>
#include <tuple>

#include <gtest/gtest.h>

#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/cmdlinetest.h"
#include "testutils/simulationdatabase.h"
#include "testutils/stdiohelper.h"
#include "testutils/textblockmatchers.h"
#include "testutils/trajectoryreader.h"

namespace gmx
{
namespace test
{
namespace
{

class TrjconvWithDifferentInputFormats :
    public gmx::test::CommandLineTestBase,
    public ::testing::WithParamInterface<const char*>
{
};

TEST_P(TrjconvWithDifferentInputFormats, WithIndexGroupSubset)
{
    if (!GMX_USE_TNG && std::strstr(GetParam(), ".tng") != nullptr)
    {
        GTEST_SKIP() << "Cannot test TNG reading if TNG support is not configured";
    }
    auto& cmdline = commandLine();

    setInputFile("-s", "spc2.gro");
    setInputFile("-f", GetParam());
    setInputFile("-n", "spc2.ndx");
    std::string outputFile = setOutputFile("-o", "spc-traj.trr", NoTextMatch());

    StdioTestHelper stdioHelper(&fileManager());
    stdioHelper.redirectStringToStdin("SecondWaterMolecule\n");

    ASSERT_EQ(0, gmx_trjconv(cmdline.argc(), cmdline.argv()));

    TrajectoryFrameReader reader(outputFile);
    int                   frameIndex = 0;
    while (reader.readNextFrame())
    {
        TrajectoryFrame frame = reader.frame();
        EXPECT_EQ(frame.x().size(), 3) << "One water molecule should be present in frame " << frameIndex;
        frameIndex++;
    }
}

TEST_P(TrjconvWithDifferentInputFormats, WithoutTopologyFile)
{
    if (!GMX_USE_TNG && std::strstr(GetParam(), ".tng") != nullptr)
    {
        GTEST_SKIP() << "Cannot test TNG reading if TNG support is not configured";
    }
    auto& cmdline = commandLine();

    setInputFile("-f", GetParam());
    setInputFile("-n", "spc2.ndx");
    std::string outputFile = setOutputFile("-o", "spc-traj.trr", NoTextMatch());

    StdioTestHelper stdioHelper(&fileManager());
    stdioHelper.redirectStringToStdin("SecondWaterMolecule\n");

    ASSERT_EQ(0, gmx_trjconv(cmdline.argc(), cmdline.argv()));

    TrajectoryFrameReader reader(outputFile);
    int                   frameIndex = 0;
    while (reader.readNextFrame())
    {
        TrajectoryFrame frame = reader.frame();
        EXPECT_EQ(frame.x().size(), 3) << "One water molecule should be present in frame " << frameIndex;
        frameIndex++;
    }
}

/*! \brief Helper array of input files present in the source repo
 * database. These all have two identical frames of two SPC water
 * molecules, which were generated via trjconv from the .gro
 * version. */
const char* const trajectoryFileNames[] = { "spc2-traj.trr", "spc2-traj.tng", "spc2-traj.xtc",
                                            "spc2-traj.gro", "spc2-traj.pdb", "spc2-traj.g96" };
//! Help GoogleTest name our test cases
std::string nameOfTrjconvWithDifferentInputFormatsTest(const testing::TestParamInfo<const char*>& info)
{
    std::string testName = formatString("file_%s", info.param);

    // Note that the returned names must be unique and may use only
    // alphanumeric ASCII characters. It's not supposed to contain
    // underscores (see the GoogleTest FAQ
    // why-should-test-suite-names-and-test-names-not-contain-underscore),
    // but doing so works for now, is likely to remain so, and makes
    // such test names much more readable.
    testName = replaceAll(testName, "-", "_");
    testName = replaceAll(testName, ".", "_");
    testName = replaceAll(testName, " ", "_");
    return testName;
}

INSTANTIATE_TEST_SUITE_P(Works,
                         TrjconvWithDifferentInputFormats,
                         ::testing::ValuesIn(trajectoryFileNames),
                         nameOfTrjconvWithDifferentInputFormatsTest);

using DumpTestParameters = std::tuple<const char*, double>;

class TrjconvDumpTest :
    public gmx::test::CommandLineTestBase,
    public ::testing::WithParamInterface<DumpTestParameters>
{
};

TEST_P(TrjconvDumpTest, DumpsFrame)
{
    if (!GMX_USE_TNG && std::strstr(std::get<0>(GetParam()), ".tng") != nullptr)
    {
        GTEST_SKIP() << "Cannot test TNG reading if TNG support is not configured";
    }
    auto& cmdline = commandLine();

    setInputFile("-f", std::get<0>(GetParam()));
    const real dumpTime = std::get<1>(GetParam());
    cmdline.addOption("-dump", std::to_string(dumpTime));
    std::string outputFile = setOutputFile("-o", "dumped-frame.trr", gmx::test::NoTextMatch());

    ASSERT_EQ(0, gmx_trjconv(cmdline.argc(), cmdline.argv()));

    // This relies on the input trajectories having frames with times
    // including 0 and 1, so that we test the logic meaningfully.
    real                  expectedFrameTime = (dumpTime < 0.5) ? 0 : 1;
    TrajectoryFrameReader reader(outputFile);
    int                   frameIndex = 0;
    while (reader.readNextFrame())
    {
        TrajectoryFrame frame = reader.frame();
        EXPECT_EQ(frame.time(), expectedFrameTime) << "Dumped frame has expected time";
        frameIndex++;
    }
    EXPECT_EQ(frameIndex, 1) << "A single frame should be dumped";
}

//! Help GoogleTest name our test cases
std::string nameOfTrjconvDumpTest(const testing::TestParamInfo<DumpTestParameters>& info)
{
    std::string testName = formatString(
            "file_%s_"
            "dump_time_%.2f",
            std::get<0>(info.param),
            std::get<1>(info.param));

    // Note that the returned names must be unique and may use only
    // alphanumeric ASCII characters. It's not supposed to contain
    // underscores (see the GoogleTest FAQ
    // why-should-test-suite-names-and-test-names-not-contain-underscore),
    // but doing so works for now, is likely to remain so, and makes
    // such test names much more readable.
    testName = replaceAll(testName, "-", "_");
    testName = replaceAll(testName, ".", "_");
    testName = replaceAll(testName, " ", "_");
    return testName;
}

INSTANTIATE_TEST_SUITE_P(Works,
                         TrjconvDumpTest,
                         ::testing::Combine(::testing::ValuesIn(trajectoryFileNames),
                                            ::testing::Values(-1, 0, 0.3, 1, 999999)),
                         nameOfTrjconvDumpTest);
} // namespace
} // namespace test
} // namespace gmx
