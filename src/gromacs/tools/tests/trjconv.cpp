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

#include "gromacs/fileio/filetypes.h"
#include "gromacs/topology/index.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/message_string_collector.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/cmdlinetest.h"
#include "testutils/simulationdatabase.h"
#include "testutils/stdiohelper.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"
#include "testutils/testmatchers.h"
#include "testutils/textblockmatchers.h"
#include "testutils/tprfilegenerator.h"
#include "testutils/trajectoryreader.h"

//! \brief Helper function for GTest to print tested input and output format.
static void PrintTo(const std::tuple<const char*, GromacsFileType>& info, std::ostream* os)
{
    const auto [inputFile, outputExt] = info;
    *os << "(" << inputFile << ", " << ftp2ext(outputExt) << ")";
}

namespace gmx
{
namespace test
{
namespace
{

using ::testing::Pointwise;

class TrjconvFormatToFormatTest :
    public CommandLineTestBase,
    public ::testing::WithParamInterface<std::tuple<const char*, GromacsFileType>>
{
public:
    TrjconvFormatToFormatTest() :
        inputFile_{ std::get<0>(GetParam()) },
        inputExt_{ static_cast<GromacsFileType>(fn2ftp(inputFile_)) },
        outputExt_{ std::get<1>(GetParam()) }
    {
    }

    /*! \brief Compare data in \p inputFrame to that of \p outputFrame.
     *
     * \param[in] inputFrame Data for frame in original trajectory.
     * \param[in] outputFrame Data for frame in converted trajectory.
     * \param[in] frameIndex  Frame index.
     * \param[in] tolerance   Tolerance for real-value comparisons.
     * \param[in] groupInds   (Optional) Particle indices to compare data for,
     *                        if empty all particle data is compared.
     */
    void compareFrames(const TrajectoryFrame&        inputFrame,
                       const TrajectoryFrame&        outputFrame,
                       const int                     frameIndex,
                       const FloatingPointTolerance& tolerance,
                       const ArrayRef<const int>     groupInds = {})
    {
        EXPECT_EQ(outputFrame.step(), inputFrame.step()) << "Step mismatch in frame " << frameIndex;
        // While time() is returned as a `double`, it may have been read
        // from trajectories which stored it as `float`, so compare it in
        // float precision.
        EXPECT_FLOAT_EQ(outputFrame.time(), inputFrame.time()) << "Time mismatch in frame " << frameIndex;

        EXPECT_THAT(outputFrame.x(), Pointwise(RVecEq(tolerance), valuesForGroupInds(inputFrame.x(), groupInds)))
                << "Particle position mismatch in frame " << frameIndex;

        if (outputFormatSupportsVelocities())
        {
            EXPECT_THAT(outputFrame.v(), Pointwise(RVecEq(tolerance), valuesForGroupInds(inputFrame.v(), groupInds)))
                    << "Particle velocity mismatch in frame " << frameIndex;
        }

        EXPECT_THAT(outputFrame.box().toConstArrayRef(),
                    Pointwise(RealEq(tolerance), inputFrame.box().toConstArrayRef()))
                << "Box mismatch in frame " << frameIndex;
    }

    //! \brief Return the input file name for the current test.
    const char* inputFile() const { return inputFile_; }

    //! \brief Return the input file extension for the current test.
    GromacsFileType inputExt() const { return inputExt_; }

    //! \brief Return the output file extension for the current test.
    GromacsFileType outputExt() const { return outputExt_; }

    /*! \brief Return any common reasons to skip this test if necessary.
     *
     * Common reasons include whether the currently compiled build can read
     * or write the input or output files.
     */
    MessageStringCollector getSkipMessagesIfNecessary() const
    {
        MessageStringCollector messages;
        MessageStringContext   messagesContext(&messages, "Test is being skipped because:");

        messages.appendIf(!GMX_USE_TNG && (inputExt_ == efTNG || outputExt_ == efTNG),
                          "TNG support is not configured");
        messages.appendIf(!GMX_USE_HDF5 && (inputExt_ == efH5MD || outputExt_ == efH5MD),
                          "H5MD support is not configured");
        messages.appendIf(
                GMX_DOUBLE && inputExt_ == efH5MD,
                "Cannot yet read single-precision H5MD trajectory with double-precision build");
        messages.appendIf(outputExt_ == efH5MD, "H5MD output is not yet implemented");
        // TODO: Remove this conditional once a fix is in (see #5747).
        messages.appendIf(inputExt_ == efTNG && outputExt_ == efTNG,
                          "TNG to TNG conversion is currently broken");

        return messages;
    }

    //! \brief Return true if the current parameter output format supports velocities.
    bool outputFormatSupportsVelocities() const
    {
        return outputExt_ == efGRO || outputExt_ == efTRR || outputExt_ == efH5MD
               || outputExt_ == efTNG || outputExt_ == efG96;
    }

private:
    /*! \brief Return a view into \p values for the given \p groupInds (unchanged if empty).
     */
    ArrayRef<const RVec> valuesForGroupInds(const ArrayRef<const RVec> values,
                                            const ArrayRef<const int>  groupInds)
    {
        if (values.empty() || groupInds.empty())
        {
            return values;
        }

        valuesForGroupInds_.resize(groupInds.size());
        for (int i = 0; i < gmx::ssize(groupInds); ++i)
        {
            valuesForGroupInds_[i] = values[groupInds[i]];
        }

        return valuesForGroupInds_;
    }

    //! Input file name for the current test.
    const char* inputFile_;

    //! Input file extension for the current test.
    const GromacsFileType inputExt_;

    //! Output file extension for the current test.
    const GromacsFileType outputExt_;

    //! Persistent buffer to store values for group indices in.
    std::vector<RVec> valuesForGroupInds_;
};

TEST_P(TrjconvFormatToFormatTest, FullSystemWithTprTopology)
{
    if (const MessageStringCollector skipMessages = getSkipMessagesIfNecessary(); !skipMessages.isEmpty())
    {
        GTEST_SKIP() << skipMessages.toString();
    }
    auto&             cmdline = commandLine();
    TprAndFileManager tprManager("spc2-traj");

    setInputFile("-f", inputFile());
    setInputFile("-s", tprManager.tprName());
    const std::filesystem::path outputFile =
            setOutputFile("-o", formatString("trjconv-out.%s", ftp2ext(outputExt())), NoTextMatch());

    StdioTestHelper stdioHelper(&fileManager());
    stdioHelper.redirectStringToStdin("System\n");

    ASSERT_EQ(0, gmx_trjconv(cmdline.argc(), cmdline.argv()));

    TrajectoryFrameReader inputFileReader(fileManager().getInputFilePath(inputFile()));
    TrajectoryFrameReader outputFileReader(outputFile);
    int                   frameIndex = 0;
    while (inputFileReader.readNextFrame())
    {
        ASSERT_TRUE(outputFileReader.readNextFrame()) << "Missing frame " << frameIndex << " in output";
        compareFrames(
                inputFileReader.frame(), outputFileReader.frame(), frameIndex, absoluteTolerance(0.001));
        ++frameIndex;
    }
}

TEST_P(TrjconvFormatToFormatTest, FullSystemWithGroTopology)
{
    if (const MessageStringCollector skipMessages = getSkipMessagesIfNecessary(); !skipMessages.isEmpty())
    {
        GTEST_SKIP() << skipMessages.toString();
    }
    auto& cmdline = commandLine();

    setInputFile("-f", inputFile());
    setInputFile("-s", "spc2.gro");
    const std::filesystem::path outputFile =
            setOutputFile("-o", formatString("trjconv-out.%s", ftp2ext(outputExt())), NoTextMatch());

    StdioTestHelper stdioHelper(&fileManager());
    stdioHelper.redirectStringToStdin("System\n");

    ASSERT_EQ(0, gmx_trjconv(cmdline.argc(), cmdline.argv()));

    TrajectoryFrameReader inputFileReader(fileManager().getInputFilePath(inputFile()));
    TrajectoryFrameReader outputFileReader(outputFile);
    int                   frameIndex = 0;
    while (inputFileReader.readNextFrame())
    {
        ASSERT_TRUE(outputFileReader.readNextFrame()) << "Missing frame " << frameIndex << " in output";
        compareFrames(
                inputFileReader.frame(), outputFileReader.frame(), frameIndex, absoluteTolerance(0.001));
        ++frameIndex;
    }
}

TEST_P(TrjconvFormatToFormatTest, FullSystemWithoutTopology)
{
    MessageStringCollector skipMessages = getSkipMessagesIfNecessary();
    MessageStringContext skipMessagesExtraContext(&skipMessages, "Test is being skipped because:");
    skipMessages.appendIf(outputExt() == efPDB || outputExt() == efGRO || outputExt() == efTNG,
                          "PDB, GRO and TNG output requires a present topology");
    if (!skipMessages.isEmpty())
    {
        GTEST_SKIP() << skipMessages.toString();
    }
    auto& cmdline = commandLine();

    setInputFile("-f", inputFile());
    const std::filesystem::path outputFile =
            setOutputFile("-o", formatString("trjconv-out.%s", ftp2ext(outputExt())), NoTextMatch());

    StdioTestHelper stdioHelper(&fileManager());
    stdioHelper.redirectStringToStdin("System\n");

    ASSERT_EQ(0, gmx_trjconv(cmdline.argc(), cmdline.argv()));

    TrajectoryFrameReader inputFileReader(fileManager().getInputFilePath(inputFile()));
    TrajectoryFrameReader outputFileReader(outputFile);
    int                   frameIndex = 0;
    while (inputFileReader.readNextFrame())
    {
        ASSERT_TRUE(outputFileReader.readNextFrame()) << "Missing frame " << frameIndex << " in output";
        compareFrames(
                inputFileReader.frame(), outputFileReader.frame(), frameIndex, absoluteTolerance(0.001));
        ++frameIndex;
    }
}

TEST_P(TrjconvFormatToFormatTest, WithIndexGroupSubset)
{
    if (const MessageStringCollector skipMessages = getSkipMessagesIfNecessary(); !skipMessages.isEmpty())
    {
        GTEST_SKIP() << skipMessages.toString();
    }
    auto&             cmdline = commandLine();
    TprAndFileManager tprManager("spc2-traj");

    setInputFile("-f", inputFile());
    setInputFile("-s", tprManager.tprName());
    const std::filesystem::path outputFile =
            setOutputFile("-o", formatString("trjconv-out.%s", ftp2ext(outputExt())), NoTextMatch());

    constexpr char indexGroupFile[] = "spc2.ndx";
    constexpr char indexGroupName[] = "SecondWaterMolecule";
    setInputFile("-n", indexGroupFile);
    StdioTestHelper stdioHelper(&fileManager());
    stdioHelper.redirectStringToStdin(formatString("%s\n", indexGroupName).c_str());

    const std::vector<IndexGroup> indexGroups =
            init_index(fileManager().getInputFilePath(indexGroupFile));
    const int groupIndex = find_group(indexGroupName, indexGroups);
    ASSERT_NE(groupIndex, -1) << "Sanity check failed: could not find index group in file";

    ASSERT_EQ(0, gmx_trjconv(cmdline.argc(), cmdline.argv()));

    TrajectoryFrameReader inputFileReader(fileManager().getInputFilePath(inputFile()));
    TrajectoryFrameReader outputFileReader(outputFile);
    int                   frameIndex = 0;
    while (inputFileReader.readNextFrame())
    {
        ASSERT_TRUE(outputFileReader.readNextFrame()) << "Missing frame " << frameIndex << " in output";
        compareFrames(inputFileReader.frame(),
                      outputFileReader.frame(),
                      frameIndex,
                      absoluteTolerance(0.001),
                      indexGroups[groupIndex].particleIndices);
        ++frameIndex;
    }
}

//! Help GoogleTest name our test cases.
std::string nameOfTest(const testing::TestParamInfo<std::tuple<const char*, GromacsFileType>>& info)
{
    const auto [inputFile, outputExt] = info.param;
    std::string testName              = formatString("%s_to_%s", inputFile, ftp2ext(outputExt));

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

/*! \brief Helper array of input files present in the source repo
 * database. These all have two identical frames of two SPC water
 * molecules, which were generated via trjconv from the .gro
 * version. */
const char* const g_fromTrajectoryFiles[] = { "spc2-traj.trr",  "spc2-traj-double.trr",
                                              "spc2-traj.xtc",  "spc2-traj.gro",
                                              "spc2-traj.pdb",  "spc2-traj.g96",
                                              "spc2-traj.h5md", "spc2-traj.tng" };

//! List of file types to convert to.
static const GromacsFileType g_toTrajectoryFormats[] = {
    efXTC, efTRR, efPDB, efGRO, efG96, efH5MD, efTNG,
};

INSTANTIATE_TEST_SUITE_P(Works,
                         TrjconvFormatToFormatTest,
                         ::testing::Combine(testing::ValuesIn(g_fromTrajectoryFiles),
                                            testing::ValuesIn(g_toTrajectoryFormats)),
                         nameOfTest);

using TrjconvFormatToG96Test = TrjconvFormatToFormatTest;

TEST_P(TrjconvFormatToG96Test, WithoutTopologyFileG96InvalidCharacters)
{
    if (const MessageStringCollector skipMessages = getSkipMessagesIfNecessary(); !skipMessages.isEmpty())
    {
        GTEST_SKIP() << skipMessages.toString();
    }
    auto& cmdline = commandLine();

    // Skip lines with only numbers, a lazy way to avoid handling floating point
    auto matcher = FilteringExactTextMatch({ "[0-9\\. ]+" }, false, true);
    setInputFile("-f", inputFile());
    setOutputFile("-o", "spc-traj.g96", matcher);

    ASSERT_EQ(0, gmx_trjconv(cmdline.argc(), cmdline.argv()));

    checkOutputFiles();
}

INSTANTIATE_TEST_SUITE_P(Works,
                         TrjconvFormatToG96Test,
                         ::testing::Combine(::testing::ValuesIn(g_fromTrajectoryFiles),
                                            ::testing::Values(efG96)),
                         nameOfTest);

using DumpTestParameters = std::tuple<const char*, double>;

class TrjconvDumpTest : public CommandLineTestBase, public ::testing::WithParamInterface<DumpTestParameters>
{
};

TEST_P(TrjconvDumpTest, DumpsFrame)
{
    if (!GMX_USE_TNG && std::strstr(std::get<0>(GetParam()), ".tng") != nullptr)
    {
        GTEST_SKIP() << "Cannot test TNG reading if TNG support is not configured";
    }
    if (!GMX_USE_HDF5 && std::strstr(std::get<0>(GetParam()), ".h5md") != nullptr)
    {
        GTEST_SKIP() << "Cannot test H5MD reading if H5MD support is not configured";
    }
    if (GMX_USE_HDF5 && GMX_DOUBLE && std::strstr(std::get<0>(GetParam()), ".h5md") != nullptr)
    {
        GTEST_SKIP() << "Cannot yet test H5MD reading from double precision since trajectory file "
                        "is single precision";
    }
    auto& cmdline = commandLine();

    setInputFile("-f", std::get<0>(GetParam()));
    const real dumpTime = std::get<1>(GetParam());
    cmdline.addOption("-dump", std::to_string(dumpTime));
    const std::filesystem::path outputFile = setOutputFile("-o", "dumped-frame.trr", NoTextMatch());

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
                         ::testing::Combine(::testing::ValuesIn(g_fromTrajectoryFiles),
                                            ::testing::Values(-1, 0, 0.3, 1, 999999)),
                         nameOfTrjconvDumpTest);

//! \brief Test fixture for gmx trjconv -pbc (option) tests
class TrjconvPbcTest : public CommandLineTestBase, public ::testing::WithParamInterface<const char*>
{
};

TEST_P(TrjconvPbcTest, MoleculeCrossingBox)
{
    auto&             cmdline = commandLine();
    TprAndFileManager tprManager("spc2-traj");
    StdioTestHelper   stdioHelper(&fileManager());

    SCOPED_TRACE("Test -pbc <GetParam()> option for trajectory with molecule broken across PBC");
    setInputFile("-f", fileManager().getInputFilePath("spc2-traj-shifted-across-pbc-boundaries.gro"));
    setInputFile("-s", tprManager.tprName());
    setInputFile("-n", "spc2.ndx");
    setOutputFile("-o", "out.gro", ExactTextMatch());
    cmdline.addOption("-pbc", GetParam());
    cmdline.addOption("-dump", 0); // Only test the first frame

    stdioHelper.redirectStringToStdin("System");
    ASSERT_EQ(0, gmx_trjconv(cmdline.argc(), cmdline.argv()));
    checkOutputFiles();
}

//! Help GoogleTest name our test cases.
std::string nameOfPbcTest(const testing::TestParamInfo<const char*>& info)
{
    std::string testName = info.param;

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
                         TrjconvPbcTest,
                         ::testing::Values("mol", "res", "atom", "nojump", "whole"),
                         nameOfPbcTest);

} // namespace
} // namespace test
} // namespace gmx
