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
 * Tests for the .mdp nst*out functionality
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include "config.h"

#include <filesystem>
#include <string>

#include <gtest/gtest.h>

#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/h5md/binary_information.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/cmdlinetest.h"
#include "testutils/mpitest.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"
#include "testutils/trajectoryreader.h"

#include "programs/mdrun/tests/mdruncomparison.h"
#include "programs/mdrun/tests/moduletest.h"
#include "programs/mdrun/tests/trajectorycomparison.h"

namespace gmx
{
namespace test
{
namespace
{

// TODO configure these tests so they test all formats for mdrun trajectory writing
#if GMX_USE_TNG

//! Test fixture for mdrun trajectory writing
class TrajectoryWritingTest :
    public gmx::test::MdrunTestFixture,
    public ::testing::WithParamInterface<std::string>
{
public:
    //! The file name of the MDP file
    std::string theMdpFile;

    //! Execute the trajectory writing test
    void runTest(int maxwarn)
    {
        runner_.useStringAsMdpFile(theMdpFile);
        runner_.useTopGroAndNdxFromDatabase("spc-and-methanol");
        runner_.setMaxWarn(maxwarn);
        EXPECT_EQ(0, runner_.callGrompp());

        runner_.fullPrecisionTrajectoryFileName_ =
                fileManager_.getTemporaryFilePath("spc-and-methanol.tng").string();
        runner_.reducedPrecisionTrajectoryFileName_ =
                fileManager_.getTemporaryFilePath("spc-and-methanol-reduced.tng").string();
        ASSERT_EQ(0, runner_.callMdrun());
        // TODO When there is a way to sense something like the
        // output of gmx check, compare the result with that from
        // writing .trr and .xtc and assert the behaviour is
        // correct. Note that TNG will always write the box, even
        // when constant - this will be a source of
        // trajectory-file differences.
    }
};

//! Helper typedef for naming test cases like sentences
typedef TrajectoryWritingTest Trajectories;

/* This test ensures mdrun can write various quantities at various
   frequencies */
TEST_P(Trajectories, ThatDifferInNstxout)
{
    const auto nstxout = GetParam();
    theMdpFile         = gmx::formatString(
            "integrator = md\n"
                    "nsteps = 6\n"
                    "nstxout = %s\n"
                    "nstvout = 2\n"
                    "nstfout = 4\n"
                    "nstxout-compressed = 5\n"
                    "tcoupl = v-rescale\n"
                    "tc-grps = System\n"
                    "tau-t = 1\n"
                    "ref-t = 298\n"
                    "compressed-x-grps = Sol\n",
            nstxout.c_str());
    runTest(0);
}

//! Helper typedef for naming test cases like sentences
typedef TrajectoryWritingTest NptTrajectories;

/* This test ensures mdrun can write trajectories in TNG format from NPT ensembles. */
TEST_P(NptTrajectories, WithDifferentPcoupl)
{
    const auto& pcouple = GetParam();
    int         maxwarn = (pcouple == "Berendsen") ? 1 : 0;
    theMdpFile          = gmx::formatString(
            "integrator = md\n"
                     "nsteps = 2\n"
                     "nstxout = 2\n"
                     "nstvout = 1\n"
                     "pcoupl = %s\n"
                     "tau-p = 1\n"
                     "ref-p = 1\n"
                     "compressibility = 4.5e-5\n"
                     "tcoupl = v-rescale\n"
                     "tc-grps = System\n"
                     "tau-t = 1\n"
                     "ref-t = 298\n",
            pcouple.c_str());
    runTest(maxwarn);
}

// TODO Consider spamming more of the parameter space when we don't
// have to write .mdp and .tpr files to do it.
INSTANTIATE_TEST_SUITE_P(MdrunCanWrite, Trajectories, ::testing::Values("1", "2", "3"));

INSTANTIATE_TEST_SUITE_P(MdrunCanWrite,
                         NptTrajectories,
                         ::testing::Values("no", "Berendsen", "Parrinello-Rahman"));

#endif // GMX_USE_TNG

#if GMX_USE_HDF5
//! Test fixture for mdrun H5md trajectory writing
class H5mdTrajectoryWritingTest : public MdrunTestFixture
{
public:
    void compareTrajectories(const std::filesystem::path& referenceTrajectory,
                             const std::filesystem::path& testTrajectory)
    {
        TrajectoryTolerances tolerances = TrajectoryComparison::s_defaultTrajectoryTolerances;
        if (fn2ftp(referenceTrajectory) == efXTC)
        {
            // XTC stores time as a float
            tolerances.time = defaultFloatTolerance();
        }

        // Build the manager that will present matching pairs of frames to compare
        FramePairManager<TrajectoryFrameReader> trajectoryManager(
                std::make_unique<TrajectoryFrameReader>(referenceTrajectory),
                std::make_unique<TrajectoryFrameReader>(testTrajectory));

        trajectoryManager.compareAllFramePairs<TrajectoryFrame>(
                TrajectoryComparison{ trajectoryMatchSettings_, tolerances });
    }

private:
    const TrajectoryFrameMatchSettings trajectoryMatchSettings_{
        false,
        false,
        false,
        ComparisonConditions::CompareIfReferenceFound,
        ComparisonConditions::CompareIfReferenceFound,
        ComparisonConditions::CompareIfReferenceFound,
        MaxNumFrames::compareAllFrames(),
    };
};

TEST_F(H5mdTrajectoryWritingTest, FullAndReducedPrecisionOutputMatchesTrrAndXtc)
{
    const std::string mdpContent =
            "integrator = md\n"
            "nsteps = 6\n"
            "nstxout = 5\n"
            "nstvout = 2\n"
            "nstfout = 4\n"
            "nstxout-compressed = 3\n";

    runner_.useStringAsMdpFile(mdpContent);
    runner_.useTopGroAndNdxFromDatabase("spc-and-methanol");
    ASSERT_EQ(0, runner_.callGrompp());

    const std::filesystem::path referenceFullPrecisionTrajectory =
            fileManager_.getTemporaryFilePath("spc-and-methanol.trr");
    const std::filesystem::path h5mdFullPrecisionTrajectory =
            fileManager_.getTemporaryFilePath("spc-and-methanol.h5md");
    const std::filesystem::path referenceReducedPrecisionTrajectory =
            fileManager_.getTemporaryFilePath("spc-and-methanol-reduced.xtc");
    const std::filesystem::path h5mdReducedPrecisionTrajectory =
            fileManager_.getTemporaryFilePath("spc-and-methanol-reduced.h5md");

    runner_.fullPrecisionTrajectoryFileName_    = h5mdFullPrecisionTrajectory.string();
    runner_.reducedPrecisionTrajectoryFileName_ = h5mdReducedPrecisionTrajectory.string();
    if (!sz3CompressionIsAvailable())
    {
        GMX_EXPECT_DEATH_IF_SUPPORTED(runner_.callMdrun(), "SZ3");
    }
    else
    {
        ASSERT_EQ(0, runner_.callMdrun());

        runner_.fullPrecisionTrajectoryFileName_    = referenceFullPrecisionTrajectory.string();
        runner_.reducedPrecisionTrajectoryFileName_ = referenceReducedPrecisionTrajectory.string();
        ASSERT_EQ(0, runner_.callMdrun());

        compareTrajectories(referenceFullPrecisionTrajectory, h5mdFullPrecisionTrajectory);
        compareTrajectories(referenceReducedPrecisionTrajectory, h5mdReducedPrecisionTrajectory);
    }
}

#endif // GMX_USE_HDF5

} // namespace
} // namespace test
} // namespace gmx
