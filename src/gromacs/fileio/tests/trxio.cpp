/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2026- The GROMACS Authors
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
 * Tests for trajectory file I/O routines
 *
 * \author Petter Johansson <pettjoha@kth.se>
 * \ingroup module_fileio
 */
#include "gmxpre.h"

#include "gromacs/fileio/trxio.h"

#include <algorithm>

#include <gtest/gtest.h>

#include "gromacs/coordinateio/coordinatefile.h"
#include "gromacs/coordinateio/requirements.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/selection/selection.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/vec.h"
#include "gromacs/utility/vectypes.h"

#include "testutils/generate_frame_data.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"
#include "testutils/testmatchers.h"
#include "testutils/trajectoryreader.h"

struct t_fileio;

// TODO: These should really appear somewhere centralized.
/*! \brief
 * Google Test formatter for GromacsFileType values.
 */
static void PrintTo(const GromacsFileType& ftp, std::ostream* os)
{
    *os << "'" << ftp2ext(ftp) << "'";
}

namespace gmx
{
namespace test
{
namespace
{

using ::testing::Pointwise;

class TrxWriteFrameTest : public ::testing::Test, public ::testing::WithParamInterface<GromacsFileType>
{
public:
    gmx::test::TestFileManager fileManager_;
};

TEST_P(TrxWriteFrameTest, WriteTrxFrame_TrajectoryDataOnly)
{
    const std::filesystem::path trajectoryFileName =
            fileManager_.getTemporaryFilePath(formatString("traj.%s", ftp2ext(GetParam())));

    constexpr int     nframes             = 3;
    constexpr int     natoms              = 10;
    constexpr int64_t steps[nframes]      = { 42, -5, 0 };
    constexpr real    times[nframes]      = { 50.5, -501.0, 0.0 };
    constexpr real    precisions[nframes] = { 100.0, 1000.0, 10000.0 };

    // Generate per-frame data for box and trajectory
    std::array<matrix, nframes>                   boxes;
    std::array<std::array<rvec, natoms>, nframes> positions, velocities, forces;

    std::for_each(boxes.begin(), boxes.end(), MatrixFrameDataGenerator());
    std::for_each(positions.begin(),
                  positions.end(),
                  TrajectoryFrameDataGenerator(TrajectoryFrameMode::Positive));
    std::for_each(velocities.begin(),
                  velocities.end(),
                  TrajectoryFrameDataGenerator(TrajectoryFrameMode::Negative));
    std::for_each(forces.begin(), forces.end(), TrajectoryFrameDataGenerator(TrajectoryFrameMode::Alternating));

    {
        SCOPED_TRACE("Create a new file and write the frames into it");

        t_trxframe fr;
        clear_trxframe(&fr, true);
        fr.natoms = natoms;
        fr.bStep  = true;
        fr.bTime  = true;
        fr.bPrec  = true;
        fr.bBox   = true;
        fr.bX     = true;
        fr.bV     = true;
        fr.bF     = true;

        auto writer = createTrajectoryFrameWriter(nullptr, {}, trajectoryFileName, nullptr, {});
        for (int i = 0; i < nframes; ++i)
        {
            SCOPED_TRACE(formatString("write_trxframe for frame %d", i));
            copy_mat(boxes[i], fr.box);
            fr.time = times[i];
            fr.step = steps[i];
            fr.prec = precisions[i];
            // make shallow copies of the reference data to avoid extra allocations
            fr.x = positions[i].data();
            fr.v = velocities[i].data();
            fr.f = forces[i].data();

            writer->prepareAndWriteFrame(i, fr);
        }
    }
    {
        SCOPED_TRACE("Open the file and read the frames from it");
        TrajectoryFrameReader reader(trajectoryFileName);

        for (int frame = 0; frame < nframes; ++frame)
        {
            ASSERT_TRUE(reader.readNextFrame());
            const TrajectoryFrame fr = reader.frame();

            EXPECT_REAL_EQ(times[frame], fr.time());
            EXPECT_EQ(steps[frame], fr.step());

            switch (GetParam())
            {
                case efXTC:
                    ASSERT_THAT(fr.positionPrecision(), ::testing::Optional(precisions[frame]));
                    EXPECT_THAT(fr.x(),
                                Pointwise(RVecEq(absoluteTolerance(1.0 / fr.positionPrecision().value())),
                                          positions[frame]));
                    break;
                case efTRR:
                    EXPECT_FALSE(fr.positionPrecision().has_value());
                    EXPECT_THAT(fr.x(),
                                Pointwise(RVecEq(relativeToleranceAsFloatingPoint(1.0, 1e-7)),
                                          positions[frame]));
                    EXPECT_THAT(fr.v(),
                                Pointwise(RVecEq(relativeToleranceAsFloatingPoint(1.0, 1e-7)),
                                          velocities[frame]));
                    EXPECT_THAT(fr.f(),
                                Pointwise(RVecEq(relativeToleranceAsFloatingPoint(1.0, 1e-7)),
                                          forces[frame]));
                    break;
                default: ASSERT_TRUE(false) << "Sanity check failed: Uncovered file type for test";
            }

            ASSERT_TRUE(fr.hasBox());
            EXPECT_THAT(fr.box().toConstArrayRef(),
                        Pointwise(RealEq(relativeToleranceAsFloatingPoint(1.0, 1e-7)),
                                  arrayRefFromArray(boxes[frame][0], DIM * DIM)));
        }
        {
            SCOPED_TRACE("Reading subsequent frames should now fail gracefully");
            ASSERT_FALSE(reader.readNextFrame());
        }
    }
}

INSTANTIATE_TEST_SUITE_P(WithDifferentFormats, TrxWriteFrameTest, ::testing::Values(efXTC, efTRR));

//! \brief Helper function for to construct test and reference data file names.
static std::string nameOfTest(const ::testing::TestParamInfo<const char*>& info)
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

class TrxReadFrameTest : public ::testing::Test, public ::testing::WithParamInterface<const char*>
{
public:
    gmx::test::TestFileManager fileManager_;
};

// We here check that read_first/next_frame works and read trajectories
// in our simulation data base to ensure that the functions remain stable
// over time.
TEST_P(TrxReadFrameTest, ReadTrajectoryFromSimulationDataBase)
{
    TestReferenceData    data;
    TestReferenceChecker checker(data.rootChecker());
    checker.setDefaultTolerance(relativeToleranceAsFloatingPoint(1.0, 1e-7));

    t_trxstatus* status;
    t_trxframe*  fr;
    snew(fr, 1);
    gmx_output_env_t* oenv;
    output_env_init_default(&oenv);

    ASSERT_TRUE(read_first_frame(
            oenv, &status, fileManager_.getInputFilePath(GetParam()), fr, TRX_READ_X | TRX_READ_V | TRX_READ_F));
    checker.checkInteger(fr->natoms, "natoms");

    int frame = 0;
    do
    {
        // Checking the t_atoms structure is beyond this test, only check for its (non)presence
        checker.checkBoolean(fr->bAtoms, formatString("bAtoms (frame %d)", frame).c_str());
        if (fr->bStep)
        {
            checker.checkInteger(fr->step, formatString("step (frame %d)", frame).c_str());
        }
        if (fr->bTime)
        {
            checker.checkReal(fr->time, formatString("time (frame %d)", frame).c_str());
        }
        if (fr->bLambda)
        {
            checker.checkReal(fr->lambda, formatString("lambda (frame %d)", frame).c_str());
        }
        if (fr->bFepState)
        {
            checker.checkInteger(fr->fep_state, formatString("fep_state (frame %d)", frame).c_str());
        }
        if (fr->bPrec)
        {
            checker.checkReal(fr->prec, formatString("prec (frame %d)", frame).c_str());
        }
        if (fr->bX)
        {
            checker.checkSequence(constArrayRefFromArray(fr->x, fr->natoms),
                                  formatString("x (frame %d)", frame).c_str());
        }
        if (fr->bV)
        {
            checker.checkSequence(constArrayRefFromArray(fr->v, fr->natoms),
                                  formatString("v (frame %d)", frame).c_str());
        }
        if (fr->bF)
        {
            checker.checkSequence(constArrayRefFromArray(fr->f, fr->natoms),
                                  formatString("f (frame %d)", frame).c_str());
        }
        if (fr->bBox)
        {
            checker.checkSequence(constArrayRefFromArray(fr->box, 3),
                                  formatString("box (frame %d)", frame).c_str());
        }
        if (fr->bPBC)
        {
            checker.checkInteger(static_cast<int>(fr->pbcType),
                                 formatString("pbcType (frame %d)", frame).c_str());
        }
        if (fr->bIndex)
        {
            checker.checkSequence(constArrayRefFromArray(fr->index, fr->natoms),
                                  formatString("index (frame %d)", frame).c_str());
        }
        ++frame;
    } while (read_next_frame(oenv, status, fr));

    output_env_done(oenv);
    done_frame(fr);
    sfree(fr);
    close_trx(status);
}

const char* g_trajectoryFilesToTest[] = {
    "spc2-traj.pdb",
    "spc2-traj.trr",
    "spc2-traj-double.trr",
    "spc2-traj.xtc",
    "spc2-traj.gro",
    "spc2-traj.g96",
#if GMX_USE_HDF5 && !GMX_DOUBLE
    // NOTE: Currently the build precision must match that of the file (see #5474)
    "spc2-traj.h5md",
#endif
#if GMX_USE_TNG
    "spc2-traj.tng",
#endif
};
INSTANTIATE_TEST_SUITE_P(WithDifferentFormats,
                         TrxReadFrameTest,
                         ::testing::ValuesIn(g_trajectoryFilesToTest),
                         nameOfTest);

} // namespace
} // namespace test
} // namespace gmx
