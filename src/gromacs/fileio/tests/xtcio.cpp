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
 * Tests for XTC file I/O routines
 *
 * \author Petter Johansson <pettjoha@kth.se>
 * \ingroup module_fileio
 */
#include "gmxpre.h"

#include "gromacs/fileio/xtcio.h"

#include <gtest/gtest.h>

#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/vectypes.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"
#include "testutils/testmatchers.h"

struct t_fileio;

namespace gmx
{
namespace test
{
namespace
{

using ::testing::Pointwise;

class XtcTest : public ::testing::Test
{
public:
    XtcTest() {}
    gmx::test::TestFileManager fileManager_;
};

TEST_F(XtcTest, WriteAndRead)
{
    const std::filesystem::path fn = fileManager_.getTemporaryFilePath("traj.xtc");

    constexpr int     nframes = 3;
    constexpr int     natoms  = 10; // must be > 9 to activate compression (lets us check precision)
    constexpr int64_t steps[nframes]      = { 42, -5, 0 };
    constexpr real    times[nframes]      = { 50.5, -501.0, 0.0 };
    constexpr real    precisions[nframes] = { 1000, 100, 1000 };
    constexpr matrix  boxes[nframes]      = {
        { { 1.0, 2.0, 3.0 }, { 4.0, 5.0, 6.0 }, { 7.0, 8.0, 9.0 } },
        { { 10.0, 20.0, 30.0 }, { 40.0, 50.0, 60.0 }, { 70.0, 80.0, 90.0 } },
        { { 100.0, 200.0, 300.0 }, { 400.0, 500.0, 600.0 }, { 700.0, 800.0, 900.0 } }
    };

    // Create positions for frames to write and then read
    std::array<std::array<RVec, natoms>, nframes> positions;
    for (int frame = 0; frame < nframes; ++frame)
    {
        for (int i = 0; i < natoms; ++i)
        {
            const real f        = frame + (0.1 * i);
            positions[frame][i] = { f + 0.01_real, f + 0.02_real, f + 0.03_real };
        }
    }

    {
        SCOPED_TRACE("Create a new file and write the frames into it");
        t_fileio* fio = open_xtc(fn, "w");
        for (int frame = 0; frame < nframes; ++frame)
        {
            SCOPED_TRACE(formatString("write_xtc for frame %d", frame));
            ASSERT_TRUE(write_xtc(fio,
                                  natoms,
                                  steps[frame],
                                  times[frame],
                                  boxes[frame],
                                  as_rvec_array(positions[frame].data()),
                                  precisions[frame]));
        }
        close_xtc(fio);
    }
    {
        SCOPED_TRACE("Open the file and read the frames from it");
        t_fileio* fio = open_xtc(fn, "r");
        int       natomsRead;
        int64_t   step;
        real      time;
        real      precision;
        matrix    box;
        rvec*     x;
        gmx_bool  bOk;

        for (int frame = 0; frame < nframes; ++frame)
        {
            SCOPED_TRACE(frame == 0 ? "read_first_xtc for frame 0"
                                    : formatString("read_next_xtc for frame %d", frame));
            if (frame == 0)
            {
                ASSERT_TRUE(read_first_xtc(fio, &natomsRead, &step, &time, box, &x, &precision, &bOk));
            }
            else
            {
                ASSERT_TRUE(read_next_xtc(fio, natomsRead, &step, &time, box, x, &precision, &bOk));
            }

            ASSERT_TRUE(bOk);
            ASSERT_EQ(natomsRead, natoms);
            EXPECT_EQ(steps[frame], step);
            EXPECT_REAL_EQ(times[frame], time);
            // if the read precision is -1, compression was not activated (see note about natoms)
            EXPECT_REAL_EQ(precisions[frame], precision);
            EXPECT_THAT(arrayRefFromArray(box, DIM),
                        Pointwise(RVecEq(defaultRealTolerance()), arrayRefFromArray(boxes[frame], DIM)));
            EXPECT_THAT(arrayRefFromArray(x, natoms),
                        Pointwise(RVecEq(absoluteTolerance(1.0 / precision)), positions[frame]));
        }
        {
            SCOPED_TRACE("Reading subsequent frames should now fail gracefully");
            ASSERT_FALSE(read_next_xtc(fio, natomsRead, &step, &time, box, x, &precision, &bOk));
        }

        sfree(x);
        close_xtc(fio);
    }
}

TEST_F(XtcTest, ReadTrajectoryFromSimulationDataBase)
{
    TestReferenceData    data;
    TestReferenceChecker checker(data.rootChecker());
    checker.setDefaultTolerance(relativeToleranceAsFloatingPoint(1.0, 1e-7));

    t_fileio* fio = open_xtc(fileManager_.getInputFilePath("spc2-traj.xtc"), "r");

    int      natoms;
    int64_t  step;
    real     time;
    real     prec;
    matrix   box;
    rvec*    x;
    gmx_bool bOk;

    ASSERT_TRUE(read_first_xtc(fio, &natoms, &step, &time, box, &x, &prec, &bOk));
    checker.checkInteger(natoms, "natoms");

    int frame = 0;
    do
    {
        checker.checkInteger(step, formatString("step (frame %d)", frame).c_str());
        checker.checkReal(time, formatString("time (frame %d)", frame).c_str());
        checker.checkReal(prec, formatString("prec (frame %d)", frame).c_str());
        checker.checkSequence(constArrayRefFromArray(box, 3),
                              formatString("box (frame %d)", frame).c_str());
        checker.checkSequence(constArrayRefFromArray(x, natoms),
                              formatString("x (frame %d)", frame).c_str());
        ++frame;
    } while (read_next_xtc(fio, natoms, &step, &time, box, x, &prec, &bOk));

    sfree(x);
    close_xtc(fio);
}

} // namespace
} // namespace test
} // namespace gmx
