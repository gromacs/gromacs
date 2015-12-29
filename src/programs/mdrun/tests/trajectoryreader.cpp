/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2017, by the GROMACS development team, led by
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
 * \brief Implementions of related classes for tests that want to
 * inspect trajectories produced by mdrun.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include "trajectoryreader.h"

#include <memory>
#include <string>

#include "gromacs/fileio/oenv.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/testasserts.h"

namespace gmx
{
namespace test
{

//! Helper function to obtain resources
t_trxframe *make_trxframe()
{
    t_trxframe *frame;

    snew(frame, 1);
    clear_trxframe(frame, true);

    return frame;
}

//! Helper function to clean up resources
void done_trxframe(t_trxframe *fr)
{
    // Free the contents, then the pointer itself
    sfree(fr->x);
    sfree(fr->v);
    sfree(fr->f);
    sfree(fr->index);
    sfree(fr);
}

// === TrajectoryFrameReader ===

TrajectoryFrameReader::TrajectoryFrameReader(const std::string &filename)
    : filename_(filename),
      trajectoryFileGuard_(nullptr),
      trxframeGuard_(make_trxframe()),
      haveReadFirstFrame_(false),
      haveProbedForNextFrame_(false),
      nextFrameExists_(false)
{
    gmx_output_env_t *oenv;
    output_env_init_default(&oenv);
    oenvGuard_.reset(oenv);
}

bool
TrajectoryFrameReader::readNextFrame()
{
    if (haveProbedForNextFrame_)
    {
        if (nextFrameExists_)
        {
            GMX_THROW(APIError("This frame has already been probed for, it should be used before probing again."));
        }
        else
        {
            GMX_THROW(APIError("This frame has already been probed for, it doesn't exist, so there should not be subsequent attempts to probe for it."));
        }
    }
    haveProbedForNextFrame_ = true;
    // If there's a next frame, read it into trxframe_, and report the result.
    if (!haveReadFirstFrame_)
    {
        t_trxstatus *trajectoryFile;
        int          flags = TRX_READ_X | TRX_READ_V | TRX_READ_F;
        nextFrameExists_ = read_first_frame(oenvGuard_.get(),
                                            &trajectoryFile,
                                            filename_.c_str(),
                                            trxframeGuard_.get(),
                                            flags);
        if (!trajectoryFile)
        {
            GMX_THROW(FileIOError("Could not open trajectory file " + filename_ + " for reading"));
        }
        trajectoryFileGuard_.reset(trajectoryFile);
        haveReadFirstFrame_ = true;
    }
    else
    {
        nextFrameExists_ = read_next_frame(oenvGuard_.get(),
                                           trajectoryFileGuard_.get(),
                                           trxframeGuard_.get());
    }
    return nextFrameExists_;
}

TrajectoryFrame
TrajectoryFrameReader::frame()
{
    TrajectoryFrame frame;

    if (!haveProbedForNextFrame_)
    {
        readNextFrame();
    }
    if (!nextFrameExists_)
    {
        GMX_THROW(APIError("There is no next frame, so there should have been no attempt to use the data, e.g. by reacting to a call to readNextFrame()."));
    }

    // Prepare for reading future frames
    haveProbedForNextFrame_ = false;
    nextFrameExists_        = false;

    // The probe filled trxframeGuard_ with new data, so return it
    frame.frame_ = trxframeGuard_.get();

    if (!frame.frame_->bStep)
    {
        GMX_THROW(APIError("Cannot handle trajectory frame that lacks a step number"));
    }

    if (!frame.frame_->bTime)
    {
        GMX_THROW(APIError("Cannot handle trajectory frame that lacks a time"));
    }

    return frame;
}

// === TrajectoryFrame ===

TrajectoryFrame::TrajectoryFrame() : frame_(nullptr) {};

std::string TrajectoryFrame::getFrameName() const
{
    GMX_RELEASE_ASSERT(frame_, "Cannot get name of invalid frame");
    return formatString("Time %f Step %" GMX_PRId64, frame_->time, frame_->step);
}

// === Free functions ===

static void compareBox(const TrajectoryFrame              &reference,
                       const TrajectoryFrame              &test,
                       const TrajectoryFrameMatchSettings &matchSettings,
                       const FloatingPointTolerance        tolerance)
{
    if (!matchSettings.mustCompareBox)
    {
        return;
    }
    bool canCompareBox = true;
    if (!reference.frame_->bBox)
    {
        ADD_FAILURE() << "Comparing the box was required, but the reference frame did not have one";
        canCompareBox = false;
    }
    if (!test.frame_->bBox)
    {
        ADD_FAILURE() << "Comparing the box was required, but the test frame did not have one";
        canCompareBox = false;
    }
    if (!canCompareBox)
    {
        return;
    }

    // Do the comparing.
    for (int d = 0; d < DIM; ++d)
    {
        for (int dd = 0; dd < DIM; ++dd)
        {
            EXPECT_REAL_EQ_TOL(reference.frame_->box[d][dd], test.frame_->box[d][dd], tolerance)
            << " box[" << d <<"][" << dd << " didn't match between reference frame " << reference.getFrameName() << " and test frame " << test.getFrameName();
        }
    }
}

static void comparePositions(const TrajectoryFrame              &reference,
                             const TrajectoryFrame              &test,
                             const TrajectoryFrameMatchSettings &matchSettings,
                             const FloatingPointTolerance        tolerance)
{
    // Note we don't need to compare bPBC because put_atoms_in_box
    // implements a fallback if nothing specific was set in the
    // trajectory frame.
    bool canHandlePbc = true;
    if (!reference.frame_->bBox)
    {
        if (matchSettings.mustComparePositions)
        {
            ADD_FAILURE() << "Comparing positions required PBC handling, but the reference frame did not have a box";
        }
        canHandlePbc = false;
    }
    if (!test.frame_->bBox)
    {
        if (matchSettings.mustComparePositions)
        {
            ADD_FAILURE() << "Comparing positions required PBC handling, but the test frame did not have a box";
        }
        canHandlePbc = false;
    }

    bool positionsAreMissing = false;
    if (!reference.frame_->bX)
    {
        if (matchSettings.mustComparePositions)
        {
            ADD_FAILURE() << "Comparing positions requires that the reference frame have positions, which it does not";
        }
        positionsAreMissing = true;
    }
    if (!test.frame_->bX)
    {
        if (matchSettings.mustComparePositions)
        {
            ADD_FAILURE() << "Comparing positions requires that the test frame have positions, which it does not";
        }
        positionsAreMissing = true;
    }

    if (!canHandlePbc)
    {
        if (matchSettings.requirePbcHandling)
        {
            ADD_FAILURE() << "Cannot compare positions for the above reason(s)";
            return;
        }
    }
    if (positionsAreMissing)
    {
        if (matchSettings.mustComparePositions)
        {
            ADD_FAILURE() << "Cannot compare positions for the above reason(s)";
        }
        return;
    }

    // If we get here, we have positions, and we know whether we can
    // handle PBC, so we can get to work.
    if ((matchSettings.handlePbcIfPossible || matchSettings.requirePbcHandling) && canHandlePbc)
    {
        // Put all atoms in the box.
        put_atoms_in_box(reference.frame_->ePBC, reference.frame_->box, reference.frame_->natoms, reference.frame_->x);
        put_atoms_in_box(test.frame_->ePBC, test.frame_->box, test.frame_->natoms, test.frame_->x);
    }
    for (int i = 0; i < reference.frame_->natoms && i < test.frame_->natoms; ++i)
    {
        for (int d = 0; d < DIM; ++d)
        {
            EXPECT_REAL_EQ_TOL(reference.frame_->x[i][d], test.frame_->x[i][d], tolerance)
            << " x[" << i << "][" << d <<"] didn't match between reference frame " << reference.getFrameName() << " and test frame " << test.getFrameName();
        }
    }
}

static void compareVelocities(const TrajectoryFrame              &reference,
                              const TrajectoryFrame              &test,
                              const TrajectoryFrameMatchSettings &matchSettings,
                              const FloatingPointTolerance        tolerance)
{
    if (!matchSettings.mustCompareVelocities)
    {
        return;
    }
    bool canCompareBox = true;
    if (!reference.frame_->bF)
    {
        ADD_FAILURE() << "Comparing the velocities was required, but the reference frame did not have any";
        canCompareBox = false;
    }
    if (!test.frame_->bF)
    {
        ADD_FAILURE() << "Comparing the velocities was required, but the test frame did not have any";
        canCompareBox = false;
    }
    if (!canCompareBox)
    {
        return;
    }

    // Do the comparing.
    for (int i = 0; i < reference.frame_->natoms && i < test.frame_->natoms; ++i)
    {
        for (int d = 0; d < DIM; ++d)
        {
            EXPECT_REAL_EQ_TOL(reference.frame_->v[i][d], test.frame_->v[i][d], tolerance)
            << " v[" << i << "][" << d <<"] didn't match between reference frame " << reference.getFrameName() << " and test frame " << test.getFrameName();
        }
    }
}

static void compareForces(const TrajectoryFrame              &reference,
                          const TrajectoryFrame              &test,
                          const TrajectoryFrameMatchSettings &matchSettings,
                          const FloatingPointTolerance        tolerance)
{
    if (!matchSettings.mustCompareForces)
    {
        return;
    }
    bool canCompareBox = true;
    if (!reference.frame_->bV)
    {
        ADD_FAILURE() << "Comparing the forces was required, but the reference frame did not have any";
        canCompareBox = false;
    }
    if (!test.frame_->bV)
    {
        ADD_FAILURE() << "Comparing the forces was required, but the test frame did not have any";
        canCompareBox = false;
    }
    if (!canCompareBox)
    {
        return;
    }

    // Do the comparing.
    for (int i = 0; i < reference.frame_->natoms && i < test.frame_->natoms; ++i)
    {
        for (int d = 0; d < DIM; ++d)
        {
            EXPECT_REAL_EQ_TOL(reference.frame_->f[i][d], test.frame_->f[i][d], tolerance)
            << " f[" << i << "][" << d <<"] didn't match between reference frame " << reference.getFrameName() << " and test frame " << test.getFrameName();
        }
    }
}


void compareFrames(const TrajectoryFrame              &reference,
                   const TrajectoryFrame              &test,
                   const TrajectoryFrameMatchSettings &matchSettings,
                   const FloatingPointTolerance        tolerance)
{
    // The trajectory reader ensures that bStep and bTime are set
    GMX_RELEASE_ASSERT(reference.frame_->bStep, "Reference trajectory must have bStep set");
    GMX_RELEASE_ASSERT(reference.frame_->bTime, "Reference trajectory must have bTime set");
    GMX_RELEASE_ASSERT(test.frame_->bStep, "Test trajectory must have bStep set");
    GMX_RELEASE_ASSERT(test.frame_->bTime, "Test trajectory must have bTime set");

    EXPECT_EQ(reference.frame_->step, test.frame_->step)
    << "step didn't match between reference frame " << reference.getFrameName() << " and test frame " << test.getFrameName();

    EXPECT_EQ(reference.frame_->time, test.frame_->time)
    << "time didn't match between reference frame " << reference.getFrameName() << " and test frame " << test.getFrameName();

    compareBox(reference, test, matchSettings, tolerance);
    comparePositions(reference, test, matchSettings, tolerance);
    compareVelocities(reference, test, matchSettings, tolerance);
    compareForces(reference, test, matchSettings, tolerance);
}

} // namespace
} // namespace
