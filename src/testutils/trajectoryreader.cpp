/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2016- The GROMACS Authors
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
 * \brief Implementions of related classes for tests that want to
 * inspect trajectories produced by mdrun.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include "testutils/trajectoryreader.h"

#include <filesystem>
#include <memory>
#include <string>

#include "gromacs/fileio/oenv.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/smalloc.h"

#include "testutils/testasserts.h"
#include "testutils/testmatchers.h"

struct gmx_output_env_t;

namespace gmx
{
namespace test
{

//! Helper function to obtain resources
static t_trxframe* make_trxframe()
{
    t_trxframe* frame;

    snew(frame, 1);
    clear_trxframe(frame, true);

    return frame;
}

//! Helper function to clean up resources
void done_trxframe(t_trxframe* fr)
{
    // Free the contents, then the pointer itself
    sfree(fr->x);
    sfree(fr->v);
    sfree(fr->f);
    sfree(fr->index);
    sfree(fr);
}

TrajectoryFrameReader::TrajectoryFrameReader(const std::string& filename) :
    filename_(filename),
    trajectoryFileGuard_(nullptr),
    trxframeGuard_(make_trxframe()),
    haveReadFirstFrame_(false),
    haveProbedForNextFrame_(false),
    nextFrameExists_(false)
{
    gmx_output_env_t* oenv;
    output_env_init_default(&oenv);
    oenvGuard_.reset(oenv);
}

bool TrajectoryFrameReader::readNextFrame()
{
    if (haveProbedForNextFrame_)
    {
        if (nextFrameExists_)
        {
            GMX_THROW(
                    APIError("This frame has already been probed for, it should be used before "
                             "probing again."));
        }
        else
        {
            GMX_THROW(
                    APIError("This frame has already been probed for, it doesn't exist, so there "
                             "should not be subsequent attempts to probe for it."));
        }
    }
    haveProbedForNextFrame_ = true;
    // If there's a next frame, read it into trxframe_, and report the result.
    if (!haveReadFirstFrame_)
    {
        t_trxstatus* trajectoryFile;
        int          flags = TRX_READ_X | TRX_READ_V | TRX_READ_F;
        nextFrameExists_   = read_first_frame(
                oenvGuard_.get(), &trajectoryFile, filename_.c_str(), trxframeGuard_.get(), flags);
        if (!trajectoryFile)
        {
            GMX_THROW(FileIOError("Could not open trajectory file " + filename_ + " for reading"));
        }
        trajectoryFileGuard_.reset(trajectoryFile);
        haveReadFirstFrame_ = true;
    }
    else
    {
        nextFrameExists_ =
                read_next_frame(oenvGuard_.get(), trajectoryFileGuard_.get(), trxframeGuard_.get());
    }
    return nextFrameExists_;
}

TrajectoryFrame TrajectoryFrameReader::frame()
{
    if (!haveProbedForNextFrame_)
    {
        readNextFrame();
    }
    if (!nextFrameExists_)
    {
        GMX_THROW(
                APIError("There is no next frame, so there should have been no attempt to get it. "
                         "Perhaps the return value of readNextFrame() was misused."));
    }

    // Prepare for reading future frames
    haveProbedForNextFrame_ = false;
    nextFrameExists_        = false;

    // The probe filled trxframeGuard_ with new data, so return it
    return TrajectoryFrame(*trxframeGuard_.get());
}

} // namespace test
} // namespace gmx
