/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2018,2019, by the GROMACS development team, led by
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
 * \brief Declares interface of a class for tests that want to inspect
 * trajectories produced by mdrun.
 *
 * Intended usage is to create a TrajectoryFrameReader. Successive
 * calls to its readNextFrameStub and frame methods will return a handle
 * to a t_trxframe object for each frame.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdrun_integration_tests
 */
#ifndef GMX_PROGRAMS_MDRUN_TESTS_TRAJECTORYREADER_H
#define GMX_PROGRAMS_MDRUN_TESTS_TRAJECTORYREADER_H

#include <memory>
#include <string>

#include "gromacs/fileio/oenv.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/unique_cptr.h"

namespace gmx
{

class TrajectoryFrame;

namespace test
{

//! Convenience smart pointer typedef
typedef unique_cptr<gmx_output_env_t, output_env_done> oenv_ptr;
//! Convenience smart pointer typedef
typedef unique_cptr<t_trxstatus, close_trx> trxstatus_file_ptr;
//! Helper function to free all resources
void done_trxframe(t_trxframe* fr);
//! Convenience smart pointer typedef
typedef unique_cptr<t_trxframe, done_trxframe> trxframe_ptr;

/*! \internal
 * \brief Manages returning a t_trxframe whose contents were read from
 * successive frames of an trajectory file. */
class TrajectoryFrameReader
{
public:
    /*! \brief Attempt to read the next frame from the trajectory file.
     *
     * \return Whether a frame was available to read.
     *
     * This call wraps the read_first_frame()/read_next_frame()
     * API, which does the file opening as a side effect of
     * reading the first frame.
     *
     * If true is returned, then TrajectoryFrame frame() should be called
     * to get access to the data. If false is returned, then no
     * further data exists and no further call to
     * readNextFrameStub or TrajectoryFrame frame() should occur.
     *
     * \throws FileIOError upon reading the first frame, if the trajectory file cannot be opened
     * \throws APIError    if an earlier probe has not been properly handled
     *                     (by calling frame(), or stopping trying to read
     *                     from the file). */
    bool readNextFrame();
    /*! \brief Return the next frame from the trajectory file.
     *
     * If the next frame has not been probed for, then probe for
     * it. If no next frame exists, then throw APIError, because
     * user code should have called readNextFrameStub itself if this
     * is possible. (This permits user code to avoid making calls
     * to readNextFrameStub in a case where it already knows that
     * the frame exists.)
     *
     * \throws APIError  if no next frame exists, or if it lacks either time or step number. */
    TrajectoryFrame frame();
    /*! \brief Constructor
     *
     * \param[in] filename  Name of trajectory file to open and read. */
    explicit TrajectoryFrameReader(const std::string& filename);

private:
    //! Name of trajectory file to open and read
    std::string filename_;
    //! Owning handle of output environment object
    oenv_ptr oenvGuard_;
    //! Owning handle of an open trajectory file ready to read frames.
    trxstatus_file_ptr trajectoryFileGuard_;
    //! Owning handle of contents of trajectory file frame after reading.
    const trxframe_ptr trxframeGuard_;
    //! Whether the first frame has been read
    bool haveReadFirstFrame_;
    //! Whether the API has been used properly (ie. probe before reading).
    bool haveProbedForNextFrame_;
    //! Whether there has been a probe that found a next frame.
    bool nextFrameExists_;

    // Multiple owners of these resources isn't very sensible, so prevent it
    GMX_DISALLOW_COPY_AND_ASSIGN(TrajectoryFrameReader);
};

} // namespace test
} // namespace gmx

#endif
