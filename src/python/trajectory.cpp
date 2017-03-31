/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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

#include "gmxpre.h"

#include "trajectory.h"

#include "gromacs/trajectory/trajectoryframe.h"

namespace gmx
{
namespace pyapi
{

using std::shared_ptr;
using std::unique_ptr;
using std::make_shared;

PyTrajectoryFrame::PyTrajectoryFrame(shared_ptr<t_trxframe> frame) :
    frame_ {frame}
{
}

PyTrajectoryFrame::PyTrajectoryFrame(const t_trxframe &frame) :
    frame_ {gmx::trajectory::trxframe_copy(frame)}
{
}


// Implementation to retrieve a read handle
template<>
unique_ptr<Data3Handle> PyTrajectoryFrame::get_read_handle<trjvectorfield>(const trjvectorfield &t) const
{
    // return value
    unique_ptr<Data3Handle> handle;

    if (t == trjvectorfield::POSITION)
    {
        // Currently, individual arrays are not separable from the frame object without copy.
        if (position_cache_)
        {
            // if we've already got a copy floating around, use it
        }
        else
        {
            // stash the data safely
            position_cache_ = make_shared< Data3 >(frame_->x[0], frame_->natoms);
        }
        handle.reset(new LocalTrajDataHandle(position_cache_));
    }
    else
    {
        //TODO: a switch expression may be fine...
        return nullptr;
    }
    return handle;
}


unique_ptr<Data3Handle> PyTrajectoryFrame::get_positions() const
{
    return get_read_handle<trjvectorfield>(trjvectorfield::POSITION);
}


} // end namespace pyapi
} // end namespace gmx
