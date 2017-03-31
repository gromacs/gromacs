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

/// \cond internal
/*! \internal \file
 * \brief Declares trajectory and frame structures
 *
 * \ingroup module_python
 */
#ifndef PYGMX_TRAJECTORY_H
#define PYGMX_TRAJECTORY_H

#include <memory>
#include <stdexcept>

#include "gromacs/options/options.h"
#include "gromacs/trajectoryanalysis/analysismodule.h"
#include "gromacs/trajectoryanalysis/runner.h"
#include "python/data.h"

namespace gmx
{
class OptionsVisitor;

namespace pyapi
{
using std::shared_ptr;
using std::unique_ptr;
using std::make_shared;

/// \cond internal
/*! \brief Minimal wrapper for t_trxframe.
 *
 * Hopefully very temporary. Wrapping and exposing t_trxframe is silly.
 * The next step is probably to provide a flexible wrapper to arbitrary
 * TrajectoryAnalysisDataModule data, trajectory or derived.
 * \internal \ingroup module_python
 */
class PyTrajectoryFrame
{
    public:
        /*! \brief Share ownership of a t_trxframe
         *
         * \param frame a managed frame pointer (with an appropriate deleter) that we can safely share.
         *
         * These shared pointers must originate from gmx::trajectory::trxframe_copy
         * where a sensible deleter is provided. Unfortunately, this does not allow
         * the lifetime of a member array to be decoupled from the rest of the frame.
         * \internal
         */
        explicit PyTrajectoryFrame(shared_ptr<t_trxframe> frame);

        /*! \brief Copy a t_trxframe
         *
         * \param frame a t_trxframe object to make a deep copy of
         *
         * The copy is performed by gmx::trajectory::trxframe_copy, which provides
         * a sensible deleter, but cannot allow the lifetime of member arrays to
         * be decoupled from the whole frame.
         * \internal
         */
        explicit PyTrajectoryFrame(const t_trxframe &frame);

        /// With current t_trxframe usage, we have to be careful. \internal
        PyTrajectoryFrame() = delete;
        /// TODO: this should be fine now, but need to test.
        const PyTrajectoryFrame &operator=(const PyTrajectoryFrame &) = delete;

        // TODO: doxygen doesn't seem to like that template...
        /// Return a handle to a buffer (of positions)
        /*! Ideally, this buffer's life is not tied to the frame it is in, but it is
         * while we are using t_trxframe. We can copy the arrays for now to be on
         * the safe side, which happens in the TrajDataArray constructor used.
         * There is also something weird about lifetime management if a unique_ptr
         * is used to provide the buffer interface to Python. Maybe the py::buffer_info
         * does not keep alive the object that provides it, but a quick stepping
         * through the code looks like the object does not survive the conversion
         * from unique_ptr to shared_ptr in the export, so we can just use a shared_ptr
         * return value for now.
         *
         * Requests of certain request type imply a return value type,
         * while the value of the request determines the code path.
         * To specialize for an additional Request type, also specialize
         * t_field<> with a new typedef for the return type.
         *
         * \param t the requested data field. Typing the requests allows
         * optimization and templated behavior. ref trjvectorfield and ref t_field
         * \return shared_ptr to position data of knowable shape
         * \internal
         */
        template<class Request>
        unique_ptr < PyDataHandle < typename t_field<Request>::type > > get_read_handle(const Request &t) const;

        /// A more explicit alternative. Not sure about this one... \internal
        /// \return handle to position data object
        unique_ptr<Data3Handle> get_positions() const;

    private:
        /// Handle to a shareable t_trxframe object \internal
        shared_ptr<t_trxframe>      frame_;
        /// Reduce unnecessary operations and copies by caching results \internal
        mutable shared_ptr< Data3 > position_cache_;
};

/// Declare specialization for Nx3 trajectory data \internal
template<>
unique_ptr<Data3Handle> PyTrajectoryFrame::get_read_handle(const trjvectorfield &t) const;

/// \endcond
}      // end namespace pyapi
}      // end namespace gmx

#endif // header guard
/// \endcond
