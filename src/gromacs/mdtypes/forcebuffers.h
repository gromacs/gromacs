/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2020, by the GROMACS development team, led by
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

/*! \libinternal \file
 *
 * \brief
 * This file contains the definition of a container for force buffers.
 *
 * \author Berk Hess
 *
 * \inlibraryapi
 * \ingroup module_mdtypes
 */

#ifndef GMX_MDTYPES_FORCEBUFFERS_H
#define GMX_MDTYPES_FORCEBUFFERS_H

#include "gromacs/gpu_utils/hostallocator.h"
#include "gromacs/math/arrayrefwithpadding.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/classhelpers.h"

namespace gmx
{

enum class PinningPolicy : int;

/*! \libinternal \brief A view of the force buffer
 *
 * This is used for read and/or write access to the force buffer.
 */
class ForceBuffersView
{
public:
    //! Constructor, creates a view to \p force
    ForceBuffersView(const ArrayRefWithPadding<RVec>& force) : force_(force) {}

    //! Copy constructor deleted to avoid creating non-const from const
    ForceBuffersView(const ForceBuffersView& o) = delete;

    //! Move constructor deleted, but could be implemented
    ForceBuffersView(ForceBuffersView&& o) = delete;

    ~ForceBuffersView() = default;

    //! Copy assignment deleted to avoid creating non-const from const
    ForceBuffersView& operator=(const ForceBuffersView& v) = delete;

    //! Move assignment, moves the view
    ForceBuffersView& operator=(ForceBuffersView&& v) = default;

    //! Returns a const arrayref to the force buffer without padding
    ArrayRef<const RVec> force() const { return force_.unpaddedConstArrayRef(); }

    //! Returns an arrayref to the force buffer without padding
    ArrayRef<RVec> force() { return force_.unpaddedArrayRef(); }

    //! Returns an ArrayRefWithPadding to the force buffer
    ArrayRefWithPadding<RVec> forceWithPadding() { return force_; }

private:
    //! The force buffer
    ArrayRefWithPadding<RVec> force_;
};

/*! \libinternal \brief Object that holds the force buffers
 *
 * More buffers can be added when needed. Those should also be added
 * to ForceBuffersView.
 * Can be pinned for efficient transfer to/from GPUs.
 * All access happens through the ForceBuffersView object.
 */
class ForceBuffers
{
public:
    //! Constructor, creates an empty force buffer with pinning not active
    ForceBuffers(PinningPolicy pinningPolicy = PinningPolicy::CannotBePinned);

    //! Copy constructor deleted, but could be implemented
    ForceBuffers(const ForceBuffers& o) = delete;

    //! Move constructor deleted, but could be implemented
    ForceBuffers(ForceBuffers&& o) = delete;

    ~ForceBuffers();

    //! Copy assignment operator, sets the pinning policy to CannotBePinned
    ForceBuffers& operator=(ForceBuffers const& o);

    //! Move assignment operator, deleted but could be implemented
    ForceBuffers& operator=(ForceBuffers&& o) = delete;

    //! Returns a const view to the force buffers
    const ForceBuffersView& view() const { return view_; }

    //! Returns a view to the force buffer
    ForceBuffersView& view() { return view_; }

    //! Resize the force buffer
    void resize(int numAtoms);

    /*! \brief Returns the active pinning policy.
     *
     * Does not throw.
     */
    PinningPolicy pinningPolicy() const;

private:
    //! The force buffer
    PaddedHostVector<RVec> force_;
    //! The view to the force buffer
    ForceBuffersView view_;
};

} // namespace gmx

#endif
