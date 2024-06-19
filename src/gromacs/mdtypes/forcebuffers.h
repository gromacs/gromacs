/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
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

#include <memory>
#include <utility>

#include "gromacs/gpu_utils/hostallocator.h"
#include "gromacs/math/arrayrefwithpadding.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/gmxassert.h"

namespace gmx
{
template<typename T>
class ArrayRef;
template<typename T>
class ArrayRefWithPadding;

enum class PinningPolicy : int;

/*! \libinternal \brief A view of the force buffer
 *
 * This is used for read and/or write access to the force buffer.
 */
class ForceBuffersView
{
public:
    //! Constructor, creates a view to \p force
    ForceBuffersView(const ArrayRefWithPadding<RVec>& force,
                     const ArrayRefWithPadding<RVec>& forceMtsCombined,
                     const bool                       useForceMtsCombined) :
        force_(force), forceMtsCombined_(forceMtsCombined), useForceMtsCombined_(useForceMtsCombined)
    {
    }

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

    //! Returns a const arrayref to the MTS force buffer without padding
    ArrayRef<const RVec> forceMtsCombined() const
    {
        GMX_ASSERT(useForceMtsCombined_, "Need the MTS buffer");
        return forceMtsCombined_.unpaddedConstArrayRef();
    }

    //! Returns an arrayref to the MTS force buffer without padding
    ArrayRef<RVec> forceMtsCombined()
    {
        GMX_ASSERT(useForceMtsCombined_, "Need the MTS buffer");
        return forceMtsCombined_.unpaddedArrayRef();
    }

    //! Returns an ArrayRefWithPadding to the MTS force buffer
    ArrayRefWithPadding<RVec> forceMtsCombinedWithPadding()
    {
        GMX_ASSERT(useForceMtsCombined_, "Need the MTS buffer");
        return forceMtsCombined_;
    }

private:
    //! The force buffer
    ArrayRefWithPadding<RVec> force_;
    //! The force buffer for combined fast and slow forces with MTS
    ArrayRefWithPadding<RVec> forceMtsCombined_;
    // GCC 9 complains about unused attribute "unused" as it never warns about unused members,
    // while clang requires it to avoid -Wunused
    GCC_DIAGNOSTIC_IGNORE("-Wattributes")
    //! Whether we use forceMtsCombined_
    gmx_used_in_debug bool useForceMtsCombined_;
    GCC_DIAGNOSTIC_RESET
};

/*! \libinternal \brief Object that holds the force buffers
 *
 * Contains a normal force buffer and optionally a force buffer for combined fast and slow
 * forces for use with multiple time stepping.
 * More buffers can be added when needed. Those should also be added
 * to ForceBuffersView.
 * The force buffer (not forceMtsCombined) can be pinned for efficient transfer to/from GPUs.
 * All access happens through the ForceBuffersView object.
 */
class ForceBuffers
{
public:
    //! Constructor, creates an empty force buffer with pinning not active and no MTS force buffer
    ForceBuffers();

    /*! \brief Constructor, with options for using the MTS force buffer and the pinning policy
     *
     * \param[in] useForceMtsCombined  Whether to enable use of the forceMtsCombined buffer
     * \param[in] pinningPolicy        The pinning policy for the force (not MTS) buffer
     */
    ForceBuffers(bool useForceMtsCombined, PinningPolicy pinningPolicy);

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
    //! Force buffer with combined fast and slow forces for use with multiple time stepping
    PaddedHostVector<RVec> forceMtsCombined_;
    //! The view to the force buffer
    ForceBuffersView view_;
    //! Whether we use forceMtsCombined_
    bool useForceMtsCombined_;
};

} // namespace gmx

#endif
