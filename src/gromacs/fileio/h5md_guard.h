/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2024- The GROMACS Authors
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

/*! \brief Helper functions for managing scope guards in the H5md module.
 *
 * \author Petter Johansson <pettjoha@kth.se>
 */

#ifndef GMX_FILEIO_H5MD_GUARD_H
#define GMX_FILEIO_H5MD_GUARD_H

#include "gmxpre.h"

#include "config.h"

#include <utility>

#include "external/scope_guard/scope_guard.h"

#include "h5md.h"

#if GMX_USE_HDF5

#    include <hdf5.h>

namespace gmx
{

//! \brief Function signature for HDF5 closing functions (H5Gclose, etc.)
using H5mdCloseType = herr_t (*)(hid_t);

/*! \brief Bound closing function created by scope guard helper functions.
 *
 * Required since we cannot return anonymous types from helper functions,
 * they must be bound to a concrete type.
 */
class H5mdCloser
{
public:
    H5mdCloser(hid_t object, H5mdCloseType deleter) :
        object_{ std::move(object) }, deleter_{ std::move(deleter) }
    {
    }
    void operator()() { deleter_(object_); }

private:
    //! Object to be deleted by the callback function
    hid_t object_;
    //! Callback function
    H5mdCloseType deleter_;
};

//! \brief Alias for scope guards used in the H5md module.
using H5mdGuard = sg::detail::scope_guard<H5mdCloser>;

/*! \brief Helper function for guarding an HDF5 object with a closing function (H5Gclose, etc.).
 *
 *  When the returned scope guard goes out of scope the closing function is
 *  called with the object identifier as its argument, ensuring that the resource
 *  is freed. If the identifier needs to be stored for later use it should not
 *  be guarded using these functions.
 *
 *  Intended use is for the object to be the result of an opening or creating
 *  function call, which is why it is returned along with the guard. E.g. to
 *  open a group:
 *
 *    auto [object, guard] = makeH5mdGuard(openOrCreateGroup(<...>), H5Gclose);
 *
 *  \param[in] object Object to guard.
 *  \param[in] unboundDeleter
 *  \returns Pair of the (object, guard), where object is the id and guard
 *  is the scope guard which closes the object once it exits scope.
 */
std::pair<hid_t, H5mdGuard> makeH5mdGuard(hid_t object, H5mdCloseType unboundCloser);

/*! \brief Helper function for guarding an HDF5 group with H5Gclose.
 *
 *  \param[in] object Object to guard.
 *  \returns Pair of the (object, guard), where object is the id and guard
 *  is the scope guard which closes the object once it exits scope.
 */
std::pair<hid_t, H5mdGuard> makeH5mdGroupGuard(hid_t object);

/*! \brief Helper function for guarding an HDF5 data set with H5Dclose.
 *
 *  \param[in] object Object to guard.
 *  \returns Pair of the (object, guard), where object is the id and guard
 *  is the scope guard which closes the object once it exits scope.
 */
std::pair<hid_t, H5mdGuard> makeH5mdDataSetGuard(hid_t object);

/*! \brief Helper function for guarding an HDF5 data space with H5Sclose.
 *
 *  \param[in] object Object to guard.
 *  \returns Pair of the (object, guard), where object is the id and guard
 *  is the scope guard which closes the object once it exits scope.
 */
std::pair<hid_t, H5mdGuard> makeH5mdDataSpaceGuard(hid_t object);

/*! \brief Helper function for guarding an HDF5 property list with H5Pclose.
 *
 *  \param[in] object Object to guard.
 *  \returns Pair of the (object, guard), where object is the id and guard
 *  is the scope guard which closes the object once it exits scope.
 */
std::pair<hid_t, H5mdGuard> makeH5mdPropertyListGuard(hid_t object);

/*! \brief Helper function for guarding an HDF5 attribute with H5Aclose.
 *
 *  \param[in] object Object to guard.
 *  \returns Pair of the (object, guard), where object is the id and guard
 *  is the scope guard which closes the object once it exits scope.
 */
std::pair<hid_t, H5mdGuard> makeH5mdAttributeGuard(hid_t object);

/*! \brief Helper function for guarding an HDF5 type with H5Tclose.
 *
 *  \param[in] object Object to guard.
 *  \returns Pair of the (object, guard), where object is the id and guard
 *  is the scope guard which closes the object once it exits scope.
 */
std::pair<hid_t, H5mdGuard> makeH5mdTypeGuard(hid_t object);

} // namespace gmx

#endif // GMX_USE_HDF5
#endif // GMX_FILEIO_H5MD_GUARD_H
