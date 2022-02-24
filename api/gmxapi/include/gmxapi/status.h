/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
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

#ifndef GMXAPI_STATUS_H
#define GMXAPI_STATUS_H
/*! \file
 * \brief Declare Status class for API operations with no other results.
 *
 * \author M. Eric Irrgang <ericirrgang@gmail.com>
 * \ingroup gmxapi
 */

#include <memory>

namespace gmxapi
{
/*! \brief Trivial API operation return value.
 *
 * Returned by some API operations when a richer return value type does not
 * exist.
 *
 * \cond internal
 * In general, API operations should return an object and not to allow
 * unsuccessful operations to go undetected. In some cases, there is either no
 * obvious result object or the API needs a placeholder for a more elaborate
 * result to be implemented in a future version. This class allows a return
 * value and potential proxy object as a fall-back where no other option exists.
 * \endcond
 *
 * \ingroup gmxapi
 */
class Status final
{
public:
    /*!
     * \brief Status is "unsuccessful" until set otherwise.
     *
     * \internal
     * Default constructor can be used for convenience when preparing the
     * return value for an operation that should be assumed unsuccessful
     * until proven otherwise.
     */
    Status();
    /*!
     * \brief Can be copy-initialized
     *
     * \param status
     */
    Status(const Status& status);
    /*!
     * \brief Can be moved.
     */
    Status(Status&& /*unused*/) noexcept;
    /*!
     * \brief Can be copy-assigned.
     *
     * \param status
     * \return reference to lhs.
     */
    Status& operator=(const Status& status);
    /*!
     * \brief Transfer ownership by move assignment.
     *
     * \param status
     * \return reference to lhs
     */
    Status& operator=(Status&& status) noexcept;
    /*!
     * \brief Set success status from boolean.
     *
     * \param success true to indicate successful operation.
     * \return reference to lhs
     */
    Status& operator=(bool success);

    /*!
     * \brief Initialize with success set true or false.
     *
     * \param success
     */
    explicit Status(bool success);

    /*!
     * \brief Clean up resources.
     *
     * Note non-virtual destructor. This class is not heritable.
     */
    ~Status();

    /*!
     * \brief Check success status.
     *
     * Forces evaluation of any pending computation. The operation that
     * returned the Status object is guaranteed to complete, successfully
     * or unsuccessfully, before this function returns.
     *
     * \return true if the operation described was successful.
     *
     * \internal
     * Unsuccessful operations should have more useful information associated
     * than just a boolean. Future versions should behave more like a
     * `gmx::future<gmx::expected<T>>`, but this functionality is not yet
     * available.
     */
    bool success() const;

private:
    //! \cond
    class Impl;
    std::unique_ptr<Impl> impl_;
    //! \endcond
};

} // end namespace gmxapi

#endif // GMXAPI_STATUS_H
