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

#ifndef GMXAPI_SESSION_RESOURCES_IMPL_H
#define GMXAPI_SESSION_RESOURCES_IMPL_H

/*! \file
 * \brief Implementation details for SessionResources infrastructure.
 *
 * Define the library interface for classes with opaque external interfaces.
 *
 * These are specifically details of the gmx Mdrunner session implementation.
 *
 * \author M. Eric Irrgang <ericirrgang@gmail.com>
 * \ingroup gmxapi
 */

#include <string>

#include "gmxapi/md/mdsignals.h"
#include "gmxapi/session.h"
#include "gmxapi/session/resources.h"

namespace gmxapi
{

/*!
 * \brief Consumer-specific access to Session resources.
 *
 * Each element of work that is managed by a Session and which may need access to Session resources
 * is uniquely identified. SessionResources objects allow client code to be identified by the
 * Session so that appropriate resources can be acquired when needed.
 *
 * Resources are configured at Session launch by SessionImpl::createResources()
 *
 * \ingroup gmxapi
 */
class SessionResources final
{
public:
    /*!
     * \brief Construct a resources object for the named operation.
     *
     * \param session implementation object backing these resources.
     * \param name Unique name of workflow operation.
     */
    SessionResources(SessionImpl* session, std::string name);

    /*!
     * \brief no default constructor.
     *
     * \see SessionResources(SessionImpl* session, std::string name)
     */
    SessionResources() = delete;

    ///@{
    /*!
     * \brief Not moveable or copyable.
     *
     * Objects of this type should only exist in their Session container.
     * If necessary, ownership can be transferred by owning through a unique_ptr handle.
     */
    SessionResources(const SessionResources&) = delete;
    SessionResources& operator=(const SessionResources&) = delete;
    SessionResources(SessionResources&&)                 = delete;
    SessionResources& operator=(SessionResources&&) = delete;
    ///@}

    ~SessionResources();

    /*!
     * \brief Get the name of the gmxapi operation for which these resources exist.
     *
     * \return workflow element name
     */
    std::string name() const;

    /*!
     * \brief Get a Signal instance implementing the requested MD signal.
     *
     * The caller is responsible for ensuring that the session is still active.
     * Unfortunately, there isn't really a way to do that right now. This needs improvemnt
     * in a near future version.
     *
     * Also, this is an external interface that should avoid throwing exceptions for ABI compatibility.
     *
     * \param signal currently must be gmxapi::md::signals::STOP
     * \return callable object.
     *
     * Example:
     *
     *     auto signal = sessionResources->getMdrunnerSignal(md::signals::STOP);
     *     signal();
     *
     * \throws gmxapi::MissingImplementationError if an implementation is not available for the requested signal.
     * \throws gmxapi::ProtocolError if the Session or Signaller is not available.
     */
    Signal getMdrunnerSignal(md::signals signal);

private:
    /*!
     * \brief pointer to the session owning these resources
     */
    SessionImpl* sessionImpl_ = nullptr;

    /*!
     * \brief name of the associated gmxapi operation
     */
    std::string name_;
};

} // end namespace gmxapi

#endif // GMXAPI_SESSION_RESOURCES_IMPL_H
