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
/*! \file
 * \brief Launch and manage a GROMACS session.
 *
 * \author M. Eric Irrgang <ericirrgang@gmail.com>
 *
 * \ingroup gmxapi
 */

#ifndef GROMACS_SESSION_H
#define GROMACS_SESSION_H

/*! \file
 * \brief Declarations for a workflow execution session and helper functions.
 *
 * \ingroup gmxapi
 */

#include <memory>
#include <string>

namespace gmxapi
{

// forward declarations
class Context; // defined in gmxapi/context.h
class MDModule;
class Status;   // defined in gmxapi/status.h
class Workflow; // implementation detail

/*!
 * \brief Private implementation class for a session.
 *
 * Actual implementation class may depend on the execution context, but this should be irrelevant to
 * a client. The implementation details are not exposed in the high-level API, but may be in the
 * extension API. See developer documentation for details.
 *
 * \ingroup gmxapi
 */
class SessionImpl;

/*!
 * \brief Workflow execution session.
 *
 * When a workflow is launched in an execution context, the result is a Session object
 * that serves as a handle to interact with the running workflow. The handle allows dynamic changes
 * to the workflow, or control or data crossing the API boundary.
 *
 * Separating run() from construction allows the client to examine the running execution
 * environment or to retrieve the communicator before beginning long-running computation.
 *
 * The session should be explicitly `close()`ed before destruction to allow exception
 * handling during shutdown. The destructor will close() the session if it is not
 * closed already, but errors may be hidden.
 *
 * The session owns a Status object that can be shared and retained by client code
 * for further inspection.
 *
 * \ingroup gmxapi
 */
class Session
{
public:
    /// A session must be created by launching a workflow in an execution context.
    Session() = delete;

    /*! \brief A session cannot be copied, only moved.
     *
     * For shared ownership of a session, use a shared pointer.
     * \{
     */
    Session(const Session&) = delete;
    Session& operator=(const Session&) = delete;
    //! \}

    /*!
     * \brief Pass ownership of a Session.
     *
     * \{
     */
    Session(Session&&) noexcept = default;
    Session& operator=(Session&&) noexcept = default;
    //! \}

    /*!
     * \brief Construct by taking ownership of an implementation object.
     *
     * \param impl Concrete object to take ownership of.
     */
    explicit Session(std::unique_ptr<SessionImpl> impl) noexcept;

    /*!
     * \brief Destroy Session.
     *
     * If the session is still active (has not been closed) then it will be closed
     * with exceptions suppressed. If possible, problems encountered will be
     * noted in the Status object, which the client may have retained shared
     * ownership of.
     */
    ~Session();

    /*!
     * \brief Close a running session.
     *
     * close() should be called before destroying the Session object so that the
     * client can catch any exceptions thrown during shut down that may be
     * unavoidable in the parallel computing environment.
     *
     * \return status of close() operation.
     */
    Status close();

    /*!
     * \brief Run the current workflow to completion.
     *
     * The client should examine the Status object for errors and resulting workflow state.
     * \return the Session's Status after the run.
     */
    Status run() noexcept;

    /*!
     * \brief Check if session is running.
     *
     * A Session will be "open" from the point of launch until "close" is
     * called. If this is not true, either an error has occurred or there is
     * a bug in the implementation.
     *
     * \return `true` if the session is running, else `false`
     */
    bool isOpen() const noexcept;

    /*! \cond internal
     * \brief Get a non-owning handle to the implementation object.
     *
     * Get a raw pointer to the implementation object. The pointer is valid only during the lifetime of the Session,
     * so retain a shared pointer to this Session object or only hold the pointer for the duration of a code block
     * guaranteed to exist entirely within the lifetime of a Session object.
     *
     * \return opaque pointer used by gmxapi implementation and extension code.
     */
    SessionImpl* getRaw() const noexcept;
    //! \endcond

private:
    //! \brief opaque pointer to implementation
    std::unique_ptr<SessionImpl> impl_;
};

/*!
 * \brief Add a uniquely identifiable restraint instance to the MD simulator.
 *
 * \param session Session with an MD runner to attach a restraint to.
 * \param restraint wrapped restraint force calculation object.
 * \return success if restraint was attached successfully, else failure.
 *
 * \todo Clarify whether we are adding modules generally, adding restraint modules, or adding restraints.
 * \todo Update for new scheme in which all restraints are managed by a single Restraint module.
 * \todo Figure out what this object should return to provide more utility and rigorous error checking.
 */
Status addSessionRestraint(Session* session, std::shared_ptr<gmxapi::MDModule> restraint);

/*!
 * \brief Launch a workflow in the provided execution context.
 *
 * \param context Execution environment
 * \param work Directed acyclic graph defining workflow and data flow.
 * \return non-unique ownership of the session or nullptr in failure.
 *
 * The provided context maintains a weak reference to the executing session, while
 * the session extends the life of the context.
 *
 * If session fails to launch (nullptr returned) examine the Status of the context
 * for details.
 * \ingroup gmxapi
 */
std::shared_ptr<Session> launchSession(Context* context, const Workflow& work) noexcept;


} // end namespace gmxapi

#endif // GROMACS_SESSION_H
