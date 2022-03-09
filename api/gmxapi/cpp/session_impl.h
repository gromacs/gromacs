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

#ifndef GMXAPI_SESSION_IMPL_H
#define GMXAPI_SESSION_IMPL_H
/*! \file
 * \brief Declare implementation interface for Session API class(es).
 *
 * \ingroup gmxapi
 */

#include <map>

#include "gromacs/mdrun/runner.h"
#include "gromacs/mdrun/simulationcontext.h"
#include "gromacs/mdrunutility/logging.h"

// Above are the public headers from other modules.
// Following are public headers for the current module.
#include "gmxapi/context.h"
#include "gmxapi/md.h"
#include "gmxapi/md/mdmodule.h"
#include "gmxapi/session/resources.h"
#include "gmxapi/status.h"


namespace gmxapi
{

// Forward declaration
class ContextImpl;   // locally defined in context.cpp
class SignalManager; // defined in mdsignals_impl.h

/*!
 * \brief Implementation class for executing sessions.
 *
 * Since 0.0.3, there is only one context and only one session type. This may
 * change at some point to allow templating on different resource types or
 * implementations provided by different libraries.
 * \ingroup gmxapi
 */
class SessionImpl
{
public:
    //! Use create() factory to get an object.
    SessionImpl() = delete;
    ~SessionImpl();

    /*!
     * \brief Check if the session is (still) running.
     *
     * When a session is launched, it should be returned in an "open" state by the launcher
     * function. \return True if running, false if already closed.
     */
    bool isOpen() const noexcept;

    /*!
     * \brief Explicitly close the session.
     *
     * Sessions should be explicitly `close()`ed to allow for exceptions to be caught by the client
     * and because closing a session involves a more significant state change in the program than
     * implied by a typical destructor. If close() can be shown to be exception safe, this protocol may be removed.
     *
     * \return On closing a session, a status object is transferred to the caller.
     */
    Status close();

    /*!
     * \brief Run the configured workflow to completion or error.
     *
     * \return copy of the resulting status.
     *
     * \internal
     * By the time we get to the run() we shouldn't have any unanticipated exceptions.
     * If there are, they can be incorporated into richer future Status implementations
     * or some other more appropriate output type.
     */
    Status run() noexcept;

    /*!
     * \brief Create a new implementation object and transfer ownership.
     *
     * \param context Shared ownership of a Context implementation instance.
     * \param runnerBuilder MD simulation builder to take ownership of.
     * \param simulationContext Take ownership of the simulation resources.
     * \param logFilehandle Take ownership of filehandle for MD logging
     *
     * \todo Log file management will be updated soon.
     *
     * \return Ownership of new Session implementation instance.
     */
    static std::unique_ptr<SessionImpl> create(std::shared_ptr<ContextImpl> context,
                                               gmx::MdrunnerBuilder&&       runnerBuilder,
                                               gmx::SimulationContext&&     simulationContext,
                                               gmx::LogFilePtr              logFilehandle);

    /*!
     * \brief Add a restraint to the simulation.
     *
     * \param module
     * \return
     */
    Status addRestraint(std::shared_ptr<gmxapi::MDModule> module);

    /*!
     * \brief Get a handle to the resources for the named session operation.
     *
     * \param name unique name of element in workflow
     * \return temporary access to the resources.
     *
     * If called on a non-const Session, creates the resource if it does not yet exist.
     * If called on a const Session,
     * returns nullptr if the resource does not exist.
     */
    gmxapi::SessionResources* getResources(const std::string& name) const noexcept;

    /*!
     * \brief Create SessionResources for a module and bind the module.
     *
     * Adds a new managed resources object to the Session for the uniquely named module.
     * Allows the module to bind to the SignalManager and to the resources object.
     *
     * \param module
     * \return non-owning pointer to created resources or nullptr for error.
     *
     * If the named module is already registered, calling createResources again is considered an
     * error and nullptr is returned.
     */
    gmxapi::SessionResources* createResources(std::shared_ptr<gmxapi::MDModule> module) noexcept;

    /*! \internal
     * \brief API implementation function to retrieve the current runner.
     *
     * \return non-owning pointer to the current runner or nullptr if none.
     */
    gmx::Mdrunner* getRunner();

    /*!
     * \brief Get a non-owning handle to the SignalManager for the active MD runner.
     *
     * Calling code is responsible for ensuring that the SessionImpl is kept alive and "open"
     * while the returned SignalManager handle is in use.
     *
     * \return non-owning pointer if runner and signal manager are active, else nullptr.
     */
    SignalManager* getSignalManager();

    /*!
     * \brief Constructor for use by create()
     *
     * \param context specific context to keep alive during session.
     * \param runnerBuilder ownership of the MdrunnerBuilder object.
     * \param simulationContext take ownership of a SimulationContext
     * \param logFilehandle Take ownership of filehandle for MD logging
     *
     */
    SessionImpl(std::shared_ptr<ContextImpl> context,
                gmx::MdrunnerBuilder&&       runnerBuilder,
                gmx::SimulationContext&&     simulationContext,
                gmx::LogFilePtr              logFilehandle);

private:
    /*!
     * \brief Manage session resources for named workflow elements.
     */
    std::map<std::string, std::unique_ptr<SessionResources>> resources_;

    /*!
     * \brief Extend the life of the owning context.
     *
     * The session will get handles for logging, UI status messages,
     * and other facilities through this interface. This is a facility
     * provided by the client to the Session implementation during
     * Context.launch().
     */
    std::shared_ptr<ContextImpl> context_;

    /*!
     * \brief Simulation runner object.
     *
     * If a simulation Session is active, points to a valid Mdrunner object.
     * Null if simulation is inactive.
     */
    std::unique_ptr<gmx::Mdrunner> runner_;

    /*!
     * \brief An active session owns the resources it is using.
     *
     * This encapsulate details of the run time context that the
     * Session makes available to the simulator, tied to the
     * lifetime of the Session.
     */
    gmx::SimulationContext simulationContext_;

    /*! \brief Handle to file used for logging.
     *
     * \todo Move to RAII filehandle management; open and close in one place.
     */
    gmx::LogFilePtr logFilePtr_;

    /*!
     * \brief Own and manager the signalling pathways for the current session.
     *
     * Registers a stop signal issuer with the stopConditionBuilder that is
     * passed to the Mdrunner at launch. Session members issuing stop signals
     * are proxied through this resource.
     */
    std::unique_ptr<SignalManager> signalManager_;

    /*!
     * \brief Restraints active in this session.
     *
     * Client owns these restraint objects, but session has the ability to
     * lock the resource to take temporary ownership in case the client
     * releases its handle.
     * \todo clarify and update object lifetime management
     * A restraint module manager and / or a mapping of factory functions with
     * which the runner can get objects at run time can encapsulate object management.
     */
    std::map<std::string, std::weak_ptr<gmx::IRestraintPotential>> restraints_;
};

} // end namespace gmxapi

#endif // GMXAPI_SESSION_IMPL_H
