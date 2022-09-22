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
#ifndef GMXAPI_SYSTEM_H
#define GMXAPI_SYSTEM_H
/*! \file
 * \brief Declare container for molecular systems
 *
 * \ingroup gmxapi
 */
#include <memory>
#include <string>

#include "gmxapi/gmxapi.h"
#include "gmxapi/status.h"

namespace gmxapi
{

// Forward declarations.
class Context; // full specification in gmxapi/context.h
class Session; // full specification in gmxapi/session.h

/*!
 * \brief Deprecated: A wrapper for gmx::Mdrunner
 *
 * Container for molecular model and simulation parameters.
 *
 * \internal
 * It was not intended as such, but gmxapi::System has ended up as basically a
 * wrapper for the gmx::Mdrunner simulator object and serves as the aspect of a
 * gmxapi Session that performs simulation work. As such, this class does not fit the
 * current gmxapi paradigm and will be removed, reworked, or renamed soon.
 *
 * # Protocol
 *
 * As of gmxapi 0.0.6, a simulation is configured and launched as follows.
 *
 * 1. Caller gets a System handle with gmxapi::fromTprFile().
 * 2. Caller optionally attaches additional MD Modules with
 *    getSpec()->addModule(std::shared_ptr<gmxapi::MDModule> module).
 *    See gmxapi::MDHolder
 * 3. Caller gets a runnable object by passing a Context to System::launch()
 *
 * During launch() configured gmxapi::MDModules are attached to the simulator, which is
 * then run by calling run() on the object returned by launch().
 *
 * \ingroup gmxapi
 */
class System final
{
public:
    /*! \brief Private implementation class.
     *
     * System::Impl does not have a public interface and is only exposed in opaque pointers.
     */
    class Impl;

    /*! \brief No copy.
     *
     * The semantics of copying a System are ambiguous, so disallow implicit
     * copy.
     * \{
     * \internal
     * Some sort of prototype or clone idiom could be useful, but
     * needs to explicitly identify any expensive operations.
     */
    System(const System&) = delete;
    System& operator=(const System&) = delete;
    /*! \} */

    /*! \brief Allow move.
     *
     * \{
     */
    System(System&& /*unused*/) noexcept;
    System& operator=(System&& /*unused*/) noexcept;
    /*! \} */

    /*!
     * \brief Create by taking ownership of an implementation object.
     *
     * \param implementation
     */
    explicit System(std::unique_ptr<Impl> implementation);

    /*! \cond internal
     *  Destructor defined later to allow unique_ptr members of partially-defined types.
     */
    ~System();
    /*! \endcond */

    /*!
     * \brief Configure the computing environment for the specified workflow.
     *
     * \return Ownership of a ready-to-run workflow or nullptr if there were errors.
     *
     * If errors occur, they will be stored in the context object. If run without
     * an argument, launch() uses the current context of the System object. If a
     * context argument is given, the system and its configured workflow are
     * translated to the provided context and launched.
     *
     * \param context (optional) execution context in which to launch.
     *
     * \note The Session object does not "own" the Context, but must be able
     * to extend the lifetime of the Context in which it is running.
     *
     * \todo Policy: does System then track the (potentially remote) context or should
     * it be considered to have "forked", and the new session object retrieved from
     * the session handle if needed?
     *
     * \cond internal
     * # Protocol
     *
     * The current implementation of System::launch() performs the following actions.
     *
     * When launch() is called, a new gmxapi::Session is created by passing a
     * gmxapi::Workflow to context->launch().
     * The Workflow basically just contains the TPR filename.
     * 1. A new Mdrunner is created from the information in the gmxapi::Workflow
     * 2. A new Session is created using the ContextImpl and the runner
     *
     * Then, for each module available through getSpec()->getModules(),
     * the session and module are passed to gmxapi::addSessionRestraint().
     * 1. A gmx::IRestraintPotential is retrieved from the module.
     * 2. A unique, named SessionResources is created for the module and attached to the SessionImpl.
     *     1. The module is added as a signaller to the session SignalManager
     *        with getSignalManager()->addSignaller(module->name())
     *     2. The SessionResources is passed to IRestraintPotential::bindSession().
     *        Currently, the only thing the restraint could do at this point is to
     *        save a copy of the pointer and later pass it to gmxapi::getMdrunnerSignal().
     *     3. The restraint is passed to gmx::Mdrunner::addPullPotential(),
     *        which adds the restraint to the global gmx::restraint::Manager,
     *        which then needs to be `clear()`ed after the runner completes.
     *
     * Shared ownership of the Session is returned to the caller of launch().
     *
     * \endcond
     */
    std::shared_ptr<Session> launch(const std::shared_ptr<Context>& context);

    [[nodiscard]] Impl* get() const;

private:
    /*!
     * \brief Opaque pointer to implementation.
     */
    std::unique_ptr<Impl> impl_;
};


/*! \brief Defines an MD workflow from a TPR file.
 *
 * The TPR file has sufficient information to fully specify an MD run, though
 * various parameters are implicit until the work is launched. The TPR filename
 * provided must refer to identical TPR files at the API client and at the
 * main rank of the execution host.
 *
 * \param filename Filesystem path of TPR file.
 * \returns gmxapi::System object with the specified workflow.
 * \ingroup gmxapi
 */
System fromTprFile(const std::string& filename);

class Workflow;
std::shared_ptr<Workflow> getWork(const System::Impl& system);

} // end namespace gmxapi

#endif // include guard
