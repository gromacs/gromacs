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

#ifndef GMXAPI_MDSIGNALS_IMPL_H
#define GMXAPI_MDSIGNALS_IMPL_H

/*! \file
 * \brief Implementation details for gmxapi::Signal and related gmxapi::md infrastructure.
 *
 * \ingroup gmxapi_md
 */

#include <atomic>
#include <functional>
#include <map>
#include <memory>
#include <string>

#include "gromacs/mdlib/simulationsignal.h"
#include "gromacs/mdlib/stophandler.h"
#include "gromacs/mdrun/runner.h"

// Public gmxapi headers.
#include "gmxapi/md/mdsignals.h"
#include "gmxapi/session.h"

// Internal gmxapi headers.
#include "session_impl.h"

namespace gmxapi
{


/*! \internal
 * \brief The Signal Implementation interface.
 *
 * A SignalImpl concrete class must implement a `call()` method that issues the signal.
 *
 * \ingroup gmxapi_md
 */
class Signal::SignalImpl
{
public:
    //! Required functor behavior.
    virtual void call() = 0;

    //! May be subclassed.
    virtual ~SignalImpl() = default;
};

/*!
 * \brief Manage signal paths exposed through session resources to gmxapi operations.
 *
 * Manages signals for a single gmx::Mdrunner. Currently only supports a stop signal that
 * is required to be issued by all registered possible issuers before the signal is sent to
 * the associated runner. This is not what we want in the long run.
 * \todo This class should handle signal inputs to operations that take signals as input (like
 * Mdrunner) and \todo should allow multiple subscribers. For additional signal processing, such as
 * boolean operations, additional operations should be inserted in a chain.
 *
 * SignalManager objects are created during Session launch and are owned exclusively by session
 * implementation objects. If Session::isOpen() is true, the SignalManager should still be valid,
 * but the intended use case is that SignalManager handles should be retrieved immediately before
 * use by implementation code within the library with SessionImpl::getSignalManager().
 *
 * A SignalManager should be created for each signal consumer (each gmx::Mdrunner) in a Session.
 * This occurs in the SessionImpl::create() function.
 *
 * \ingroup gmxapi_md
 */
class SignalManager
{
public:
    /*!
     * \brief Set up a manager to mediate access to an upcoming MD stop handler.
     *
     * \param mdStopHandlerBuilder access to a builder that can be used during construction.
     */
    explicit SignalManager(gmx::StopHandlerBuilder* mdStopHandlerBuilder);

    //! \cond
    ~SignalManager();
    //! \endcond

    /*!
     * \brief Add a name to the list of operations that will be using this signal.
     */
    void addSignaller(const std::string& name);

    /*!
     * \brief Allow a registered signaller to retrieve a functor.
     *
     * \param name Registered signal issuer.
     * \param signal type of signal the client would like to issue.
     * \return Generic Signal object.
     *
     * \throws gmxapi::ProtocolError if named signaller was not previously registered.
     */
    Signal getSignal(const std::string& name, md::signals signal);

    /*!
     * \brief Signal operation that issues only when all sources have issued.
     *
     * Implemented as a member class that can access SignalManager's private members.
     * \todo Decouple logical operations from SignalManager class definition.
     */
    class LogicalAND;

private:
    //! Non-owning handle to the associated runner.
    gmx::Mdrunner* runner_;


    /*!
     * \brief State of the stop condition to be returned by the registered MD signaller.
     *
     * Ownership is shared by the function objects in the StopConditionHandler
     * (owned by the simulator), which read the value, and the
     * SessionImpl SignalManager, which mediates write access.
     *
     * The signal state is either gmx::StopSignal::noSignal or gmx::StopSignal::stopAtNextNSStep,
     * so atomicity is not important, and we share the state across
     * threads in a tMPI simulation.
     */
    std::shared_ptr<gmx::StopSignal> state_;

    /*!
     * \brief Track whether the signal has been issued by each registrant.
     *
     * \todo This is an implementation detail of LogicalAND that should not be here.
     */
    std::map<std::string, std::atomic_bool> called_;
};


} // end namespace gmxapi

#endif // GMXAPI_MDSIGNALS_IMPL_H
