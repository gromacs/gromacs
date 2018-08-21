/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
//
// Created by Eric Irrgang on 7/10/18.
//

#ifndef GMXAPI_MDSIGNALS_IMPL_H
#define GMXAPI_MDSIGNALS_IMPL_H

/*! \file
 * \brief Implementation details for gmxapi::Signal and related gmxapi::md infrastructure.
 *
 * \ingroup gmxapi_md
 */

#include <atomic>

#include "session-impl.h"
#include "gmxapi/session.h"
#include "gmxapi/md/mdsignals.h"
#include "gromacs/compat/make_unique.h"
#include "gromacs/mdlib/simulationsignal.h"
#include "gromacs/mdrun/runner.h"

namespace gmxapi
{


/*!
 * \brief The Signal Implementation interface.
 *
 * A SignalImpl concrete class must implement a `call()` method that issues the signal.
 */
class Signal::SignalImpl
{
    public:
        virtual void call() = 0;
};

/*!
 * \brief Signal implementation for MD simulation stop signals.
 *
 * Provides a call() operator that sets the stop condition for the MD simulation.
 *
 * Client code is not expected to create objects of this type directly, but to retrieve
 * one, wrapped in a gmxapi::Signal for immediate use, from a SignalManager with SignalManager::getSignal()
 *
 * It is the responsibility of the client code to make sure that the address of the Mdrunner
 * remains valid for the lifetime of this object.
 */
class StopSignal : public Signal::SignalImpl
{
    public:
        /*!
         * \brief Create short-lived signal implementation.
         *
         * \param runner non-owning handle
         *
         * The object is constructed with a handle to the runner associated with the SignalManager and
         * owned by the owner of the SignalManager.
         */
        explicit StopSignal(gmx::Mdrunner* runner);

        /*!
         * \brief Set a stop condition for the attached runner.
         */
        void call() override;

    private:
        /// non-owning handle to a runner owned by the owner of the SignalManager.
        gmx::Mdrunner* runner_;
};

/*!
 * \brief Manage signal paths exposed through session resources to gmxapi operations.
 *
 * Manages signals for a single gmx::Mdrunner. Currently only supports a stop signal that
 * is required to be issued by all registered possible issuers before the signal is sent to
 * the associated runner. This is not what we want in the long run. This class should handle
 * signal inputs to operations that take signals as input (like Mdrunner) and should allow
 * multiple subscribers. For additional signal processing, such as boolean operations,
 * additional operations should be inserted in a chain.
 *
 * SignalManager objects are created during Session launch and are owned exclusively by session
 * implementation objects. If Session::isOpen() is true, the SignalManager should still be valid,
 * but the intended use case is for SignalManager handles to be retrieved immediately before use
 * by implementation code within the library with SessionImpl::getSignalManager().
 *
 * A SignalManager should be created for each consumer (each gmx::Mdrunner) in a Session.
 * This occurs in the SessionImpl::create() function.
 *
 * \ingroup gmxapi_md
 */
class SignalManager
{
    public:
        explicit SignalManager(gmx::Mdrunner* runner);
        ~SignalManager();

        /*!
         * \brief Add a name to the list of operations that will be using this signal.
         */
        void addSignaller(std::string name);

        /*!
         * \brief Allow a registered signaller to retrieve a functor.
         *
         * \param name Registered signal issuer.
         * \param signal type of signal the client would like to issue.
         * \return Generic Signal object.
         *
         * \throws gmxapi::ProtocolError if named signaller was not previously registered.
         */
        Signal getSignal(std::string name,
                         md::signals signal);

        /*!
         * \brief A member class that can access SignalManager's private members.
         */
        class LogicalAND;

    private:
        /// Non-owning handle to the associated runner.
        gmx::Mdrunner* runner_;
        /*!
         * \brief Track whether the signal has been issued by each registrant.
         */
        std::map<std::string, std::atomic_bool> called_;
};



}      //end namespace gmxapi

#endif //GMXAPI_MDSIGNALS_IMPL_H
