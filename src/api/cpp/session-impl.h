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

#ifndef GMXAPI_SESSION_IMPL_H
#define GMXAPI_SESSION_IMPL_H
/*! \file
 * \brief Declare implementation interface for Session API class(es).
 *
 * \ingroup gmxapi
 */

#include <map>

#include "gmxapi/context.h"
#include "gmxapi/status.h"

#include "gromacs/mdrun/runner.h"

namespace gmxapi
{

// Forward declaration
class MpiContextManager; // Locally defined in session.cpp
class ContextImpl;       // locally defined in context.cpp

/*!
 * \brief Implementation class for executing sessions.
 *
 * In 0.0.3, there is only one context and only one session type, but this will likely change soon.
 * \ingroup gmxapi
 */
class SessionImpl
{
    public:
        /// Use create() factory to get an object.
        SessionImpl() = delete;
        ~SessionImpl();

        /*!
         * \brief Check if the session is (still) running.
         *
         * When a session is launched, it should be returned in an "open" state by the launcher function.
         * \return True if running, false if already closed.
         */
        bool isOpen() const noexcept;

        /*!
         * \brief Get the current / most recent status.
         *
         * \return Copy of current status object reflecting most recent operation.
         */
        Status status() const noexcept;

        /*!
         * \brief Explicitly close the session.
         *
         * Sessions should be explicitly `close()`ed to allow for exceptions to be caught by the client
         * and because closing a session involves a more significant state change in the program than
         * implied by a typical destructor. If close() can be shown to be exception safe, this protocol may be removed.
         *
         * On closing a session, the status object is transfered to the caller.
         * \return
         */
        std::unique_ptr<Status> close();

        /*!
         * \brief Run the configured workflow to completion or error.
         *
         * \return copy of the resulting status.
         *
         * \internal
         * By the time we get to the run() we shouldn't have any unanticipated exceptions.
         */
        Status run() noexcept;

        /*!
         * \brief Create a new implementation object and transfer ownership.
         *
         * \param context Shared ownership of a Context implementation instance.
         * \param runner MD simulation operation to take ownership of.
         * \return Ownership of new instance.
         *
         * A gmxapi::SignalManager is created with lifetime tied to the SessionImpl. The SignalManager
         * is created with a non-owning pointer to runner. Signal issuers are registered with the
         * manager when createResources() is called.
         */
        static std::unique_ptr<SessionImpl> create(std::shared_ptr<ContextImpl>   context,
                                                   std::unique_ptr<gmx::Mdrunner> runner);

        /*! \internal
         * \brief API implementation function to retrieve the current runner.
         *
         * \return non-owning pointer to the current runner or nullptr if none.
         */
        gmx::Mdrunner* getRunner();

    private:
        /*!
         * \brief Private constructor for use by create()
         *
         * \param context specific context to keep alive during session.
         * \param runner ownership of live Mdrunner object.
         */
        SessionImpl(std::shared_ptr<ContextImpl>   context,
                    std::unique_ptr<gmx::Mdrunner> runner);

        /*!
         * \brief Current / most recent Status for the session.
         *
         * \internal
         * An open session has a valid status object. A closed session has status_ == nullptr.
         */
        std::unique_ptr<Status> status_;

        /*!
         * \brief Extend the life of the owning context. The session will get handles for logging, UI status messages, and other facilities through this interface.
         */
        std::shared_ptr<Context> context_;

        /*!
         * \brief RAII management of gmx::init() and gmx::finalize()
         */
        std::unique_ptr<MpiContextManager> mpiContextManager_;

        std::unique_ptr<gmx::Mdrunner>     runner_;
};

}      //end namespace gmxapi

#endif //GMXAPI_SESSION_IMPL_H
