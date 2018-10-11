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

#include "gmxpre.h"

#include "config.h"

#include "gromacs/compat/make_unique.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/mdlib/sighandler.h"
#include "gromacs/mdrun/logging.h"
#include "gromacs/mdrun/multisim.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/basenetwork.h"
#include "gromacs/utility/init.h"

#include "gmxapi/session.h"
#include "gmxapi/status.h"

#include "session-impl.h"

namespace gmxapi
{

/*!
 * \brief Provide RAII management of communications resource state.
 *
 * To acquire an MpiContextManager is to have assurance that the GROMACS MPI
 * environment is ready to use. When the MpiContextManager is released or
 * goes out of scope, the destructor finalizes the resources.
 *
 * \todo Figure out how to manage MPI versus tMPI.
 * \todo gmx::init should take a subcommunicator rather than use MPI_COMM_WORLD
 * \todo There is no resource for logging or reporting errors during initialization
 * \todo Clarify relationship with gmx::SimulationContext.
 *
 * \ingroup gmxapi
 */
class MpiContextManager
{
    public:
        MpiContextManager()
        {
            gmx::init(nullptr, nullptr);
#ifdef GMX_MPI
#if GMX_MPI
#if GMX_THREAD_MPI
            // With thread-MPI, gmx_mpi_initialized() is false until after
            // spawnThreads in the middle of gmx::Mdrunner::mdrunner(), but
            // without thread-MPI, MPI is initialized by the time gmx::init()
            // returns. In other words, this is not an effective context manager
            // for thread-MPI, but it should be effective for MPI.
            // \todo Distinguish scope / lifetime for comm resources from implementation details.
            // \todo Normalize scope / lifetime of comm resources.
#else
            GMX_ASSERT(gmx_mpi_initialized(), "MPI should be initialized before reaching this point.");
#endif // GMX_THREAD_MPI
#endif // GMX_MPI
#endif // defined(GMX_MPI)
        };

        ~MpiContextManager()
        {
            // This is safe to call. It is a no-op if thread-MPI, and if the
            // constructor completed then MPI is initialized.
            gmx::finalize();
        }

        /*!
         * \brief Exclusive ownership of a scoped context means copying is impossible.
         *
         * \{
         */
        MpiContextManager(const MpiContextManager &)            = delete;
        MpiContextManager &operator=(const MpiContextManager &) = delete;
        //! \}

        /*!
         * \brief Move semantics are trivial.
         *
         * \{
         */
        MpiContextManager(MpiContextManager &&) noexcept            = default;
        MpiContextManager &operator=(MpiContextManager &&) noexcept = default;
        //! \}
};

bool SessionImpl::isOpen() const noexcept
{
    // Currently, an active session is equivalent to an active Mdrunner.
    return bool(runner_);
}

Status SessionImpl::close()
{
    // Assume unsuccessful until proven otherwise.
    auto successful = Status(false);

    // When the Session is closed, we need to know that the MD output has been finalized.
    runner_.reset();
    logFilePtr_.reset();

    if (GMX_LIB_MPI)
    {
        done_commrec(simulationContext_.communicationRecord_);
    }
    // What happens to *cr otherwise?

    done_multisim(multiSim_);

    // \todo provide meaningful result.
    // We should be checking that resources were properly shut down, but
    // there isn't currently a way to do that.
    successful = true;
    return successful;
}

Status SessionImpl::run() noexcept
{
    // Status is failure until proven otherwise.
    auto successful = Status(false);
    GMX_ASSERT(runner_, "SessionImpl invariant implies valid Mdrunner handle.");
    auto rc = runner_->mdrunner();
    if (rc == 0)
    {
        successful = true;
    }
    return successful;
}

std::unique_ptr<SessionImpl> SessionImpl::create(std::shared_ptr<ContextImpl>   context,
                                                 std::unique_ptr<gmx::Mdrunner> runner,
                                                 const gmx::SimulationContext  &simulationContext,
                                                 gmx::LogFilePtr                logFilehandle,
                                                 gmx_multisim_t               * multiSim)
{
    using gmx::compat::make_unique;
    // We should be able to get a communicator (or subcommunicator) through the
    // Context.
    std::unique_ptr<SessionImpl> impl = make_unique<SessionImpl>(std::move(context),
                                                                 std::move(runner),
                                                                 simulationContext,
                                                                 std::move(logFilehandle),
                                                                 multiSim);
    return impl;
}

SessionImpl::SessionImpl(std::shared_ptr<ContextImpl>   context,
                         std::unique_ptr<gmx::Mdrunner> runner,
                         const gmx::SimulationContext  &simulationContext,
                         gmx::LogFilePtr                fplog,
                         gmx_multisim_t               * multiSim) :
    context_(std::move(context)),
    mpiContextManager_(gmx::compat::make_unique<MpiContextManager>()),
    runner_(std::move(runner)),
    simulationContext_(simulationContext),
    logFilePtr_(std::move(fplog)),
    multiSim_(multiSim)
{
    GMX_ASSERT(context_, "SessionImpl invariant implies valid ContextImpl handle.");
    GMX_ASSERT(mpiContextManager_, "SessionImpl invariant implies valid MpiContextManager guard.");
    GMX_ASSERT(runner_, "SessionImpl invariant implies valid Mdrunner handle.");
    GMX_ASSERT(simulationContext_.communicationRecord_, "SessionImpl invariant implies valid commrec.");
    // \todo Session objects can have logic specialized for the runtime environment.

    // For the libgromacs context, a session should explicitly reset global variables that could
    // have been set in a previous simulation during the same process.
    gmx_reset_stop_condition();
}

gmx::Mdrunner *SessionImpl::getRunner()
{
    gmx::Mdrunner * runner = nullptr;
    if (runner_)
    {
        runner = runner_.get();
    }
    return runner;
}

Session::Session(std::unique_ptr<SessionImpl> impl) noexcept
{
    if (impl != nullptr)
    {
        impl_ = std::move(impl);
    }
    GMX_ASSERT(impl_->isOpen(), "SessionImpl invariant implies valid Mdrunner handle.");
}

Status Session::run() noexcept
{
    GMX_ASSERT(impl_, "Session invariant implies valid implementation object handle.");

    const Status status = impl_->run();
    return status;
}

Status Session::close()
{
    GMX_ASSERT(impl_, "Session invariant implies valid implementation object handle.");

    auto status = Status(false);
    if (isOpen())
    {
        // \todo catch exceptions when we know what they might be
        status = impl_->close();
    }

    return status;
}

Session::~Session()
{
    GMX_ASSERT(impl_, "Session invariant implies valid implementation object handle.");
    if (isOpen())
    {
        try
        {
            impl_->close();
        }
        catch (const std::exception &)
        {
            // \todo find some exception-safe things to do with this via the Context interface.
        }
    }
}

bool Session::isOpen() const noexcept
{
    GMX_ASSERT(impl_, "Session invariant implies valid implementation object handle.");
    const auto result = impl_->isOpen();
    return result;
}

//! \cond internal
SessionImpl *Session::getRaw() const noexcept
{
    return impl_.get();
}
//! \endcond

std::shared_ptr<Session> launchSession(Context          * context,
                                       const std::string &filename)
{
    auto session = context->launch(filename);
    return session;
}

SessionImpl::~SessionImpl() = default;

} // end namespace gmxapi
