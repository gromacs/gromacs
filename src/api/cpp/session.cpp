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

#include <cassert>

#include "gmxapi/session.h"
#include "gmxapi/status.h"

#include "session-impl.h"
#include "gromacs/compat/make_unique.h"
#include "gromacs/utility/init.h"

namespace gmxapi
{

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
#else
            assert(gmx_mpi_initialized());
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

        // Disallow copying
        MpiContextManager(const MpiContextManager &)            = delete;
        MpiContextManager &operator=(const MpiContextManager &) = delete;

        // Trivial move
        MpiContextManager(MpiContextManager &&) noexcept            = default;
        MpiContextManager &operator=(MpiContextManager &&) noexcept = default;
};

/*!
 * \brief Check if a an object can be considered "open".
 *
 * This should be generalized to an API idiom.
 *
 * \tparam T type that can be open or closed.
 * \param object something that has a concept of "open" or "closed."
 * \return true if open, false if closed, compiler error if non-sensical.
 */
template<class T>
bool isOpen(const T &object);
//{
//    (void) object;
//    static_assert(false, "Compiler could not find open/close concept for the given object.");
//    return false;
//}

template<>
bool isOpen<SessionImpl>(const SessionImpl &object)
{
    return object.isOpen();
}

template<>
bool isOpen<Session>(const Session &object)
{
    return object.isOpen();
}

bool SessionImpl::isOpen() const noexcept
{
    return status_ != nullptr;
}

Status SessionImpl::status() const noexcept
{
    return *status_;
}

std::unique_ptr<Status> SessionImpl::close()
{
    // When the Session is closed, we need to know that the MD output has been finalized, which currently requires
    // gmx::MDrunner::~MDrunner() to be called.
    runner_.reset();
    std::unique_ptr<Status> status {
        nullptr
    };
    status.swap(status_);
    assert(status_ == nullptr);
    return status;
}

Status SessionImpl::run() noexcept
{
    // Status is failure until proven otherwise.
    Status status {
        false
    };
    assert(runner_ != nullptr);
    auto rc = runner_->mdrunner();
    if (rc == 0)
    {
        status = true;
    }
    return status;
}

std::unique_ptr<SessionImpl> SessionImpl::create(std::shared_ptr<ContextImpl>   context,
                                                 std::unique_ptr<gmx::Mdrunner> runner)
{
    std::unique_ptr<SessionImpl> impl {
        new SessionImpl(std::move(context), std::move(runner))
    };
    return impl;
}

SessionImpl::SessionImpl(std::shared_ptr<ContextImpl>   context,
                         std::unique_ptr<gmx::Mdrunner> runner) :
    status_ {gmx::compat::make_unique<Status>(true)},
context_ {
    std::make_shared<Context>(std::move(context))
},
mpiContextManager_ {
    gmx::compat::make_unique<MpiContextManager>()
},
runner_ {
    std::move(runner)
}
{
    assert(status_ != nullptr);
    assert(context_ != nullptr);
    assert(mpiContextManager_ != nullptr);
    assert(runner_ != nullptr);
}

gmx::Mdrunner *SessionImpl::getRunner()
{
    gmx::Mdrunner * runner {
        nullptr
    };
    if (runner_)
    {
        runner = runner_.get();
    }
    return runner;
}

Session::Session(std::unique_ptr<SessionImpl> &&impl) noexcept :
    impl_ {std::move(impl)}
{
    assert(impl_ != nullptr);
    assert(impl == nullptr);
    assert(impl_->isOpen());
}

Status Session::run() noexcept
{
    assert(impl_ != nullptr);

    Status status {
        impl_->run()
    };
    return status;
}

Status Session::close()
{
    assert(impl_ != nullptr);

    Status status {
        false
    };
    if (isOpen())
    {
        // \todo catch exceptions when we know what they might be
        auto status_ptr = impl_->close();
        if (status_ptr != nullptr)
        {
            status = *status_ptr;
        }
        // what to do if we get nullptr?
    }

    return status;
}

Session::~Session()
{
    assert(impl_ != nullptr);
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
    assert(impl_ != nullptr);
    const auto result = impl_->isOpen();
    return result;
}

SessionImpl *Session::getRaw() const noexcept
{
    return impl_.get();
}

std::shared_ptr<Session> launchSession(Context   * context,
                                       std::string filename)
{
    auto session = context->launch(std::move(filename));
    return session;
}

SessionImpl::~SessionImpl() = default;

} // end namespace gmxapi
