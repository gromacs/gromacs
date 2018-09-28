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
/*! \file
 * \brief Implementation details of gmxapi::Context
 *
 * \author M. Eric Irrgang <ericirrgang@gmail.com>
 * \ingroup gmxapi
 */

#include "gmxapi/context.h"

#include <cassert>
#include <cstring>

#include <memory>
#include <vector>

#include "gromacs/commandline/filenm.h"
#include "gromacs/compat/make_unique.h"
#include "gromacs/mdrun/runner.h"

#include "gmxapi/version.h"

namespace gmxapi
{

/*! \brief Temporary shim until proper messaging.
 *
 */
class warn
{
    public:
        /*! \brief Create a warning message.
         *
         * \param message must outlive warn instance...
         */
        explicit warn(const char* message) :
            message_ {message}
        { };
        //! pointer to string managed somewhere else.
        const char* message_;
};

//! Work-around for GCC bug 58265
constexpr bool BUGFREE_NOEXCEPT_STRING = std::is_nothrow_move_assignable<std::string>::value;

/*!
 * \brief Context implementation.
 *
 * \todo Separate interface and implementation. Context should allow external implementations.
 */
class ContextImpl final : public std::enable_shared_from_this<ContextImpl>
{
    public:
        /*!
         * \brief Default constructor.
         *
         * Don't use this. Use create() to get a shared pointer right away.
         * Otherwise, shared_from_this() is potentially dangerous.
         *
         * \todo Make default constructor private or otherwise reduce brittleness of construction.
         */
        ContextImpl();

        /*!
         * \brief Factory function
         *
         * Since this class provides `shared_from_this`, we need to make sure
         * that it never exists without a shared_ptr owning it.
         *
         * If we can confirm `shared_from_this` is no longer necessary, implementation may change.
         *
         * \return ownership of a new object
         */
        static std::shared_ptr<gmxapi::ContextImpl> create();

        /*!
         * \brief Copy disallowed because Session state would become ambiguous.
         *
         * The API implementation needs to unambiguously determine
         * which Sessions and Contexts are associated with each other.
         * \{
         */
        ContextImpl(const ContextImpl&) = delete;
        ContextImpl                  &operator=(const ContextImpl &) = delete;
        //! \}

        /*!
         * \brief Objects are trivial to move.
         *
         * \{
         */
        ContextImpl(ContextImpl &&) noexcept(BUGFREE_NOEXCEPT_STRING)            = default;
        ContextImpl &operator=(ContextImpl &&) noexcept(BUGFREE_NOEXCEPT_STRING) = default;
        //! \}        static std::shared_ptr<gmxapi::ContextImpl> create();

        /*!
         * \brief Retain the ability to find a launched session while it exists.
         *
         * The client owns the Session launched by a Context, but it is helpful
         * for the Context to know if it has an active Session associated with it.
         */
        std::weak_ptr<Session>  session_;

        /*!
         * \brief mdrun command line arguments.
         *
         * Store arguments provided by the client and pass them when launching
         * a simulation runner. This allows client code to access the same
         * options as are available to mdrun on the command line while the API
         * evolves.
         */
        MDArgs                  mdArgs_;
};

ContextImpl::ContextImpl()
{
    assert(session_.expired());
}

std::shared_ptr<gmxapi::ContextImpl> ContextImpl::create()
{
    std::shared_ptr<ContextImpl> impl = std::make_shared<ContextImpl>();
    return impl;
}

// As of gmxapi 0.0.3 there is only one Context type
Context::Context() :
    Context {ContextImpl::create()}
{
    assert(impl_ != nullptr);
}

Context::Context(std::shared_ptr<ContextImpl> impl) :
    impl_ {std::move(impl)}
{
    assert(impl_ != nullptr);
}

void Context::setMDArgs(const MDArgs &mdArgs)
{
    impl_->mdArgs_ = mdArgs;
}

Context::~Context() = default;

} // end namespace gmxapi
