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
 */

#include "gmxpre.h"

#include <cassert>
#include <cstring>

#include <memory>
#include <utility>
#include <vector>

#include "gmxapi/version.h"

#include "gmxapi/context.h"
#include "gmxapi/exceptions.h"
#include "gmxapi/gmxapi.h"
#include "gmxapi/status.h"
#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/compat/make_unique.h"
#include "gromacs/domdec/domdec.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/mdlib/main.h"
#include "gromacs/mdlib/mdrun.h"
#include "gromacs/mdlib/repl_ex.h"
#include "gromacs/mdrun/runner.h"
#include "gromacs/mdrunutility/handlerestart.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

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
        /// pointer to string managed somewhere else.
        const char* message_;
};

using gmxapi::MDArgs;

/*!
 * \brief Context implementation base class.
 *
 * Execution contexts have a uniform interface specified by the API. Implementations for
 * particular execution environments can specialize / derive from this base.
 *
 * \todo Separate interface and implementation.
 */
class ContextImpl final : public std::enable_shared_from_this<ContextImpl>
{
    public:
        static std::shared_ptr<gmxapi::ContextImpl> create();

        /*!
         * \brief Default constructor.
         *
         * Don't use this. Use create() to get a shared pointer right away.
         * Otherwise, shared_from_this() is potentially dangerous.
         */
        ContextImpl();

        /*!
         * \brief Get a reference to the current status object.
         *
         * \return shared ownership of the current status.
         */
        std::shared_ptr<const Status> status() const noexcept;

        /*!
         * \brief Status of the last operation in the local context.
         *
         * This pointer should always be valid while the context handle exists and
         * client code can extend the life of the object. We use a shared_ptr because
         * it may be expensive or dangerous to copy the object when it is most needed.
         */
        std::shared_ptr<Status> status_;
        std::weak_ptr<Session>  session_;

        MDArgs                  mdArgs_;
};

ContextImpl::ContextImpl() :
    status_ {std::make_shared<Status>(true)},
session_ {}
{
    assert(status_->success());
    assert(session_.expired());
}

std::shared_ptr<gmxapi::ContextImpl> ContextImpl::create()
{
    auto impl = std::make_shared<gmxapi::ContextImpl>();
    return impl;
}

std::shared_ptr<const Status> ContextImpl::status() const noexcept
{
    return status_;
}

// In 0.0.3 there is only one Context type
Context::Context() :
    Context {std::make_shared<ContextImpl>()}
{
    assert(impl_ != nullptr);
}

Context::Context(std::shared_ptr<ContextImpl> &&impl) :
    impl_ {std::move(impl)}
{
    assert(impl_ != nullptr);
}

void Context::setMDArgs(const MDArgs &mdArgs)
{
    impl_->mdArgs_ = mdArgs;
}

Context::~Context() = default;

std::unique_ptr<Context> defaultContext()
{
    auto impl    = gmx::compat::make_unique<ContextImpl>();
    auto context = gmx::compat::make_unique<Context>(std::move(impl));
    return context;
}

} // end namespace gmxapi
