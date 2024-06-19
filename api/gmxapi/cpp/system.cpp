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

#include "gmxapi/system.h"

#include <array>
#include <memory>
#include <utility>
#include <vector>

#include "gromacs/mdrun/runner.h"
#include "gromacs/utility.h"
#include "gromacs/utility/gmxassert.h"

#include "gmxapi/context.h"
#include "gmxapi/md.h"
#include "gmxapi/session.h"
#include "gmxapi/status.h"

#include "system_impl.h"
#include "workflow.h"

namespace gmxapi
{

//! \cond
System::Impl::~Impl() = default;

System::Impl::Impl(System::Impl&&) noexcept = default;

System::Impl& System::Impl::operator=(System::Impl&& source) noexcept
{
    if (this != &source)
    {
        workflow_.swap(source.workflow_);
    }
    return *this;
}
//! \endcond

std::shared_ptr<Session> System::launch(const std::shared_ptr<Context>& context)
{
    return impl_->launch(context);
}

//! \cond
System::System(std::unique_ptr<Impl> implementation) : impl_{ std::move(implementation) }
{
    GMX_ASSERT(impl_, "Constructor requires valid implementation object.");
}
System::Impl* System::get() const
{
    return impl_.get();
}

System::~System() = default;
//! \endcond

System::System(System&&) noexcept = default;

System& System::operator=(System&&) noexcept = default;

System fromTprFile(const std::string& filename)
{
    // TODO Confirm the file is readable and parseable and note unique
    // identifying information for when the work spec is used in a different
    // environment.

    // Create a new Workflow instance.
    // TODO error handling
    auto workflow = Workflow::create(filename);

    // This may produce errors or throw exceptions in the future, but as of
    // 0.0.3 only memory allocation errors are possible, and we do not have a
    // plan for how to recover from them.
    auto systemImpl = std::make_unique<System::Impl>(std::move(workflow));
    GMX_ASSERT(systemImpl, "Could not create a valid implementation object.");
    auto system = System(std::move(systemImpl));

    // TODO Separate TPR information into appropriate abstractions.
    // The TPR file has enough information for us to
    //  1. choose an MD engine
    //  2. Get structure information
    //  3. Get topology information
    //  4. Get a lot of simulation and runtime parameters, but not all.
    // It does not have enough information on its own to determine much about the
    // necessary computation environment. That comes from environment
    // introspection and user runtime options.

    return system;
}

std::shared_ptr<Workflow> getWork(const System::Impl& system)
{
    return system.workflow_;
}

System::Impl::Impl(std::unique_ptr<gmxapi::Workflow> workflow) noexcept :
    workflow_(std::move(workflow)), spec_(std::make_shared<MDWorkSpec>())
{
    GMX_ASSERT(workflow_, "Class invariant implies non-null workflow_ member");
    GMX_ASSERT(spec_, "Class invariant implies non-null work specification member.");
}

std::shared_ptr<Session> System::Impl::launch(const std::shared_ptr<Context>& context)
{
    std::shared_ptr<Session> session = nullptr;
    if (context != nullptr)
    {
        // TODO: gmxapi::Workflow, gmxapi::MDWorkSpec, and gmxapi::MDModule need sensible consolidation.
        session = context->launch(*workflow_);
        GMX_ASSERT(session, "Context::launch() expected to produce non-null session.");

        for (auto&& module : spec_->getModules())
        {
            // TODO: This should be the job of the launching code that produces the Session.
            // Configure the restraints in a restraint manager made available to the session launcher.
            addSessionRestraint(session.get(), module);
        }
    }
    else
    {
        // we should log the error and return nullptr, but we have nowhere to set
        // a status object, by the described behavior.
    }

    return session;
}

} // end namespace gmxapi
