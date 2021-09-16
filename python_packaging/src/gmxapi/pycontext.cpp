/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019,2020,2021, by the GROMACS development team, led by
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
 * \brief Wrapper code for gmxapi::Context.
 *
 * \ingroup module_python
 * \author M. Eric Irrgang <ericirrgang@gmail.com>
 */
#include "pycontext.h"

#include "gmxapi/exceptions.h"
#include "gmxapi/gmxapi.h"
#include "gmxapi/md.h"
#include "gmxapi/session.h"
#include "gmxapi/status.h"
#include "gmxapi/version.h"

namespace py = pybind11;

namespace gmxpy
{

void PyContext::setMDArgs(const MDArgs& mdArgs)
{
    assert(context_);
    context_->setMDArgs(mdArgs);
}

std::shared_ptr<gmxapi::Session> PyContext::launch(const gmxapi::Workflow& work)
{
    assert(context_);
    std::shared_ptr<gmxapi::Session> session = nullptr;

    // TODO: gmxapi::Workflow, gmxapi::MDWorkSpec, and gmxapi::MDModule need sensible consolidation.
    session = gmxapi::launchSession(context_.get(), work);
    if (!session)
    {
        throw gmxapi::ProtocolError("Context::launch() expected to produce non-null session.");
    }

    for (auto&& module : workNodes_->getModules())
    {
        // TODO: This should be the job of the launching code that produces the Session.
        // Configure the restraints in a restraint manager made available to the session launcher.
        auto status = gmxapi::addSessionRestraint(session.get(), module);
    }

    return session;
}

std::shared_ptr<gmxapi::MDWorkSpec> PyContext::getSpec() const
{
    assert(workNodes_);
    return workNodes_;
}

std::shared_ptr<gmxapi::Context> PyContext::get() const
{
    assert(context_);
    return context_;
}

PyContext::PyContext() :
    context_{ std::make_shared<gmxapi::Context>(gmxapi::createContext()) },
    workNodes_{ std::make_shared<gmxapi::MDWorkSpec>() }
{
    assert(context_);
    assert(workNodes_);
}

void PyContext::addMDModule(const pybind11::object& force_object) const
{
    if (!::gmxapi::Version::isAtLeast(0, 2, 1))
    {
        throw ::gmxapi::NotImplementedError("Feature requires gmxapi 0.2.1 with GROMACS 2021.3.");
    }
    // If force_object has a bind method, give it a PyCapsule with a pointer
    // to our C++ object.
    if (py::hasattr(force_object, "bind"))
    {
        auto spec     = getSpec();
        auto holder   = new gmxapi::MDHolder(spec);
        holder->name_ = "pygmx holder";
        auto deleter  = [](PyObject* o) {
            if (PyCapsule_IsValid(o, gmxapi::MDHolder_Name))
            {
                auto holder_ptr = (gmxapi::MDHolder*)PyCapsule_GetPointer(o, gmxapi::MDHolder_Name);
                delete holder_ptr;
                // \todo double-check whether there is something we should do to invalidate a PyCapsule.
            }
        };
        auto       capsule = py::capsule(holder, gmxapi::MDHolder_Name, deleter);
        py::object bind    = force_object.attr("bind");
        // py::capsule does not have bindings and does not implicitly convert to py::object
        py::object obj = capsule;
        bind(obj);
    }
    else
    {
        // Note: Exception behavior is likely to change.
        // Ref: https://github.com/kassonlab/gmxapi/issues/125
        throw py::value_error("Argument must provide a `bind` method.");
    }
} // namespace gmxpy

} // end namespace gmxpy
