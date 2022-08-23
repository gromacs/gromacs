/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
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
/*! \file
 * \brief Bindings for System and session launch.
 *
 * \ingroup module_python
 * \author M. Eric Irrgang <ericirrgang@gmail.com>
 */

#include "gmxapi/session.h"
#include "gmxapi/status.h"
#include "gmxapi/system.h"
#include "gmxapi/version.h"

#include "module.h"
#include "pycontext.h"
#include "pysystem.h"

// Note that PyCapsule symbols from Python.h should be imported by way of the pybind headers, so
// let's not muddy the waters by explicitly including Python.h here unless we want to get more
// particular about the CMake configuration.

namespace gmxpy
{

namespace detail
{

namespace py = pybind11;


void export_system(py::module& m)
{
    using ::gmxapi::System;

    // The Session is a resource shared by active API operations.
    // We can't completely surrender ownership to Python because other API objects may refer to it.
    // We could use our own holder class instead of shared_ptr, but developers would
    // have to keep in mind that Python may make new references in different scopes
    // and threads, and pass references into other C++ code. Using shared_ptr
    // self-documents intent. Future implementations could refactor Session as a
    // dynamically accessed facet of the Context, which the API client would be
    // required to maintain and to pass to the API.
    py::class_<::gmxapi::Session, std::shared_ptr<::gmxapi::Session>> session(m, "MDSession");
    session.def("run", &::gmxapi::Session::run, "Run the simulation workflow");
    session.def("close",
                &::gmxapi::Session::close,
                "Shut down the execution environment and close the session.");

    // Export system container class
    py::class_<System, std::shared_ptr<System>> system(m, "MDSystem");
    system.def("launch", &launch, "Launch the configured workflow in the provided context.");

    // Module-level function
    m.def("from_tpr",
          &gmxpy::from_tpr,
          "Return a system container initialized from the given input record.");
}

} // namespace detail

} // end namespace gmxpy
