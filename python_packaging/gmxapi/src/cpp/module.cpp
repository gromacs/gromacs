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

/*! \internal \file
 * \brief Exports Python bindings for gmxapi._gmxapi module.
 *
 * Defines the entry point for an importable Python extension module
 * (in accordance with Python C API), using the pybind11 template headers.
 *
 * \author M. Eric Irrgang <ericirrgang@gmail.com>
 *
 * \ingroup module_python
 */

#include "module.h"

#include <memory>

#include "pybind11/pybind11.h"

#include "gmxapi/status.h"
#include "gmxapi/version.h"

#include "gmxpy_exceptions.h"

namespace py = pybind11;

// Export Python module.

/// used to set __doc__
/// pybind11 uses const char* objects for docstrings. C++ raw literals can be used.
const char* const docstring = R"delimeter(
gmxapi core module
==================

gmxapi._gmxapi provides Python access to the GROMACS C++ API so that client code can be
implemented in Python, C++, or a mixture. The classes provided are mirrored on the
C++ side in the gmxapi namespace as best as possible.

This documentation is generated from C++ extension code. Refer to C++ source
code and developer documentation for more details.

)delimeter";

/*! \brief Export gmxapi._gmxapi Python module in shared object file.
 *
 * \ingroup module_python
 */

// Instantiate the Python module
PYBIND11_MODULE(_gmxapi, m)
{
    using namespace gmxpy::detail;
    m.doc() = docstring;

    // Provide a module level dict for internal use to support
    // API queries for feature level.
    m.attr("_named_features") = py::dict();

    // Register exceptions and catch-all exception translators. We do this early
    // to give more freedom to the other export functions. Note that bindings
    // for C++ symbols should be expressed before those symbols are referenced
    // in other bindings, and that exception translators are tried in reverse
    // order of registration for uncaught C++ exceptions.
    const auto& baseException = export_exceptions(m);

    // Export core bindings
    m.def("library_has_feature",
          &gmxapi::Version::hasFeature,
          "Check the gmxapi library for a named feature.");

    py::class_<::gmxapi::Status> gmx_status(m, "Status", "Holds status for API operations.");

    // Get bindings exported by the various components.
    // Additional exports may be conditionally added within this export functions
    // based on implementations chosen by CMake logic. (E.g. MPI bindings or features
    // requiring specific GROMACS versions)
    export_context(m, baseException);
    export_system(m);
    export_tprfile(m);

    // Module helpers and utilities
    m.def(
            "has_feature",
            [self = m](const std::string& name) {
                py::gil_scoped_acquire lock;
                bool feature_found = py::cast<py::dict>(self.attr("_named_features")).contains(name);
                if (!feature_found)
                {
                    feature_found = py::cast<bool>(self.attr("library_has_feature")(name));
                }
                return feature_found;
            },
            "Check feature *name* first with the bindings package, then the supporting library.");

} // end pybind11 module
