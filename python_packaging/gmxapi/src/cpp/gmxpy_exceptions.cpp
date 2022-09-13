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
 * \brief Basic translation from C++ to Python Exceptions.
 *
 * Provide definitions for the common ::gmxpy C++ exception classes
 * and export corresponding Python exceptions.
 *
 * Known exceptions from supporting libraries are deferred to an
 * export function specific to the API level that is detected by CMake at
 * build time. (See wrapped_exceptions.h)
 *
 * \author M. Eric Irrgang <ericirrgang@gmail.com>
 *
 * \ingroup module_python
 */

#include "gmxpy_exceptions.h"

#include "gmxapi/exceptions.h"

#include "wrapped_exceptions.h"

namespace gmxpy
{

FeatureNotAvailable::FeatureNotAvailable(const std::string& what_arg) : Exception(what_arg) {}

FeatureNotAvailable::FeatureNotAvailable() : FeatureNotAvailable(std::string()) {}

FeatureNotAvailable::FeatureNotAvailable(const char* what_arg) :
    FeatureNotAvailable(std::string(what_arg))
{
}
FeatureNotAvailable& FeatureNotAvailable::operator=(const FeatureNotAvailable&) = default;

FeatureNotAvailable::FeatureNotAvailable(const FeatureNotAvailable&) = default;

FeatureNotAvailable& FeatureNotAvailable::operator=(FeatureNotAvailable&&) noexcept = default;

FeatureNotAvailable::FeatureNotAvailable(FeatureNotAvailable&&) noexcept = default;

FeatureNotAvailable::~FeatureNotAvailable() = default;


namespace detail
{

namespace py = pybind11;


const pybind11::exception<Exception>& export_exceptions(pybind11::module& m)
{

    // These two lines could cause exceptions, but they are already handled,
    // causing an ImportError for the _gmxapi submodule raised from the
    // ImportError or AttributeError that caused the failure.
    // If we find that this is too cryptic or that there are too many places
    // that import could fail, we can use the Python C API directly
    // (in the module.cpp code block preceded by the PYBIND11 macro) to set a
    // new PyExc_ImportError manually and return a nullptr immediately.
    const auto packageModule           = py::module::import("gmxapi");
    const auto gmxapi_exceptions_error = packageModule.attr("exceptions").attr("Error");

    static py::exception<Exception> baseException(m, "Exception", gmxapi_exceptions_error.ptr());

    // Developer note: as of pybind11 2.2.4, the py::exception template argument
    // is unused internally, but required.
    struct UnknownExceptionPlaceHolder
    {
    };
    static py::exception<UnknownExceptionPlaceHolder> unknownException(
            m, "UnknownException", baseException.ptr());
    unknownException.doc() =
            R"delimeter(Catch-all exception wrapper.

            GROMACS library produced an exception that is not mapped in gmxapi
            or which should have been caught at a lower level. I.e. a bug.
            (Please report.)
            )delimeter";

    // Catch unexpected/unbound exceptions from libgromacs or libgmxapi.
    py::register_exception_translator([](std::exception_ptr p) {
        try
        {
            if (p)
            {
                std::rethrow_exception(p);
            }
        }
        catch (const gmxapi::Exception& e)
        {
            // Nothing should be throwing the base exception and all gmxapi
            // exceptions should be mapped in this module. Differences
            // between GROMACS version and Python package version could leave
            // some exceptions unmapped, but we should add an alert in case
            // an exception gets overlooked.
            std::string message = "Generic gmxapi exception caught: ";
            message += e.what();
            baseException(message.c_str());
        }
        catch (const std::exception& e)
        {
            std::string message = "Please report GROMACS bug. Unhandled C++ exception: ";
            message += e.what();
            unknownException(message.c_str());
        }
    });

    {
        auto exception = py::register_exception<FeatureNotAvailable>(
                m, "FeatureNotAvailable", baseException.ptr());
        exception.doc() =
                R"delimeter(An API feature is not available in the current installation.

                This may occur when a new gmxapi Python package is installed with an older
                GROMACS installation that does not have the library support for a newer
                feature.
                )delimeter";
    }

    export_wrapped_exceptions(m, baseException);

    return baseException;
}


} // namespace detail

} // end namespace gmxpy
