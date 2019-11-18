/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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
 * \brief Exception translation.
 *
 * Note that C++ exception objects are not actually bound to Python objects
 * because C++ cannot give ownership of a thrown exception to the interpreter.
 * Instead, we catch C++ exceptions and translate them into Python exceptions.
 *
 * As of pybind11 2.2, pybind does not provide a way to automatically express
 * C++ exception class inheritance to Python, and only supports simple Python
 * exceptions initialized with strings at the time of translation. (Automatic
 * translation uses `std::exception::what()`)
 *
 * Currently, we restrict ourselves to this simple translation. Future versions
 * may use the Python C API directly to support richer exception classes.
 *
 * \todo Determine inheritance relationship or equality between
 * gmxapi._gmxapi.Exception and gmxapi.exceptions.Error, if any.
 *
 * \author M. Eric Irrgang <ericirrgang@gmail.com>
 *
 * \ingroup module_python
 */
#include "gmxapi/exceptions.h"

#include "module.h"


namespace gmxpy
{

namespace detail
{

namespace py = pybind11;


void export_exceptions(pybind11::module& m)
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

    static py::exception<gmxapi::Exception> baseException(m, "Exception", gmxapi_exceptions_error.ptr());

    // Developer note: as of pybind11 2.2.4, the py::exception template argument
    // is unused internally, but required.
    struct UnknownExceptionPlaceHolder
    {
    };
    static py::exception<UnknownExceptionPlaceHolder> unknownException(m, "UnknownException",
                                                                       baseException.ptr());
    unknownException.doc() =
            "GROMACS library produced an exception that is "
            "not mapped in gmxapi or which should have been "
            "caught at a lower level. I.e. a bug. (Please report.)";

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

    // Map gmxapi exceptions from gmxapi/exceptions.h to Python package exceptions.
    // Note: C++ exception translation occurs in revers order of registration,
    // So derived exceptions must be registered after their base exceptions.
    // TODO: We could have more informative exception class docstrings
    //   by linking to online docs or if we had a way to reuse doxygen docs.

    {
        auto exception =
                py::register_exception<gmxapi::ProtocolError>(m, "ProtocolError", baseException.ptr());
        exception.doc() = "Behavioral protocol violated.";
    }

    {
        auto exception = py::register_exception<gmxapi::NotImplementedError>(
                m, "NotImplementedError", baseException.ptr());
        exception.doc() = "Expected feature is not implemented.";
    }

    {
        auto exception =
                py::register_exception<gmxapi::UsageError>(m, "UsageError", baseException.ptr());
        exception.doc() = "Unacceptable API usage.";
    }
}


} // namespace detail

} // end namespace gmxpy
