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
 * \brief Bindings for Context class.
 *
 *
 * \author M. Eric Irrgang <ericirrgang@gmail.com>
 *
 * \ingroup module_python
 */
#include "gmxapi/context.h"
#include "gmxapi/exceptions.h"
#include "gmxapi/version.h"

#include "gmxpy_exceptions.h"
#include "module.h"
#include "pycontext.h"


namespace gmxpy
{

namespace detail
{

namespace py = pybind11;

/*!
 * \brief Normalize argument value types and construct argv.
 *
 * \param params named parameters from gmxapi 0.0.7
 * \return Vector for mdrun argv.
 */
static std::vector<std::string> makeMDArgs_0_0_7(const py::dict& params)
{
    std::vector<std::string> mdargs;
    if (params.contains("grid"))
    {
        if (py::len(params["grid"]) == 0)
        {
            throw gmxapi::UsageError("Grid argument must describe domain decomposition grid.");
        }
        std::vector<std::string> vals;
        auto                     iterator = py::iter(params["grid"]);
        while (iterator != py::iterator::sentinel())
        {
            vals.emplace_back(py::cast<std::string>(py::str(iterator)));
            ++iterator;
        }
        mdargs.emplace_back("-dd");
        for (auto&& val : vals)
        {
            mdargs.emplace_back(val);
        }
    }
    if (params.contains("pme_ranks"))
    {
        auto val = py::cast<std::string>(py::str(params["pme_ranks"]));
        mdargs.emplace_back("-npme");
        mdargs.emplace_back(val);
    }
    if (params.contains("threads"))
    {
        auto val = py::cast<std::string>(py::str(params["threads"]));
        mdargs.emplace_back("-nt");
        mdargs.emplace_back(val);
    }
    if (params.contains("tmpi"))
    {
        auto val = py::cast<std::string>(py::str(params["tmpi"]));
        mdargs.emplace_back("-ntmpi");
        mdargs.emplace_back(val);
    }
    if (params.contains("threads_per_rank"))
    {
        auto val = py::cast<std::string>(py::str(params["threads_per_rank"]));
        mdargs.emplace_back("-ntomp");
        mdargs.emplace_back(val);
    }
    if (params.contains("pme_threads_per_rank"))
    {
        auto val = py::cast<std::string>(py::str(params["threads_per_pme_rank"]));
        mdargs.emplace_back("-ntomp_pme");
        mdargs.emplace_back(val);
    }
    if (params.contains("steps"))
    {
        auto val = py::cast<std::string>(py::str(params["steps"]));
        mdargs.emplace_back("-nsteps");
        mdargs.emplace_back(val);
    }
    if (params.contains("max_hours"))
    {
        auto val = py::cast<std::string>(py::str(params["max_hours"]));
        mdargs.emplace_back("-maxh");
        mdargs.emplace_back(val);
    }
    if (params.contains("append_output"))
    {
        if (!params["append_output"].cast<bool>())
        {
            mdargs.emplace_back("-noappend");
        }
    }
    return mdargs;
}

static std::vector<std::string> makeMDArgs_CLI(const py::dict& params)
{
    // Make sure pybind `None` has the same auto-conversion to `NoneType` as in Python.
    assert(py::isinstance<py::none>(py::none()));

    std::vector<std::string> mdargs;
    // for key, value in params, if key.startswith('-'): mdargs.append(key); mdargs.extend(*value)
    for (auto items : params)
    {
        auto key   = items.first;
        auto value = items.second;
        auto arg   = py::cast<std::string>(key);
        if (arg.front() == '-')
        {
            mdargs.emplace_back(arg);
            if (py::isinstance<py::none>(value))
            {
                continue;
            }
            if (py::isinstance<py::str>(value))
            {
                mdargs.emplace_back(py::cast<std::string>(value));
                continue;
            }
            // If value is a non-string iterator, get each value.
            try
            {
                auto it = py::iter(value);
                while (it != py::iterator::sentinel())
                {
                    mdargs.emplace_back(py::str(*it));
                    ++it;
                }
                continue;
            }
            catch (py::error_already_set& eas)
            {
                // `py::iter` throws TypeError for bad arguments, but it comes through the C API,
                // so pybind11 wraps it in py::error_already_set, not py::type_error. Ref:
                // https://pybind11.readthedocs.io/en/stable/advanced/exceptions.html#handling-exceptions-from-python-in-c
                // Only suppress TypeError.
                if (!eas.matches(PyExc_TypeError))
                {
                    throw;
                }
            }
            // Otherwise, cast the Python value to a Python string (which is implicitly convertible
            // to std::string).
            mdargs.emplace_back(py::str(value));
        }
    }
    return mdargs;
}


/*! \internal
 * \brief Update a parameter structure for a simulation execution context.
 *
 * \param mdargs [OUT] Container for parameters made available to the MD library context.
 * \param params Python dictionary mapping parameter names to argument values.
 *
 * Note that the current library infrastructure does not provide a way for the
 * simulation machinery to express human-readable parameter names with rich
 * descriptions, so a few of the most necessary mdrun command line parameters
 * are hard coded here. Ref. https://gitlab.com/gromacs/gromacs/-/issues/2877
 *
 * For reference and default values, see
 * http://manual.gromacs.org/current/onlinehelp/gmx-mdrun.html#options
 */
static void setMDArgs(std::vector<std::string>* mdargs, const py::dict& params)
{
    // Note: params is processed twice, but entries may be silently ignored if
    // neither consumer is interested.
    // TODO: Consider copying `params` and popping values as they are processed.

    // Get key-word mapped arguments from gmxapi 0.0.7
    auto args_0_0_7 = makeMDArgs_0_0_7(params);

    // Add raw hyphen-prefixed CLI args and string values without pre-checking.
    // This doesn't check for overlap between gmxapi 0.0.7 input and gmxapi 0.1+ input.
    // For the moment, the user is at the mercy of CLI input-checking behavior.
    auto args_cli = makeMDArgs_CLI(params);

    // This function takes complete control of the values in mdargs, and so we
    // clear any leftover values from a previously-used Context.
    mdargs->reserve(args_0_0_7.size() + args_cli.size());
    *mdargs = std::move(args_0_0_7);
    if (gmxapi::Version::isAtLeast(0, 3))
    {
        mdargs->insert(std::end(*mdargs), std::begin(args_cli), std::end(args_cli));
    }
    else
    {
        // Before 0.3.0, mdrun parameters were strictly curated by gmxapi. The 0.3 Python package
        // introduced the `runtime_args` key word parameter to `gmxapi.simulation.mdrun()`, but if
        // the gmxapi library is older than 0.3.0, only the gmxapi 0.0.7 mdrun parameters are
        // supported.
        if (args_cli.size() > 0)
        {
            const auto message =
                    std::string("Invalid runtime_args for libgmxapi ") + gmxapi::Version::release();
            throw gmxapi::UsageError(message);
        }
    }
}

void export_context(pybind11::module& m, const pybind11::exception<Exception>& baseException)
{
    // Add argument type before it is used for more sensible automatic bindings behavior.
    py::class_<MDArgs, std::unique_ptr<MDArgs>> mdargs(m, "MDArgs");
    mdargs.def(py::init(), "Create an empty MDArgs object.");
    mdargs.def(
            "set",
            [](MDArgs* self, const py::dict& params) { setMDArgs(self, params); },
            "Assign parameters in MDArgs from Python dict.");
    mdargs.def(
            "get_args",
            [](const MDArgs& self) {
                auto args = py::list(self.size());
                for (int i = 0; i < self.size(); ++i)
                {
                    args[i] = py::str(self[i]);
                };
                return args;
            },
            "Get an iterator of command line argument tokens, if possible and relevant.");

    // Export execution context class
    py::class_<PyContext, std::shared_ptr<PyContext>> context(m, "Context");
    context.def(py::init(), "Create a default execution context.");
    context.def("setMDArgs", &PyContext::setMDArgs, "Set MD runtime parameters.");

    context.def("add_mdmodule", &PyContext::addMDModule, "Add an MD plugin for the simulation.");

    export_create_context(m, baseException);
}

} // namespace detail

} // end namespace gmxpy
