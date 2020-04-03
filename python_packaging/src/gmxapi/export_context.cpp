/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019,2020, by the GROMACS development team, led by
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
 * \brief Bindings for Context class.
 *
 *
 * \author M. Eric Irrgang <ericirrgang@gmail.com>
 *
 * \ingroup module_python
 */
#include "module.h"

#include "gmxapi/context.h"
#include "gmxapi/exceptions.h"
#include "pycontext.h"


namespace gmxpy
{

namespace detail
{

namespace py = pybind11;


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
    mdargs->clear();
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
        mdargs->emplace_back("-dd");
        for (auto&& val : vals)
        {
            mdargs->emplace_back(val);
        }
    }
    if (params.contains("pme_ranks"))
    {
        auto val = py::cast<std::string>(py::str(params["pme_ranks"]));
        mdargs->emplace_back("-npme");
        mdargs->emplace_back(val);
    }
    if (params.contains("threads"))
    {
        auto val = py::cast<std::string>(py::str(params["threads"]));
        mdargs->emplace_back("-nt");
        mdargs->emplace_back(val);
    }
    if (params.contains("tmpi"))
    {
        auto val = py::cast<std::string>(py::str(params["tmpi"]));
        mdargs->emplace_back("-ntmpi");
        mdargs->emplace_back(val);
    }
    if (params.contains("threads_per_rank"))
    {
        auto val = py::cast<std::string>(py::str(params["threads_per_rank"]));
        mdargs->emplace_back("-ntomp");
        mdargs->emplace_back(val);
    }
    if (params.contains("pme_threads_per_rank"))
    {
        auto val = py::cast<std::string>(py::str(params["threads_per_pme_rank"]));
        mdargs->emplace_back("-ntomp_pme");
        mdargs->emplace_back(val);
    }
    if (params.contains("steps"))
    {
        auto val = py::cast<std::string>(py::str(params["steps"]));
        mdargs->emplace_back("-nsteps");
        mdargs->emplace_back(val);
    }
    if (params.contains("max_hours"))
    {
        auto val = py::cast<std::string>(py::str(params["max_hours"]));
        mdargs->emplace_back("-maxh");
        mdargs->emplace_back(val);
    }
    if (params.contains("append_output"))
    {
        try
        {
            if (!params["append_output"].cast<bool>())
            {
                mdargs->emplace_back("-noappend");
            }
        }
        catch (const py::cast_error& e)
        {
            // Couldn't cast to bool for some reason.
            // Convert to gmxapi exception (not implemented)
            // ref. https://github.com/kassonlab/gmxapi/issues/125
            throw;
        }
    }
}

void export_context(py::module& m)
{
    // Add argument type before it is used for more sensible automatic bindings behavior.
    py::class_<MDArgs, std::unique_ptr<MDArgs>> mdargs(m, "MDArgs");
    mdargs.def(py::init(), "Create an empty MDArgs object.");
    mdargs.def("set", [](MDArgs* self, const py::dict& params) { setMDArgs(self, params); },
               "Assign parameters in MDArgs from Python dict.");

    // Export execution context class
    py::class_<PyContext, std::shared_ptr<PyContext>> context(m, "Context");
    context.def(py::init(), "Create a default execution context.");
    context.def("setMDArgs", &PyContext::setMDArgs, "Set MD runtime parameters.");

    context.def("add_mdmodule", &PyContext::addMDModule, "Add an MD plugin for the simulation.");
}

} // namespace detail

} // end namespace gmxpy
