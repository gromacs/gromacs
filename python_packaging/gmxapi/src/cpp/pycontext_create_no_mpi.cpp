/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2022- The GROMACS Authors
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
 * \brief Create a PyContext for gmxapi < 0.2.0
 *
 *
 * \author M. Eric Irrgang <ericirrgang@gmail.com>
 *
 * \ingroup module_python
 */


#include "pybind11/pybind11.h"

#include "gmxpy_exceptions.h"
#include "pycontext.h"

namespace py = pybind11;

namespace gmxpy
{

PyContext create_context()
{
    auto context     = gmxapi::createContext();
    auto context_ptr = std::make_shared<gmxapi::Context>(std::move(context));
    return PyContext(std::move(context_ptr));
}

PyContext create_context(py::object /* unused */)
{
    throw FeatureNotAvailable(
            "ResourceAssignment bindings require gmxapi library level 0.2 or higher and a "
            "compatible mpi4py installation.");
}
namespace detail
{

void export_create_context(pybind11::module& m, const pybind11::exception<Exception>& /*exception*/)
{
    py::dict features          = m.attr("_named_features");
    features["create_context"] = 1;
    m.def(
            "create_context",
            []() { return create_context(); },
            "Initialize a new API Context to manage resources and software environment.");
    m.def(
            "create_context",
            [](const py::object& resource) { return create_context(resource); },
            "Initialize a new API Context to manage resources and software environment.");
}
} // namespace detail

} // namespace gmxpy
