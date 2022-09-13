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
/*! \defgroup module_python Python module for accessing Gromacs library
 * The Python module ``gmxapi`` consists of a high-level interface implemented in
 * pure Python and a low-level interface implemented as a C++ extension in the
 * submodule, gmxapi._gmxapi.
 */
/*! \file
 * \brief Declares symbols to be exported to gmxapi._gmxapi Python module.
 *
 * Declares namespace gmxpy, used internally in the C++ extension.
 * \author M. Eric Irrgang <ericirrgang@gmail.com>
 * \ingroup module_python
 */
#ifndef GMXPY_MODULE_H
#define GMXPY_MODULE_H

#include "pybind11/pybind11.h"

// Declaration for the module initialization.
// Avoids a warning for the use of PYBIND11_MODULE macro in module.cpp
// e.g. no previous prototype for function 'PyInit__gmxapi' [-Wmissing-prototypes]
extern "C" __attribute__((visibility("default"))) PyObject* PyInit__gmxapi();

/*! \brief API client code from which to export Python bindings
 *
 * gmxpy is not a public interface. It implements bindings for the public
 * Python API in the C++ Python extension it produces, and it uses the public
 * C++ Gromacs API, but is itself an API *client* and its C++ interfaces are not
 * intended to be used in external code.
 *
 * Python objects exported by the module are in the gmxpy namespace if
 * implemented in this library, or in the gmxapi namespace if directly wrapped
 * from libgmxapi in the GROMACS installation.
 *
 * \ingroup module_python
 */
namespace gmxpy
{

/*!
 * \brief Implementation details of the Python bindings library.
 *
 * Python objects (including classes/types) are expressed with pybind11.
 * At run time (when the module is imported by the Python interpreter), an
 * initialization function (defined with a pybind11 macro in module.cpp) is run
 * by the Python C API machinery to bind the code exported by this module.
 *
 * To keep module.cpp concise, we create the Python module object and then pass
 * it to various `export_...` functions defined in other source files. Those
 * export functions are declared in gmxpy::detail in module.h
 *
 * \ingroup module_python
 */
namespace detail
{

void export_context(pybind11::module& m);
void export_system(pybind11::module& m);
void export_tprfile(pybind11::module& module);

} // namespace detail

} // end namespace gmxpy

#endif // GMXPY_MODULE_H
