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
 * \brief Handlers for MPI details.
 *
 * \ingroup module_python
 * \author M. Eric Irrgang <ericirrgang@gmail.com>
 */


#ifndef GMXPY_MPI_BINDINGS_H
#define GMXPY_MPI_BINDINGS_H

#include <mpi.h>

#include "pybind11/pybind11.h"

#include "gmxapi/context.h"

namespace gmxpy
{

// Base exception forward declared for gmxpy_exceptions.h.
class Exception;

/*!
 * \brief Get a pointer to the MPI_Comm wrapped in an mpi4py Comm
 * \param communicator (:py:class:`mpi4py.MPI.Comm`): wrapped MPI communicator.
 * \return Pointer to C object.
 */
MPI_Comm* get_mpi_comm(pybind11::object communicator);

/*!
 * \brief Adapter to the gmxapi offer_comm protocol.
 * \param communicator (:py:class:`mpi4py.MPI.Comm`): wrapped MPI communicator.
 * \return gmxapi C++ Context handle.
 *
 * Implementation is selected by CMake, depending on available GROMACS library support.
 * See :file:`gmxapi/mpi/gmxapi_mpi.h` and the :cpp:func:`gmxapi::assignResource()` template.
 *
 * \throws if communicator is invalid (e.g. MPI_COMM_NULL) or not usable (e.g. too big)
 */
gmxapi::Context context_from_py_comm(pybind11::object communicator);

namespace detail
{

/*!
 * \brief Register bindings for MPI and MPI-enabled GROMACS, if possible.
 * \param m (:cpp:class:`pybind11::module`): The Python module that is in the process of being
 * imported. \param exception Module base exception from which additional exceptions should derive.
 */
void export_mpi_bindings(pybind11::module& m, const pybind11::exception<Exception>& exception);

} // end namespace detail

} // end namespace gmxpy

#endif // GMXPY_MPI_BINDINGS_H
