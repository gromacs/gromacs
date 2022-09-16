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

#include "mpi_bindings.h"

#include "mpi4py/mpi4py.h"
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "gmxpy_exceptions.h"


namespace py = pybind11;

namespace gmxpy
{

MPI_Comm* get_mpi_comm(pybind11::object communicator)
{
    py::gil_scoped_acquire lock;

    MPI_Comm* comm_ptr = PyMPIComm_Get(communicator.ptr());

    if (comm_ptr == nullptr)
    {
        throw py::error_already_set();
    }
    return comm_ptr;
}

namespace detail
{

std::array<int, 2> mpi_report(MPI_Comm comm)
{
    int size = 0;
    MPI_Comm_size(comm, &size);

    int rank = 0;
    MPI_Comm_rank(comm, &rank);

    return { rank, size };
}

void export_mpi_bindings(pybind11::module& m, const pybind11::exception<Exception>&
                         /*exception*/)
{
    py::dict features        = m.attr("_named_features");
    features["mpi_bindings"] = 1;

    if (import_mpi4py() < 0)
    {
        throw py::error_already_set();
    }
    m.def(
            "mpi_report",
            [](py::object py_comm) {
                MPI_Comm* comm = get_mpi_comm(py_comm);
                return mpi_report(*comm);
            },
            R"pbdoc(
               Parameters of the MPI context: (rank, size)
        )pbdoc");
}

} // namespace detail
} // end namespace gmxpy
