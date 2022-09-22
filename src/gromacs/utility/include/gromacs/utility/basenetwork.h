/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
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
/*! \libinternal \file
 * \brief
 * Utility functions for basic MPI and network functionality.
 *
 * \inlibraryapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_BASENETWORK_H
#define GMX_UTILITY_BASENETWORK_H

/*! \brief
 * Returns whether MPI has been initialized.
 *
 * The return value is `FALSE` if MPI_Init() has not been called, or if
 * \Gromacs has been compiled without MPI support.
 * For thread-MPI, returns `TRUE` when the threads have been started.
 *
 * Note that there is a lot of code in between MPI_Init() and the thread-MPI
 * thread start where the return value is different depending on compilation
 * options.
 */
bool gmx_mpi_initialized();

/*! \brief
 * Returns the number of nodes.
 *
 * For thread-MPI, returns one before the threads have been started.
 * This allows code between the real MPI_Init() and the thread-MPI "init" to
 * still use this function to check for serial/parallel status and work as
 * expected: for thread-MPI, at that point they should behave as if the run was
 * serial.
 */
int gmx_node_num();

/*! \brief
 * Returns the rank of the node.
 *
 * For thread-MPI, returns zero before the threads have been started.
 * This allows code between the real MPI_Init() and the thread-MPI "init" to
 * still use this function to check for main node work as expected:
 * for thread-MPI, at that point the only thread of execution should behave as
 * if it the main node.
 */
int gmx_node_rank();

/*! \brief
 * Return a non-negative hash that is, hopefully, unique for each physical
 * node.
 *
 * This hash is useful for determining hardware locality.
 */
int gmx_physicalnode_id_hash();

/** Abort the parallel run */
[[noreturn]] void gmx_abort(int errorno);

#endif
