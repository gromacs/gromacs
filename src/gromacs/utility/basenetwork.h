/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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
/*! \libinternal \file
 * \brief
 * Utility functions for basic MPI and network functionality.
 *
 * \inlibraryapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_BASENETWORK_H
#define GMX_UTILITY_BASENETWORK_H

#include <stddef.h>

#include "gromacs/utility/basedefinitions.h"

#ifdef __cplusplus
extern "C" {
#endif

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
gmx_bool gmx_mpi_initialized(void);

/*! \brief
 * Returns the number of nodes.
 *
 * For thread-MPI, returns one before the threads have been started.
 * This allows code between the real MPI_Init() and the thread-MPI "init" to
 * still use this function to check for serial/parallel status and work as
 * expected: for thread-MPI, at that point they should behave as if the run was
 * serial.
 */
int gmx_node_num(void);

/*! \brief
 * Returns the rank of the node.
 *
 * For thread-MPI, returns zero before the threads have been started.
 * This allows code between the real MPI_Init() and the thread-MPI "init" to
 * still use this function to check for master node work as expected:
 * for thread-MPI, at that point the only thread of execution should behave as
 * if it the master node.
 */
int gmx_node_rank(void);

/*! \brief
 * Return a non-negative hash that is, hopefully, unique for each physical
 * node.
 *
 * This hash is useful for determining hardware locality.
 */
int gmx_physicalnode_id_hash(void);

/*! \brief
 * Returns an integer characteristic of and, hopefully, unique to each
 * physical node in the MPI system.
 *
 * If the first part of the MPI hostname (up to the first dot) ends with a
 * number, returns this number.  If the first part of the MPI hostname does not
 * ends in a number (0-9 characters), returns 0.
 *
 * \todo
 * This function should be fully replaced by gmx_physicalnode_id_hash().
 */
int gmx_hostname_num(void);

/** Abort the parallel run */
void gmx_abort(int errorno);

#ifdef __cplusplus
}
#endif

#endif
