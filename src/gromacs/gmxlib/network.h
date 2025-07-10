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
#ifndef GMX_GMXLIB_NETWORK_H
#define GMX_GMXLIB_NETWORK_H

/*
 * This module defines the interface of the actual communication routines.
 */

#include <cstdio>

#include <memory>

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/stringutil.h"

struct t_commrec;
struct t_filenm;

//! Free memory associated with the commrec.
void done_commrec(t_commrec* cr);

//! Wait until all processes in communicator have reached the barrier
void gmx_barrier(MPI_Comm communicator);

/*! \brief Broadcast nbytes bytes from the main to communicator
 *
 * Can be called with a single rank or without MPI
 */
void gmx_bcast(std::size_t nbytes, void* b, MPI_Comm communicator);


const char* opt2fn_main(const char* opt, int nfile, const t_filenm fnm[], t_commrec* cr);
/* Return the filename belonging to cmd-line option opt, or NULL when
 * no such option or not running on main */

[[noreturn]] void gmx_fatal_collective(int                    f_errno,
                                       const char*            file,
                                       int                    line,
                                       MPI_Comm               comm,
                                       gmx_bool               bMain,
                                       gmx_fmtstr const char* fmt,
                                       ...) gmx_format(printf, 6, 7);
/* As gmx_fatal declared in utility/fatalerror.h,
 * but only the main process prints the error message.
 * This should only be called one of the following two situations:
 * 1) On all nodes in cr->mpi_comm_mysim, with cr!=NULL,dd==NULL.
 * 2) On all nodes in dd->mpi_comm_all,   with cr==NULL,dd!=NULL.
 * This will call MPI_Finalize instead of MPI_Abort when possible,
 * This is useful for handling errors in code that is executed identically
 * for all processes.
 */

#endif
