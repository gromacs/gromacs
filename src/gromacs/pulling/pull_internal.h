/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
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
 *
 *
 * \brief
 * This file contains datatypes and function declarations for internal
   use in the pull code.
 *
 * \author Berk Hess
 *
 * \inlibraryapi
 */

#ifndef GMX_PULLING_PULL_INTERNAL_H
#define GMX_PULLING_PULL_INTERNAL_H

#include "gromacs/legacyheaders/typedefs.h"

#ifdef __cplusplus
extern "C" {
#endif


struct pull_comm_t {
    gmx_bool  bParticipateAll; /* Do all ranks always participate in pulling? */
    gmx_bool  bParticipate;    /* Does our rank participate in pulling? */
#ifdef GMX_MPI
    MPI_Comm  mpi_comm_com;    /* Communicator for pulling */
#endif
    int       nparticipate;    /* The number of ranks participating */

    gmx_int64_t setup_count;   /* The number of decomposition calls */
    gmx_int64_t must_count;    /* The last count our rank needed to be part */

    gmx_bool  bSetPBCatoms;    /* Should we set x_pbc for the groups? */

    rvec     *rbuf;            /* COM calculation buffer */
    dvec     *dbuf;            /* COM calculation buffer */
    double   *dbuf_cyl;        /* cylinder ref. groups COM calculation buffer */

    FILE     *out_x;           /* output file for pull data */
    FILE     *out_f;           /* output file for pull data */
};


#ifdef __cplusplus
}
#endif

#endif
