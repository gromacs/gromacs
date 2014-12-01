/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2008,2009,2010,2012,2014, by the GROMACS development team, led by
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
 * \brief This file declares functions for (mostly) the domdec module
 * to use MPI functionality
 *
 * \todo Wrap the raw dd_bcast in md.cpp into a higher-level function
 * in the domdec module, then this file can be module-internal.
 *
 * \author Berk Hess <hess@kth.se>
 * \inlibraryapi
 * \ingroup module_domdec
 */
#ifndef GMX_DOMDEC_DOMDEC_NETWORK_H
#define GMX_DOMDEC_DOMDEC_NETWORK_H

#include "gromacs/legacyheaders/typedefs.h"

/* \brief */
enum {
    dddirForward, dddirBackward
};

/*! \brief Move integers in the communication region one cell along
 * the domain decomposition
 *
 * Moves in the dimension indexed by ddimind, either forward
 * (direction=dddirFoward) or backward (direction=dddirBackward).
 */
void
dd_sendrecv_int(const gmx_domdec_t *dd,
                int ddimind, int direction,
                int *buf_s, int n_s,
                int *buf_r, int n_r);

/*! \brief Move reals in the comm. region one cell along the domain decomposition
 *
 * Moves in the dimension indexed by ddimind, either forward
 * (direction=dddirFoward) or backward (direction=dddirBackward).
 */
void
dd_sendrecv_real(const gmx_domdec_t *dd,
                 int ddimind, int direction,
                 real *buf_s, int n_s,
                 real *buf_r, int n_r);

/*! \brief Move revc's in the comm. region one cell along the domain decomposition
 *
 * Moves in dimension indexed by ddimind, either forward
 * (direction=dddirFoward) or backward (direction=dddirBackward).
 */
void
dd_sendrecv_rvec(const gmx_domdec_t *dd,
                 int ddimind, int direction,
                 rvec *buf_s, int n_s,
                 rvec *buf_r, int n_r);


/*! \brief Move revc's in the comm. region one cell along the domain decomposition
 *
 * Moves in dimension indexed by ddimind, simultaneously in the forward
 * and backward directions.
 */
void
dd_sendrecv2_rvec(const gmx_domdec_t *dd,
                  int ddimind,
                  rvec *buf_s_fw, int n_s_fw,
                  rvec *buf_r_fw, int n_r_fw,
                  rvec *buf_s_bw, int n_s_bw,
                  rvec *buf_r_bw, int n_r_bw);


/* The functions below perform the same operations as the MPI functions
 * with the same name appendices, but over the domain decomposition
 * nodes only.
 * The DD master node is the master for these operations.
 */

/*! \brief Broadcasts \p nbytes from \p data on \p DDMASTERRANK to all PP ranks */
void
dd_bcast(gmx_domdec_t *dd, int nbytes, void *data);

/*! \brief Copies \p nbytes from \p src to \p dest on \p DDMASTERRANK
 * and then broadcasts to \p dest on all PP ranks */
void
dd_bcastc(gmx_domdec_t *dd, int nbytes, void *src, void *dest);

/*! \brief Scatters \p nbytes from \p src on \p DDMASTERRANK to all PP ranks, received in \p dest */
void
dd_scatter(gmx_domdec_t *dd, int nbytes, void *src, void *dest);

/*! \brief Gathers \p nbytes from \p src on all PP ranks, received in \p dest on \p DDMASTERRANK */
void
dd_gather(gmx_domdec_t *dd, int nbytes, void *src, void *dest);

/*! \brief Scatters \p scounts bytes from \p src on \p DDMASTERRANK to all PP ranks, receiving \p rcount bytes in \p dest.
 *
 * See man MPI_Scatterv for details of how to construct scounts and disps.
 * If rcount==0, rbuf is allowed to be NULL */
void
dd_scatterv(gmx_domdec_t *dd,
            int *scounts, int *disps, void *sbuf,
            int rcount, void *rbuf);

/*! \brief Gathers \p rcount bytes from \p src on all PP ranks, received in \p scounts bytes in \p dest on \p DDMASTERRANK.
 *
 * See man MPI_Gatherv for details of how to construct scounts and disps.
 *
 * If scount==0, sbuf is allowed to be NULL */
void
dd_gatherv(gmx_domdec_t *dd,
           int scount, void *sbuf,
           int *rcounts, int *disps, void *rbuf);

#endif
