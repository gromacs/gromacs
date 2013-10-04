/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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

#ifndef _network_h
#define _network_h


/*
 * This module defines the interface of the actual communication routines.
 */

#include <stdio.h>
#include "visibility.h"
#include "types/simple.h"
#include "types/commrec.h"
#include "typedefs.h"
#include "main.h"
#include "gmx_fatal.h"

#ifdef __cplusplus
extern "C" {
#endif

int gmx_setup(int *argc, char **argv, int *nnodes);
/* Initializes the parallel communication, return the ID of the node */

int gmx_node_num(void);
/* return the number of nodes in the ring */

int gmx_node_rank(void);
/* return the rank of the node */

GMX_LIBGMX_EXPORT
int gmx_physicalnode_id_hash(void);
/* Return a non-negative hash that is, hopefully, unique for each physical node.
 * This hash is useful for determining hardware locality.
 */

GMX_LIBGMX_EXPORT
int gmx_hostname_num(void);
/* Ostensibly, returns a integer characteristic of and unique to each
   physical node in the MPI system. If the first part of the MPI
   hostname (up to the first dot) ends with a number, returns this
   number. If the first part of the MPI hostname does not ends in a
   number (0-9 characters), returns 0.
 */

GMX_LIBGMX_EXPORT
void gmx_setup_nodecomm(FILE *fplog, t_commrec *cr);
/* Sets up fast global communication for clusters with multi-core nodes */

GMX_LIBGMX_EXPORT
void gmx_init_intranode_counters(t_commrec *cr);
/* Initializes intra-physical-node MPI process/thread counts and ID. */

gmx_bool gmx_mpi_initialized(void);
/* return TRUE when MPI_Init has been called.
 * return FALSE when MPI_Init has not been called OR
 * when GROMACS was compiled without MPI support.
 */

GMX_LIBGMX_EXPORT
void gmx_barrier(const t_commrec *cr);
/* Wait till all processes in cr->mpi_comm_mygroup have reached the barrier */

GMX_LIBGMX_EXPORT
void gmx_bcast(int nbytes, void *b, const t_commrec *cr);
/* Broadcast nbytes bytes from the master to cr->mpi_comm_mygroup */

GMX_LIBGMX_EXPORT
void gmx_bcast_sim(int nbytes, void *b, const t_commrec *cr);
/* Broadcast nbytes bytes from the sim master to cr->mpi_comm_mysim */

GMX_LIBGMX_EXPORT
void gmx_sumi(int nr, int r[], const t_commrec *cr);
/* Calculate the global sum of an array of ints */

void gmx_sumli(int nr, gmx_large_int_t r[], const t_commrec *cr);
/* Calculate the global sum of an array of large ints */

GMX_LIBGMX_EXPORT
void gmx_sumf(int nr, float r[], const t_commrec *cr);
/* Calculate the global sum of an array of floats */

GMX_LIBGMX_EXPORT
void gmx_sumd(int nr, double r[], const t_commrec *cr);
/* Calculate the global sum of an array of doubles */

void gmx_sumf_comm(int nr, float r[], MPI_Comm mpi_comm);
/* Calculate the global sum of an array of floats */

void gmx_sumd_comm(int nr, double r[], MPI_Comm mpi_comm);
/* Calculate the global sum of an array of doubles */

GMX_LIBGMX_EXPORT
void gmx_sumi_sim(int nr, int r[], const gmx_multisim_t *ms);
/* Calculate the sum over the simulations of an array of ints */

GMX_LIBGMX_EXPORT
void gmx_sumli_sim(int nr, gmx_large_int_t r[], const gmx_multisim_t *ms);
/* Calculate the sum over the simulations of an array of large ints */

GMX_LIBGMX_EXPORT
void gmx_sumf_sim(int nr, float r[], const gmx_multisim_t *ms);
/* Calculate the sum over the simulations of an array of floats */

GMX_LIBGMX_EXPORT
void gmx_sumd_sim(int nr, double r[], const gmx_multisim_t *ms);
/* Calculate the sum over the simulations of an array of doubles */

void gmx_abort(int nodeid, int nnodes, int errorno);
/* Abort the parallel run */

GMX_LIBGMX_EXPORT
void gmx_finalize_par(void);
/* Finish the parallel run in an ordered manner */

#ifdef GMX_DOUBLE
#define gmx_sum_comm  gmx_sumd_comm
#define gmx_sum       gmx_sumd
#define gmx_sum_sim   gmx_sumd_sim
#else
#define gmx_sum_comm  gmx_sumf_comm
#define gmx_sum       gmx_sumf
#define gmx_sum_sim   gmx_sumf_sim
#endif

#ifdef DEBUG_GMX
#define debug_gmx() do { FILE *fp = debug ? debug : stderr; \
                         if (bDebugMode()) { fprintf(fp, "NODEID=%d, %s  %d\n", gmx_mpi_initialized() ? gmx_node_rank() : -1, __FILE__, __LINE__); } fflush(fp); } while (0)
#else
#define debug_gmx()
#endif

#ifdef __cplusplus
}
#endif


#endif  /* _network_h */
