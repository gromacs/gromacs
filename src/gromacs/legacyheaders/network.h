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

#ifndef _network_h
#define _network_h


/*
 * This module defines the interface of the actual communication routines.
 */

#include <stdio.h>

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/fatalerror.h"

#ifdef __cplusplus
extern "C" {
#endif

struct gmx_domdec_t;
struct gmx_multisim_t;
struct t_commrec;

struct t_commrec *init_commrec(void);
/* Allocate, initialize and return the commrec. */

struct t_commrec *reinitialize_commrec_for_this_thread(const struct t_commrec *cro);
/* Initialize communication records for thread-parallel simulations.
   Must be called on all threads before any communication takes place by
   the individual threads. Copies the original commrec to
   thread-local versions (a small memory leak results because we don't
   deallocate the old shared version).  */

void gmx_fill_commrec_from_mpi(struct t_commrec *cr);
/* Continues t_commrec construction */

void gmx_setup_nodecomm(FILE *fplog, struct t_commrec *cr);
/* Sets up fast global communication for clusters with multi-core nodes */

void gmx_init_intranode_counters(struct t_commrec *cr);
/* Initializes intra-physical-node MPI process/thread counts and ID. */

void gmx_barrier(const struct t_commrec *cr);
/* Wait till all processes in cr->mpi_comm_mygroup have reached the barrier */

void gmx_bcast(int nbytes, void *b, const struct t_commrec *cr);
/* Broadcast nbytes bytes from the master to cr->mpi_comm_mygroup */

void gmx_bcast_sim(int nbytes, void *b, const struct t_commrec *cr);
/* Broadcast nbytes bytes from the sim master to cr->mpi_comm_mysim */

void gmx_sumi(int nr, int r[], const struct t_commrec *cr);
/* Calculate the global sum of an array of ints */

void gmx_sumli(int nr, gmx_int64_t r[], const struct t_commrec *cr);
/* Calculate the global sum of an array of large ints */

void gmx_sumf(int nr, float r[], const struct t_commrec *cr);
/* Calculate the global sum of an array of floats */

void gmx_sumd(int nr, double r[], const struct t_commrec *cr);
/* Calculate the global sum of an array of doubles */

void gmx_sumi_sim(int nr, int r[], const struct gmx_multisim_t *ms);
/* Calculate the sum over the simulations of an array of ints */

void gmx_sumli_sim(int nr, gmx_int64_t r[], const struct gmx_multisim_t *ms);
/* Calculate the sum over the simulations of an array of large ints */

void gmx_sumf_sim(int nr, float r[], const struct gmx_multisim_t *ms);
/* Calculate the sum over the simulations of an array of floats */

void gmx_sumd_sim(int nr, double r[], const struct gmx_multisim_t *ms);
/* Calculate the sum over the simulations of an array of doubles */

#ifdef GMX_DOUBLE
#define gmx_sum       gmx_sumd
#define gmx_sum_sim   gmx_sumd_sim
#else
#define gmx_sum       gmx_sumf
#define gmx_sum_sim   gmx_sumf_sim
#endif

gmx_bool gmx_fexist_master(const char *fname, struct t_commrec *cr);
/* Return TRUE when fname exists, FALSE otherwise, bcast from master to others */

void
gmx_fatal_collective(int f_errno, const char *file, int line,
                     const struct t_commrec *cr, struct gmx_domdec_t *dd,
                     const char *fmt, ...);
/* As gmx_fatal declared in utility/fatalerror.h,
 * but only the master process prints the error message.
 * This should only be called one of the following two situations:
 * 1) On all nodes in cr->mpi_comm_mysim, with cr!=NULL,dd==NULL.
 * 2) On all nodes in dd->mpi_comm_all,   with cr==NULL,dd!=NULL.
 * This will call MPI_Finalize instead of MPI_Abort when possible,
 * This is useful for handling errors in code that is executed identically
 * for all processes.
 */

/* This doesn't currently work if enabled (needs some header cleanup). */
#ifdef DEBUG_GMX
#define debug_gmx() do { FILE *fp = debug ? debug : stderr; \
                         if (bDebugMode()) { fprintf(fp, "rank=%d, %s  %d\n", gmx_mpi_initialized() ? gmx_node_rank() : -1, __FILE__, __LINE__); } fflush(fp); } while (0)
#else
#define debug_gmx()
#endif

#ifdef __cplusplus
}
#endif


#endif  /* _network_h */
