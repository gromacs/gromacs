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

#ifndef _partdec_h
#define _partdec_h
#include "visibility.h"
#include "vsite.h"

#ifdef __cplusplus
extern "C" {
#endif


#define GMX_LEFT     0          /* channel to the left processor  */
#define GMX_RIGHT    1          /* channel to the right processor */

/* These are the good old ring communication routines */

void gmx_tx(const t_commrec *cr, int dir, void *buf, int bufsize);
/*
 * Asynchronously sends bufsize bytes from the buffer pointed to by buf
 * over the communication channel, identified by chan. The buffer becomes
 * available after a successful call of gmx_tx_wait(dir).
 */

void gmx_tx_wait(const t_commrec *cr, int dir);
/*
 * Waits until the asynchronous send operation associated with chan has
 * succeeded. This makes the buffer of the send operation available to
 * the sending process.
 */

void gmx_rx(const t_commrec *cr, int dir, void *buf, int bufsize);
/*
 * Asynchronously receives bufsize bytes in the buffer pointed to by buf
 * from communication channel identified by chan. The buffer becomes
 * available after a successful call of gmx_rx_wait(chan).
 */

void gmx_rx_wait(const t_commrec *cr, int dir);
/*
 * Waits until the asynchronous receive operation, associated with chan,
 * has succeeded. This makes the buffer of the receive operation
 * available to the receiving process.
 */

void gmx_left_right(int nnodes, int nodeid,
                    int *left, int *right);
/* Get left and right proc id. */

void gmx_tx_rx(const t_commrec *cr,
               int send_dir, void *send_buf, int send_bufsize,
               int recv_dir, void *recv_buf, int recv_bufsize);
/* Communicate simultaneously left and right */

void gmx_tx_rx_real(const t_commrec *cr,
                    int send_dir, real *send_buf, int send_bufsize,
                    int recv_dir, real *recv_buf, int recv_bufsize);
/* Communicate simultaneously left and right, reals only */

void gmx_wait(const t_commrec *cr, int dir_send, int dir_recv);
/* Wait for communication to finish */

void pd_move_f(const t_commrec *cr, rvec f[], t_nrnb *nrnb);
/* Sum the forces over the nodes */

int *pd_cgindex(const t_commrec *cr);

int *pd_index(const t_commrec *cr);

int pd_shift(const t_commrec *cr);

int pd_bshift(const t_commrec *cr);

GMX_LIBMD_EXPORT
void pd_cg_range(const t_commrec *cr, int *cg0, int *cg1);
/* Get the range for the home charge groups */

GMX_LIBMD_EXPORT
void pd_at_range(const t_commrec *cr, int *at0, int *at1);
/* Get the range for the home particles */

GMX_LIBMD_EXPORT
gmx_localtop_t *split_system(FILE *log,
                             gmx_mtop_t *mtop, t_inputrec *inputrec,
                             t_commrec *cr);
/* Split the system over N processors. */

gmx_bool setup_parallel_vsites(t_idef *idef, t_commrec *cr,
                               t_comm_vsites *vsitecomm);

GMX_LIBMD_EXPORT
t_state *partdec_init_local_state(t_commrec *cr, t_state *state_global);
/* Generate a local state struct from the global one */

void
pd_get_constraint_range(gmx_partdec_p_t pd, int *start, int *natoms);


int *
pd_constraints_nlocalatoms(gmx_partdec_p_t pd);

/* Move x0 and also x1 if x1!=NULL */
void
pd_move_x_constraints(t_commrec *  cr,
                      rvec *       x0,
                      rvec *       x1);

#ifdef __cplusplus
}
#endif

#endif
