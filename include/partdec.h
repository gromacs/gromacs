/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Gromacs Runs On Most of All Computer Systems
 */

#ifndef _partdec_h
#define _partdec_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "vsite.h"

#define GMX_LEFT     0          /* channel to the left processor  */
#define GMX_RIGHT    1          /* channel to the right processor */

/* These are the good old ring communication routines */

extern void gmx_tx(const t_commrec *cr,int dir,void *buf,int bufsize);
     /*
      * Asynchronously sends bufsize bytes from the buffer pointed to by buf 
      * over the communication channel, identified by chan. The buffer becomes 
      * available after a successful call of gmx_tx_wait(dir).
      */

extern void gmx_tx_wait(int dir);
     /*
      * Waits until the asynchronous send operation associated with chan has 
      * succeeded. This makes the buffer of the send operation available to 
      * the sending process.
      */

extern void gmx_rx(const t_commrec *cr,int dir,void *buf,int bufsize);
     /*
      * Asynchronously receives bufsize bytes in the buffer pointed to by buf 
      * from communication channel identified by chan. The buffer becomes 
      * available after a successful call of gmx_rx_wait(chan).
      */

extern void gmx_rx_wait(int dir);
     /*
      * Waits until the asynchronous receive operation, associated with chan, 
      * has succeeded. This makes the buffer of the receive operation 
      * available to the receiving process.
      */

extern void gmx_left_right(int nnodes,int nodeid,
			   int *left,int *right);
/* Get left and right proc id. */

extern void gmx_tx_rx(const t_commrec *cr,
		      int send_dir,void *send_buf,int send_bufsize,
		      int recv_dir,void *recv_buf,int recv_bufsize);
/* Communicate simultaneously left and right */
		      
extern void gmx_tx_rx_real(const t_commrec *cr,
			   int send_dir,real *send_buf,int send_bufsize,
			   int recv_dir,real *recv_buf,int recv_bufsize);
/* Communicate simultaneously left and right, reals only */

extern void gmx_wait(int dir_send,int dir_recv);
/* Wait for communication to finish */

extern void pd_move_f(const t_commrec *cr,rvec f[],t_nrnb *nrnb);
/* Sum the forces over the nodes */

extern int *pd_cgindex(const t_commrec *cr);

extern int *pd_index(const t_commrec *cr);

extern int pd_shift(const t_commrec *cr);

extern int pd_bshift(const t_commrec *cr);

extern void pd_cg_range(const t_commrec *cr,int *cg0,int *cg1);
/* Get the range for the home charge groups */

extern void pd_at_range(const t_commrec *cr,int *at0,int *at1);
/* Get the range for the home particles */

extern gmx_localtop_t *split_system(FILE *log,
				    gmx_mtop_t *mtop,t_inputrec *inputrec,
				    t_commrec *cr);
/* Split the system over N processors. */

extern bool setup_parallel_vsites(t_idef *idef,t_commrec *cr,
				  t_comm_vsites *vsitecomm);

extern t_state *partdec_init_local_state(t_commrec *cr,t_state *state_global);
/* Generate a local state struct from the global one */

#endif
