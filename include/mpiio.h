/*
 * $Id$
 * 
 *       This source code is part of
 * 
 *        G   R   O   M   A   C   S
 * 
 * GROningen MAchine for Chemical Simulations
 * 
 *               VERSION 2.0
 * 
 * Copyright (c) 1991-1999
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 * 
 * Also check out our WWW page:
 * http://md.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 * 
 * And Hey:
 * Good ROcking Metal Altar for Chronical Sinners
 */

#ifndef _mpiio_h
#define _mpiio_h

static char *SRCID_mpiio_h = "$Id$";

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <mpi.h>
#include "typedefs.h"

extern void mpiio_tx(int nodeid,void *buf,int bufsize);
extern void mpiio_tx_wait(int nodeid);
extern void mpiio_txs(int nodeid,void *buf,int bufsize);
extern void mpiio_rx(int nodeid,void *buf,int bufsize);
extern void mpiio_rx_wait(int nodeid);
extern void mpiio_rxs(int nodeid,void *buf,int bufsize);
extern int  mpiio_setup(int *argc,char *argv[],int *nnodes);
/* If this is the master, spawn some kids. Nnodes is set to the
 * number of nodeessors.
 * Return nodeid */
 
extern int mpiio_rx_probe(int nodeid);
/* Check whether message arrived, if so message is read and return 0 */

extern void mpiio_stat(FILE *fp,char *msg);
extern int  mpinodenumber(void);
extern int  mpinodecount(void);
extern int  mpi_idle_send(void);
extern int  mpi_idle_rec(void);
extern void mpi_left_right(int nnodes,int nodeid,int *left,int *right);
extern void mpiio_tx_rx(int send_nodeid,void *send_buf,int send_bufsize,
			int rec_nodeid,void *rec_buf,int rec_bufsize);
extern void mpiio_wait(int left,int right);
extern void mpiio_sync_ring(int nodeid,int nnodes,int left,int right);
extern void mpi_reset_idle(void);
extern void mpi_abort(int nodeid,int nnodes,int errorno);

#define gmx_tx       	mpiio_tx
#define gmx_tx_wait  	mpiio_tx_wait
#define gmx_txs      	mpiio_txs
#define gmx_rx       	mpiio_rx
#define gmx_rx_wait  	mpiio_rx_wait
#define gmx_rxs      	mpiio_rxs
#define gmx_stat     	mpiio_stat
#define gmx_wait       	mpiio_wait
#define gmx_node_num   	mpinodecount
#define gmx_node_id    	mpinodenumber
#define gmx_left_right 	mpi_left_right
#define gmx_idle_send   mpi_idle_send
#define gmx_idle_rec  	mpi_idle_rec
#define gmx_reset_idle  mpi_reset_idle
#define gmx_abort       mpi_abort

#endif 
