/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.0
 * 
 * Copyright (c) 1991-2001
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
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
 * Do check out http://www.gromacs.org , or mail us at gromacs@gromacs.org .
 * 
 * And Hey:
 * Giving Russians Opium May Alter Current Situation
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
