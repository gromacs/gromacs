/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0
 * 
 * Copyright (c) 1991-1997
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 *
 * Also check out our WWW page:
 * http://rugmd0.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 *
 * And Hey:
 * Gyas ROwers Mature At Cryogenic Speed
 */

#ifndef _pvmio_h
#define _pvmio_h

static char *SRCID_pvmio_h = "$Id$";

#include "pvm3.h"

extern void pvmio_tx(int pid,void *buf,int bufsize);
extern void pvmio_tx_wait(int pid);
extern void pvmio_txs(int pid,void *buf,int bufsize);
extern void pvmio_rx(int pid,void *buf,int bufsize);
extern void pvmio_rx_wait(int pid);
extern void pvmio_rxs(int pid,void *buf,int bufsize);
extern int  pvmio_setup(char *argv[],int nprocs);
/* If this is the master, spawn some kids. Return pid */

extern void pvmio_stat(FILE *fp,char *msg);
extern int  pvmnodenumber(void);
extern int  pvmnodecount(void);
extern int  pvm_idle_send(void);
extern int  pvm_idle_rec(void);
extern void pvm_left_right(int nprocs,int pid,int *left,int *right);
extern void pvmio_tx_rx(int send_pid,void *send_buf,int send_bufsize,
			int rec_pid,void *rec_buf,int rec_bufsize);
extern void pvmio_wait(int left,int right);
extern char *pvm_error(int errorno);
extern void pvmio_sync_ring(int pid,int nprocs,int left,int right);
extern void pvm_reset_idle();
extern void pvm_abort(int pid,int nprocs,int errno);

#define gmx_tx       	pvmio_tx
#define gmx_tx_wait  	pvmio_tx_wait
#define gmx_txs      	pvmio_txs
#define gmx_rx       	pvmio_rx
#define gmx_rx_wait  	pvmio_rx_wait
#define gmx_rxs      	pvmio_rxs
#define gmx_stat     	pvmio_stat
#define gmx_wait       	pvmio_wait
#define gmx_cpu_num    	pvmnodecount
#define gmx_cpu_id     	pvmnodenumber
#define gmx_left_right 	pvm_left_right
#define gmx_idle_send   pvm_idle_send
#define gmx_idle_rec  	pvm_idle_rec
#define gmx_reset_idle  pvm_reset_idle
#define gmx_abort       pvm_abort

#endif
