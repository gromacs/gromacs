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

#ifndef _cm5io_h
#define _cm5io_h

static char *SRCID_cm5io_h = "$Id$";

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>

extern void cm5io_tx(int pid,void *buf,int bufsize);
extern void cm5io_tx_wait(int pid);
extern void cm5io_txs(int pid,void *buf,int bufsize);
extern void cm5io_rx(int pid,void *buf,int bufsize);
extern void cm5io_rx_wait(int pid);
extern void cm5io_rxs(int pid,void *buf,int bufsize);
extern void cm5io_wait(int send,int receive);
extern void cm5io_init(int pid,int nprocs);
extern void cm5io_stat(FILE *fp,char *msg);
extern void cm5io_tx_rx(int send_pid,void *send_buf,int send_bufsize,
			int rec_pid,void *rec_buf,int rec_bufsize);
extern int  cm5_idle_send();
extern int  cm5_idle_rec();
extern void cm5_left_right(int nprocs,int pid,int *left,int *right);

#define gmx_tx		cm5io_tx
#define gmx_tx_wait	cm5io_tx_wait
#define gmx_txs		cm5io_txs
#define gmx_rx		cm5io_rx
#define gmx_rx_wait	cm5io_rx_wait
#define gmx_rxs		cm5io_rxs
#define gmx_init	cm5io_init
#define gmx_stat	cm5io_stat
#define gmx_tx_rx   	cm5io_tx_rx
#define gmx_wait    	cm5io_wait
#define gmx_sync_ring   cm5_sync_ring
#define gmx_node_num    cm5_node_num
#define gmx_node_id     cm5_node_id
#define gmx_idle_send   cm5_idle_send
#define gmx_idle_rec    cm5_idle_rec
#define gmx_left_right  cm5_left_right

#endif

