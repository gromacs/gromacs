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
 * GROningen Mixture of Alchemy and Childrens' Stories
 */

#ifndef _ambio_h
#define _ambio_h

static char *SRCID_ambio_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) ambio.h 1.11 10/15/97"
#endif /* HAVE_IDENT */

extern int cpuCount();
extern int cpuNumber();

extern void network_tx(int chan,void *buf,int bufsize);
extern void network_tx_wait(int chan);
extern void network_txs(int chan,void *buf,int bufsize);
extern void network_rx(int chan,void *buf,int bufsize);
extern void network_rx_wait(int chan);
extern void network_rxs(int chan,void *buf,int bufsize);
extern void network_init(int pid,int nprocs);
extern void amb_left_right(int nprocs,int pid,int *left,int *right);

#define gmx_tx		network_tx
#define gmx_tx_wait	network_tx_wait
#define gmx_txs		network_txs
#define gmx_rx		network_rx
#define gmx_rx_wait	network_rx_wait
#define gmx_rxs		network_rxs
#define gmx_cpu_num     cpuCount
#define gmx_cpu_id      cpuNumber
#define gmx_left_right  amb_left_right

#endif	/* _ambio_h */
