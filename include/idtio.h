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
 * Green Red Orange Magenta Azure Cyan Skyblue
 */

#ifndef _idtio_h
#define _idtio_h

static char *SRCID_idtio_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) idtio.h 1.9 11/23/92"
#endif /* HAVE_IDENT */

#include "rdklib.h"


extern void idtio_tx(int pid,void *buf,int bufsize);
extern void idtio_tx_wait(int pid);
extern void idtio_txs(int pid,void *buf,int bufsize);
extern void idtio_rx(int pid,void *buf,int bufsize);
extern void idtio_rx_wait(int pid);
extern void idtio_rxs(int pid,void *buf,int bufsize);
extern void idtio_init(int pid,int nprocs);
extern void idtio_stat(FILE *fp,char *msg);
extern void idt_reset_idle();
extern void idt_left_right(int nprocs,int pid,int *left,int *right);

#define gmx_tx		idtio_tx
#define gmx_tx_wait	idtio_tx_wait
#define gmx_txs		idtio_txs
#define gmx_rx		idtio_rx
#define gmx_rx_wait	idtio_rx_wait
#define gmx_rxs		idtio_rxs
#define gmx_init	idtio_init
#define gmx_stat	idtio_stat
#define gmx_cpu_num  	rdspccount
#define gmx_cpu_id   	rdspcnumber
#define gmx_left_right  idt_left_right
#define gmx_reset_idle  idt_reset_idle

#endif	/* _idtio_h */
