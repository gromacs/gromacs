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

#ifndef _network_h
#define _network_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


/*
 * This module defines the interface of the actual communication routines.
 */

#include <stdio.h>
#include "typedefs.h"
#include "main.h"
#include "fatal.h"

#define LEFT     0          /* channel to the left processor  */
#define RIGHT    1          /* channel to the right processor */

#define record(rec)     &((rec)),sizeof(rec)
#define array(arr,nr)   (arr),((nr)*sizeof((arr)[0]))
#define arrayp(el,nr)   &((el)),((nr)*sizeof(el))
/* 
 * These macro's can be used as shown in the following examples:
 *
 * int chan=1;
 * int nr;
 * struct {float x,y} coordinate;
 * int arr[10];
 *
 * gmx_rxs(chan,record(nr));		receive data in nr
 * gmx_txs(chan,record(coordinate));	sends data from coordinate
 * gmx_rxs(chan,array(arr,10));	sends an array of 10 elements
 * gmx_rxs(chan,arrayp(arr[3],4)); 	receives an array of 4 elements
 *					and stores it starting at element 3
 */

/****************************************************** 
 *
 * Here are the communication routines to be called from GROMACS
 * programs!
 *
 * The following 9 routines MUST be overridden !!!!!!!!
 * (for parallel processing)
 *
 * For sequential processing dummies are in src/gmxlib/libnet.c
 *
 ******************************************************/
extern void gmx_tx(int chan,void *buf,int bufsize);
     /*
      * Asynchronously sends bufsize bytes from the buffer pointed to by buf 
      * over the communication channel, identified by chan. The buffer becomes 
      * available after a successful call of gmx_tx_wait(chan).
      */

extern void gmx_tx_wait(int chan);
     /*
      * Waits until the asynchronous send operation associated with chan has 
      * succeeded. This makes the buffer of the send operation available to 
      * the sending process.
      */

extern void gmx_txs(int chan,void *buf,int bufsize);
     /*
      * Synchronously sends bufsize bytes from the buffer pointed to by buf to
      * the processor/process identified by chan. This is implemented by a call
      * to gmx_tx(chan,buf,bufsize), directly followed by a call to 
      * gmx_tx_wait(chan), so the buffer is available after 
      * gmx_txs() returns.
      */

extern void gmx_rx(int chan,void *buf,int bufsize);
     /*
      * Asynchronously receives bufsize bytes in the buffer pointed to by buf 
      * from communication channel identified by chan. The buffer becomes 
      * available after a successful call of gmx_rx_wait(chan).
      */

extern void gmx_rx_wait(int chan);
     /*
      * Waits until the asynchronous receive operation, associated with chan, 
      * has succeeded. This makes the buffer of the receive operation 
      * available to the receiving process.
      */

extern void gmx_rxs(int chan,void *buf,int bufsize);
     /*
      * Synchronously receives bufsize bytes from the buffer pointed to by 
      * buf over the communication channel identified by chan. This is 
      * implemented by a call to gmx_rx(chan,buf,bufsize), directly 
      * followed by a call to gmx_rx_wait(chan), so the buffer is 
      * available after gmx_rxs() returns.
      */

extern int gmx_setup(int *argc,char **argv,int *nnodes);
/* Initializes the parallel communication, return the ID of the node */

extern int gmx_node_num(void);
/* return the number of nodes in the ring */

extern int gmx_node_id(void);
/* return the identification ID of the node */
      
extern void gmx_left_right(int nnodes,int nodeid,int *left,int *right);
/* Get left and right proc id. */

extern void gmx_stat(FILE *fp,char *msg);
/* Prints a overview of the status of the network, useful for debugging. */

extern void gmx_reset_idle(void);
/* Reset the idle count */

extern void gmx_tx_rx(int send_nodeid,void *send_buf,int send_bufsize,
		      int rec_nodeid,void *rec_buf,int rec_bufsize);
/* Communicate simultaneously left and right */
		      
extern void gmx_tx_rx_real(int send_nodeid,real *send_buf,int send_bufsize,
			   int rec_nodeid,real *rec_buf,int rec_bufsize);
/* Communicate simultaneously left and right, reals only */

extern void gmx_wait(int send,int receive);
/* Wait for communication to finish */

extern void gmx_sync_ring(int nodeid,int nnodes,int left,int right);
/* Synchronise the ring... */

extern void gmx_sumi(int nr,int r[],const t_commrec *cr);
/* Calculate the global sum of an array of ints */

extern void gmx_sumf(int nr,float r[],const t_commrec *cr);
/* Calculate the global sum of an array of floats */

extern void gmx_sumd(int nr,double r[],const t_commrec *cr);
/* Calculate the global sum of an array of doubles */

extern void gmx_abort(int nodeid,int nnodes,int errorno);
/* Abort the parallel run */

extern void gmx_finalize(t_commrec *cr);
/* Finish the parallel run in an ordered manner */

#ifdef DOUBLE
#define gmx_sum gmx_sumd
#else
#define gmx_sum gmx_sumf
#endif

#ifdef DEBUG_GMX
#define debug_gmx() do { FILE *fp=debug ? debug : (stdlog ? stdlog : stderr);\
if (bDebugMode()) fprintf(fp,"NODEID=%d, %s  %d\n",gmx_node_id(),__FILE__,__LINE__); fflush(fp); } while (0)
#else
#define debug_gmx()
#endif

#endif	/* _network_h */
