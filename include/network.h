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
 * Grunge ROck MAChoS
 */

#ifndef _network_h
#define _network_h

static char *SRCID_network_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) network.h 1.9 11/23/92"
#endif /* HAVE_IDENT */

/*
 * This module defines the interface of the actual communication routines.
 */

#include <stdio.h>
#include "typedefs.h"

#define LEFT     0          /* channel to the left processor  */
#define RIGHT    1          /* channel to the right processor */

#define record(rec)     &(rec),sizeof(rec)
#define array(arr,nr)   (arr),((nr)*sizeof((arr)[0]))
#define arrayp(el,nr)   &(el),((nr)*sizeof(el))
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
 * For sequential processing dummies are in Kernel/sys/libnet.c
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

extern int gmx_cpu_num(void);
/* return the number of cpus in the ring */

extern int gmx_cpu_id(void);
/* return the identification ID of the cpu */
      
extern void gmx_left_right(int nprocs,int pid,int *left,int *right);
/* Get left and right proc id. */

/****************************************************** 
 *
 * Here are some routines that are platform independent.
 * In principle they might be overridden by more efficient
 * routines for particular library packages (pvm, mpi?)
 *
 ******************************************************/
extern void gmx_stat(FILE *fp,char *msg);
/* Prints a overview of the status of the network, useful for debugging. */

extern void gmx_reset_idle(void);
/* Reset the idle count */

extern void gmx_tx_rx(int send_pid,void *send_buf,int send_bufsize,
		      int rec_pid,void *rec_buf,int rec_bufsize);

extern void gmx_wait(int send,int receive);

extern void gmx_sync_ring(int pid,int nprocs,int left,int right);
/* Synchronise the ring... */

extern void gmx_sumi(int nr,int r[],t_commrec *cr);
/* Calculate the global sum of an array of ints */

extern void gmx_sumf(int nr,float r[],t_commrec *cr);
/* Calculate the global sum of an array of floats */

extern void gmx_sumd(int nr,double r[],t_commrec *cr);
/* Calculate the global sum of an array of doubles */

/******************************************************
 *
 * These routines are now superseded by a macro...
 * Each of the communication libraries may override the
 * macros, hopefully the compiler will tell you!
 *
 ******************************************************/

#ifdef USE_IDT
#include "idtio.h"
#endif

#ifdef USE_CM5
#include "cm5io.h"
#endif

#ifdef USE_PVM3
#include "pvmio.h"
#endif

#ifdef USE_MPI
#include "mpiio.h"
#endif

#ifdef USE_AMB
#include "ambio.h"
#endif

/********************************************************
 *
 * Some routines, do not have an implementation everywhere,
 * for these there are defaults using our low-level routines.
 *
 *******************************************************/
#ifndef gmx_wait
extern  void def_wait(int send,int receive);
#define gmx_wait def_wait
#endif

#ifndef gmx_tx_rx
extern  void def_tx_rx(int send_pid,void *send_buf,int send_bufsize,
		       int rec_pid,void *rec_buf,int rec_bufsize);
#define gmx_tx_rx def_tx_rx
#endif 

#ifndef gmx_sync_ring
extern  void def_sync_ring(int pid,int nprocs,int left,int right);
#define gmx_sync_ring def_sync_ring
#endif

#ifndef gmx_stat
extern  void def_stat(FILE *fp,char *msg);
#define gmx_stat def_stat
#endif

#ifndef gmx_reset_idle
extern  void def_reset_idle(void);
#define gmx_reset_idle def_reset_idle
#endif

#ifndef gmx_sumf
extern  void def_sumf(int nr,float r[],t_commrec *cr);
#define gmx_sumf def_sumf
#endif

#ifndef gmx_sumd
extern  void def_sumd(int nr,double r[],t_commrec *cr);
#define gmx_sumd def_sumd
#endif

#ifndef gmx_sumi
extern  void def_sumi(int nr,int r[],t_commrec *cr);
#define gmx_sumi def_sumi
#endif

#ifdef DOUBLE
#define gmx_sum gmx_sumd
#else
#define gmx_sum gmx_sumf
#endif

#ifdef DEBUGPAR
#define debug_gmx() fprintf(stderr,"PID=%d, %s  %d\n",gmx_cpu_id(),__FILE__,__LINE__)
#else
#define debug_gmx()
#endif

#endif	/* _network_h */
