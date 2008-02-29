/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.3.3
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
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
 * Groningen Machine for Chemical Simulation
 */
/* This file is completely threadsafe - keep it that way! */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <sysstuff.h>
#include <string.h>
#include "typedefs.h"
#include "main.h"
#include "mvdata.h"
#include "network.h"
#include "smalloc.h"
#include "gmx_fatal.h"
#include "symtab.h"
#include "main.h"
#include "typedefs.h"
#include "vec.h"
#include "tgroup.h"
#include "block_tx.h"
#include "nrnb.h"

void move_rvecs(FILE *log,bool bForward,bool bSum,
		int left,int right,rvec vecs[],rvec buf[],
		int shift,t_nsborder *nsb,t_nrnb *nrnb)
{
  int    i,j,j0=137,j1=391;
  int    cur,nsum;
#define next ((cur+1) % nsb->nnodes)
#define prev ((cur-1+nsb->nnodes) % nsb->nnodes)

  if (bSum)
    cur=(nsb->nodeid+nsb->shift) % nsb->nnodes;
  else
    cur=nsb->nodeid;

  nsum=0;
  for(i=0; (i<shift); i++) {
    if (bSum) {
      if (bForward) {
	j0=nsb->index[prev];
	j1=j0+nsb->homenr[prev];
      }
      else {
	j0=nsb->index[next];
	j1=j0+nsb->homenr[next];
      }
      for(j=j0; (j<j1); j++) {
	clear_rvec(buf[j]);
      }
    }
    /* Forward pulse around the ring, to increasing NODE number */
    if (bForward) {
      if (bSum)
	gmx_tx_rx_real(right,vecs[nsb->index[cur]], nsb->homenr[cur]*DIM,
		       left,buf [nsb->index[prev]],nsb->homenr[prev]*DIM);
      else
	gmx_tx_rx_real(right,vecs[nsb->index[cur]], nsb->homenr[cur]*DIM,
		       left, vecs[nsb->index[prev]],nsb->homenr[prev]*DIM);
      /* Wait for communication to end */
      gmx_wait(right,left);
    }
    
    /* Backward pulse around the ring, to decreasing NODE number */
    else {
      if (bSum)
	gmx_tx_rx_real(left, vecs[nsb->index[cur]], nsb->homenr[cur]*DIM,
		       right,buf [nsb->index[next]],nsb->homenr[next]*DIM);
      else
	gmx_tx_rx_real(left, vecs[nsb->index[cur]], nsb->homenr[cur]*DIM,
		       right,vecs[nsb->index[next]],nsb->homenr[next]*DIM);
      /* Wait for communication to end */
      gmx_wait(left,right);
    }

    /* Actual summation */
    if (bSum) {
      for(j=j0; (j<j1); j++) {
	rvec_inc(vecs[j],buf[j]);
      }
      nsum+=(j1-j0);
    }
    if (bForward) 
      cur=prev;
    else
      cur=next;
  }  
  if (nsum > 0)
    inc_nrnb(nrnb,eNR_FSUM,nsum);
#undef next
#undef prev
}

void move_x(FILE *log,
	    int left,int right,rvec x[],t_nsborder *nsb,
	    t_nrnb *nrnb)
{
  move_rvecs(log,FALSE,FALSE,left,right,x,NULL,nsb->shift,nsb,nrnb);
  move_rvecs(log,TRUE, FALSE,left,right,x,NULL,nsb->bshift,nsb,nrnb);

  where();
}

void move_f(FILE *log,
	    int left,int right,rvec f[],rvec fadd[],
	    t_nsborder *nsb,t_nrnb *nrnb)
{
  move_rvecs(log,TRUE, TRUE,left,right,f,fadd,nsb->shift,nsb,nrnb);
  move_rvecs(log,FALSE,TRUE,left,right,f,fadd,nsb->bshift,nsb,nrnb);

  where();
}

void move_cgcm(FILE *log,t_commrec *cr,rvec cg_cm[],int nload[])
{
  int i,start,nr;
  int cur=cr->nodeid;
#define next ((cur+1) % cr->nnodes)
  
  for(i=0; (i<cr->nnodes-1); i++) {
    start = (cur == 0) ? 0 : nload[cur-1];
    nr    = nload[cur] - start;
    gmx_tx(cr->left, cg_cm[start], nr*sizeof(cg_cm[0]));
#ifdef DEBUG
    fprintf(log,"move_cgcm: TX start=%d, nr=%d\n",start,nr);
#endif    
    start = (next == 0) ? 0 : nload[next-1];
    nr    = nload[next] - start;
    gmx_rx(cr->right,cg_cm[start], nr*sizeof(cg_cm[0]));
#ifdef DEBUG
    fprintf(log,"move_cgcm: RX start=%d, nr=%d\n",start,nr);
#endif    
    gmx_tx_wait(cr->left);
    gmx_rx_wait(cr->right);
    
    cur=next;
  }
#undef next
}


