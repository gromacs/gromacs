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
 * Green Red Orange Magenta Azure Cyan Skyblue
 */
static char *SRCID_mvxvf_c = "$Id$";

#include <sysstuff.h>
#include <string.h>
#include "typedefs.h"
#include "main.h"
#include "assert.h"
#include "mvdata.h"
#include "network.h"
#include "smalloc.h"
#include "fatal.h"
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
#define next ((cur+1) % nsb->nprocs)
#define prev ((cur-1+nsb->nprocs) % nsb->nprocs)

  if (bSum)
    cur=(nsb->pid+nsb->shift) % nsb->nprocs;
  else
    cur=nsb->pid;

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
    /* Forward pulse around the ring, to increasing CPU number */
    if (bForward) {
      gmx_tx(right,  arrayp(vecs[nsb->index[cur]], nsb->homenr[cur]));
      /* If we want to sum these arrays, we have to store the rvecs
       * we receive in a temp buffer.
       */
      if (bSum)
	gmx_rx(left, arrayp(buf [nsb->index[prev]],nsb->homenr[prev]));
      else
	gmx_rx(left, arrayp(vecs[nsb->index[prev]],nsb->homenr[prev]));
    }
    
    /* Backward pulse around the ring, to decreasing CPU number */
    else {
      gmx_tx(left,    arrayp(vecs[nsb->index[cur]], nsb->homenr[cur]));
      /* See above */
      if (bSum)
	gmx_rx(right, arrayp(buf [nsb->index[next]],nsb->homenr[next]));
      else
	gmx_rx(right, arrayp(vecs[nsb->index[next]],nsb->homenr[next]));
    }
    /* Wait for communication to end */
    if (bForward) {
      gmx_tx_wait(right);
      gmx_rx_wait(left);
    }
    else {
      gmx_tx_wait(left);
      gmx_rx_wait(right);
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
#ifdef TEST_MPI_X
  static int *recvcounts=NULL,*displs;
  int i;

  if (recvcounts == NULL) {
    snew(recvcounts,nsb->nprocs);
    snew(displs,nsb->nprocs);
    for(i=0; i<nsb->nprocs; i++) {
      recvcounts[i] = nsb->homenr[i]*sizeof(x[0]);
      displs[i]     = nsb->index[i]*sizeof(x[0]);
    }
  }
  MPI_Allgatherv(arrayp(x[nsb->index[nsb->pid]],nsb->homenr[nsb->pid]),
		 MPI_BYTE,x,recvcounts,displs,MPI_BYTE,MPI_COMM_WORLD);
#else
  move_rvecs(log,FALSE,FALSE,left,right,x,NULL,nsb->shift,nsb,nrnb);
  where();
  move_rvecs(log,TRUE, FALSE,left,right,x,NULL,nsb->bshift,nsb,nrnb);
#endif
  where();
}

void move_f(FILE *log,
	    int left,int right,rvec f[],rvec fadd[],
	    t_nsborder *nsb,t_nrnb *nrnb)
{
#ifdef MPI_TEST_F
#ifdef DOUBLE
#define mpi_type MPI_DOUBLE
#else
#define mpi_type MPI_FLOAT
#endif
  MPI_Allreduce(f,f,nsb->natoms*DIM,mpi_type,MPI_SUM,MPI_COMM_WORLD);
#else
  move_rvecs(log,TRUE, TRUE,left,right,f,fadd,nsb->shift,nsb,nrnb);
  where();
  move_rvecs(log,FALSE,TRUE,left,right,f,fadd,nsb->bshift,nsb,nrnb);
#endif
  where();
}

void move_cgcm(FILE *log,t_commrec *cr,rvec cg_cm[],int nload[])
{
  int i,start,nr;
  int cur=cr->pid;
#define next ((cur+1) % cr->nprocs)
  
  for(i=0; (i<cr->nprocs-1); i++) {
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


