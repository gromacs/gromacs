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
 * Gromacs Runs On Most of All Computer Systems
 */
static char *SRCID_network_c = "$Id$";
#include "typedefs.h"
#include "smalloc.h"
#include "network.h"

void def_sync_ring(int nodeid,int nnodes,int left,int right)
{
  int i;
  int tag=0;

  for (i=0; (i<nnodes); i++) {
    if (nodeid == 0) {
      gmx_txs(right,&tag,sizeof(tag));
      gmx_rxs(left,&tag,sizeof(tag));
    }
    else {
      gmx_rxs(left,&tag,sizeof(tag));
      gmx_txs(right,&tag,sizeof(tag));
    }
  }
}

void def_tx_rx(int send_nodeid,void *send_buf,int send_bufsize,
	       int rec_nodeid,void *rec_buf,int rec_bufsize)
{
  gmx_tx(send_nodeid,send_buf,send_bufsize);
  gmx_rx(rec_nodeid,rec_buf,rec_bufsize);
  gmx_tx_wait(send_nodeid);
  gmx_rx_wait(rec_nodeid);
}

void def_wait(int send,int receive)
{
  gmx_tx_wait(send);
  gmx_rx_wait(receive);
}

void def_stat(FILE *fp,char *msg)
{
  fprintf(fp,"def_stat: %s (from %s, %d)\n",msg,__FILE__,__LINE__);
}

void def_reset_idle(void)
{
}

void def_sumd(int nr,double r[],t_commrec *cr)
{
  /*#define TEST_MPI_SUM*/
#ifdef TEST_MPI_SUM
  static double *buf;
  static int nalloc=0;
  int i;
  
  if (nr > nalloc) {
    nalloc = nr;
    srenew(buf,nalloc);
  }
  MPI_Allreduce(r,buf,nr,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  for(i=0; i<nr; i++)
    r[i] = buf[i];
#else
  double *buf[2];
  int  NR,bufs,j,i,cur=0;
#define next (1-cur)

  bufs=nr*sizeof(buf[0][0]);
  NR=nr;

  snew(buf[0],NR);
  snew(buf[1],NR);

  for(i=0; (i<nr); i++)
    buf[cur][i]=r[i];
  for(j=0; (j<cr->nnodes-1); j++) {
    gmx_tx(cr->left,buf[cur],bufs);
    gmx_rx(cr->right,buf[next],bufs);
    gmx_tx_wait(cr->left);
    gmx_rx_wait(cr->right);
    for(i=0; (i<nr); i++)
      r[i]+=buf[next][i];

    cur=next;
  }
  sfree(buf[1]);
  sfree(buf[0]);
#endif
}

void def_sumf(int nr,float r[],t_commrec *cr)
{
  float *buf[2];
  int  NR,bufs,j,i,cur=0;
#define next (1-cur)

  bufs=nr*sizeof(float);
#ifdef _amb_
  bufs+=23;
  bufs=bufs-(bufs % 24);
  NR=bufs/sizeof(float);
#else
  NR=nr;
#endif

  snew(buf[0],NR);
  snew(buf[1],NR);

  for(i=0; (i<nr); i++)
    buf[cur][i]=r[i];
  for(j=0; (j<cr->nnodes-1); j++) {
    gmx_tx(cr->left,buf[cur],bufs);
    gmx_rx(cr->right,buf[next],bufs);
    gmx_wait(cr->left,cr->right);
    for(i=0; (i<nr); i++)
      r[i]+=buf[next][i];

    cur=next;
  }
  sfree(buf[1]);
  sfree(buf[0]);
}

void def_sumi(int nr,int r[],t_commrec *cr)
{
  int *buf[2];
  int  NR,bufs,j,i,cur=0;
#define next (1-cur)

  bufs=nr*sizeof(int);
#ifdef _amb_
  bufs+=23;
  bufs=bufs-(bufs % 24);
  NR=bufs/sizeof(int);
#else
  NR=nr;
#endif

  snew(buf[0],NR);
  snew(buf[1],NR);

  for(i=0; (i<nr); i++)
    buf[cur][i]=r[i];
  for(j=0; (j<cr->nnodes-1); j++) {
    gmx_tx(cr->left,buf[cur],bufs);
    gmx_rx(cr->right,buf[next],bufs);
    gmx_wait(cr->left,cr->right);
    for(i=0; (i<nr); i++)
      r[i]+=buf[next][i];

    cur=next;
  }
  sfree(buf[1]);
  sfree(buf[0]);
}

#ifdef GSUM
#include "main.h"

int i860main(int argc,char *argv[],FILE *log,
	     int nnodes,int nodeid,int left,int right)
{
#define MAXBUF 2737
  real      buf[MAXBUF];
  int       i;
  t_commrec cr;

  for(i=0; (i<MAXBUF); i++)
    buf[i]=i;
  cr.nnodes=nnodes;
  cr.nodeid=nodeid;
  cr.left=left;
  cr.right=right;

  gmx_sum(MAXBUF,buf,&cr);
  
  /*
  for(i=0; (i<MAXBUF); i++)
    fprintf(log,"buf[%2d]=%g\n",i,buf[i]);
    */
  return 0;
}
#endif
