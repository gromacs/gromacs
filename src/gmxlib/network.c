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
static char *SRCID_network_c = "$Id$";

#include "typedefs.h"
#include "smalloc.h"
#include "network.h"

void def_sync_ring(int pid,int nprocs,int left,int right)
{
  int i;
  int tag=0;

  for (i=0; (i<nprocs); i++) {
    if (pid == 0) {
      gmx_txs(right,record(tag));
      gmx_rxs(left,record(tag));
    }
    else {
      gmx_rxs(left,record(tag));
      gmx_txs(right,record(tag));
    }
  }
}

void def_tx_rx(int send_pid,void *send_buf,int send_bufsize,
	       int rec_pid,void *rec_buf,int rec_bufsize)
{
  gmx_tx(send_pid,send_buf,send_bufsize);
  gmx_rx(rec_pid,rec_buf,rec_bufsize);
  gmx_tx_wait(send_pid);
  gmx_rx_wait(rec_pid);
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
  double *buf[2];
  int  NR,bufs,j,i,cur=0;
#define next (1-cur)

  bufs=nr*sizeof(buf[0][0]);
#ifdef _amb_
  bufs+=23;
  bufs=bufs-(bufs % 24);
  NR=bufs/sizeof(double);
#else
  NR=nr;
#endif

  snew(buf[0],NR);
  snew(buf[1],NR);

  for(i=0; (i<nr); i++)
    buf[cur][i]=r[i];
  for(j=0; (j<cr->nprocs-1); j++) {
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
  for(j=0; (j<cr->nprocs-1); j++) {
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
  for(j=0; (j<cr->nprocs-1); j++) {
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
	     int nprocs,int pid,int left,int right)
{
#define MAXBUF 2737
  real      buf[MAXBUF];
  int       i;
  t_commrec cr;

  for(i=0; (i<MAXBUF); i++)
    buf[i]=i;
  cr.nprocs=nprocs;
  cr.pid=pid;
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
