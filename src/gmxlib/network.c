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
 * GROningen Mixture of Alchemy and Childrens' Stories
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string.h>
#include "fatal.h"
#include "main.h"
#include "smalloc.h"
#include "network.h"
#ifdef USE_MPI
#include <mpi.h>
#endif

#ifdef USE_MPI

static MPI_Request mpi_req_tx=MPI_REQUEST_NULL,mpi_req_rx;
#else
#define MYFATAL(str) fatal_error(0,"Routine %s called in %s, %d",str,__FILE__,__LINE__)
#endif

/* Try setting MPI_TEST when you experience unexplainable crashes, *
 * up til now these crashes have only occured with IRIX 6.5        */
/* #define MPI_TEST */

void gmx_tx(int nodeid,void *buf,int bufsize)
{
#ifndef USE_MPI
  MYFATAL("gmx_tx"); 
#else
  int        tag,flag;
  MPI_Status status;
  
#ifdef DEBUG
  fprintf(stderr,"gmx_tx: nodeid=%d, buf=%x, bufsize=%d\n",
	  nodeid,buf,bufsize);
#endif
#ifdef MPI_TEST
  /* workaround for crashes encountered with MPI on IRIX 6.5 */
  if (mpi_req_tx != MPI_REQUEST_NULL) {
    MPI_Test(&mpi_req_tx,&flag,&status);
    if (flag==FALSE) {
      fprintf(stdlog,"gmx_tx called before previous send was complete: nodeid=%d, buf=%x, bufsize=%d\n",
	      nodeid,buf,bufsize);
      gmx_tx_wait(nodeid);
    }
  }
#endif
  tag = 0;
  if (MPI_Isend(buf,bufsize,MPI_BYTE,nodeid,tag,MPI_COMM_WORLD,&mpi_req_tx) != 0)
    fatal_error(0,"MPI_Isend Failed !");
#endif
}

void gmx_tx_wait(int nodeid)
{
#ifndef USE_MPI
  MYFATAL("gmx_tx_wait");
#else
  MPI_Status  status;
  int mpi_result;
  
  if ((mpi_result=MPI_Wait(&mpi_req_tx,&status)) != 0)
    fatal_error(0,"MPI_Wait: result=%d",mpi_result);
#endif
}

void gmx_txs(int nodeid,void *buf,int bufsize)
{
#ifndef USE_MPI
  MYFATAL("gmx_txs");
#else
  int tag;

#ifdef DEBUG
  fprintf(stderr,"gmx_txs: nodeid=%d, buf=%x, bufsize=%d\n",
	  nodeid,buf,bufsize);
#endif
  tag = 0;
  if (MPI_Send(buf,bufsize,MPI_BYTE,nodeid,tag,MPI_COMM_WORLD) != 0)
    fatal_error(0,"MPI_Send Failed !");
#endif
}

void gmx_rx(int nodeid,void *buf,int bufsize)
{
#ifndef USE_MPI
  MYFATAL("gmx_rx");
#else
  int        tag;

#ifdef DEBUG
  fprintf(stderr,"gmx_rx: nodeid=%d, buf=%x, bufsize=%d\n",
	  nodeid,buf,bufsize);
#endif
  tag = 0;
  if (MPI_Irecv( buf, bufsize, MPI_BYTE, nodeid, tag, MPI_COMM_WORLD, &mpi_req_rx) != 0 )
    fatal_error(0,"MPI_Recv Failed !");
#endif
}

void gmx_rx_wait(int nodeid)
{
#ifndef USE_MPI
  MYFATAL("gmx_rx_wait");
#else
  MPI_Status  status;
  int mpi_result;
  
  if ((mpi_result=MPI_Wait(&mpi_req_rx,&status)) != 0)
    fatal_error(0,"MPI_Wait: result=%d",mpi_result);
#endif
}

int gmx_rx_probe(int nodeid)
{
#ifndef USE_MPI
  MYFATAL("gmx_rx_probe");
  return 0;
#else
  MPI_Status  status;
  int mpi_result,flag=0;
  
  if ((mpi_result = MPI_Test(&mpi_req_rx,&flag,&status)) != MPI_SUCCESS)
    fatal_error(0,"MPI_Test: result=%d",mpi_result);
    
  return flag;
#endif
}

void gmx_rxs(int nodeid,void *buf,int bufsize)
{
#ifndef USE_MPI
  MYFATAL("gmx_rxs");
#else
  MPI_Status stat;
  int        tag;

#ifdef DEBUG
  fprintf(stderr,"gmx_rxs: nodeid=%d, buf=%x, bufsize=%d\n",
	  nodeid,buf,bufsize);
#endif
  tag = 0;
  if (MPI_Recv( buf, bufsize, MPI_BYTE, nodeid, tag, MPI_COMM_WORLD, &stat) != 0 )
    fatal_error(0,"MPI_Recv Failed !");
#endif
}

int gmx_setup(int *argc,char **argv,int *nnodes)
{
#ifndef USE_MPI
  MYFATAL("gmx_setup");
  return 0;
#else
  char   buf[256];
  int    resultlen;               /* actual length of node name      */
  int    i,flag;
  int  mpi_num_nodes;
  int  mpi_my_rank;
  char mpi_hostname[MPI_MAX_PROCESSOR_NAME];

  /* Call the MPI routines */
  (void) MPI_Init(argc,&argv);
  (void) MPI_Comm_size( MPI_COMM_WORLD, &mpi_num_nodes );
  (void) MPI_Comm_rank( MPI_COMM_WORLD, &mpi_my_rank );
  (void) MPI_Get_processor_name( mpi_hostname, &resultlen );

  fprintf(stderr,"NNODES=%d, MYRANK=%d, HOSTNAME=%s\n",
	  mpi_num_nodes,mpi_my_rank,mpi_hostname);
  
  *nnodes=mpi_num_nodes;
  
  return mpi_my_rank;
#endif
}

int  gmx_node_num(void)
{
#ifndef USE_MPI
  return 1;
#else
  int i;
  return MPI_Comm_size(MPI_COMM_WORLD, &i);
  return i;
#endif
}

int gmx_node_id(void)
{
#ifndef USE_MPI
  return 0;
#else
  int i;
  return MPI_Comm_rank(MPI_COMM_WORLD, &i);
  return i;
#endif
}

int gmx_idle_send(void)
{
  return 0;
}

int gmx_idle_rec(void)
{
  return 0;
}

void gmx_left_right(int nnodes,int nodeid,int *left,int *right)
{
  *left  = (nnodes+nodeid-1) % nnodes;
  *right = (nodeid+1) % nnodes;
}

void gmx_tx_rx(int send_nodeid,void *send_buf,int send_bufsize,
		 int rec_nodeid,void *rec_buf,int rec_bufsize)
{
#ifndef USE_MPI
  MYFATAL("gmx_tx_rx");
#else
  int tx_tag = 0,rx_tag = 0;
  MPI_Status stat;
  
  MPI_Sendrecv(send_buf,send_bufsize,MPI_BYTE,send_nodeid,tx_tag,
	       rec_buf,rec_bufsize,MPI_BYTE,rec_nodeid,rx_tag,
	       MPI_COMM_WORLD,&stat);
#endif
}
		 
void gmx_tx_rx_real(int send_nodeid,real *send_buf,int send_bufsize,
		      int rec_nodeid,real *rec_buf,int rec_bufsize)
{
#ifndef USE_MPI
  MYFATAL("gmx_tx_rx_real");
#else
  int tx_tag = 0,rx_tag = 0;
  MPI_Status stat;
#ifdef DOUBLE
#define mpi_type MPI_DOUBLE
#else
#define mpi_type MPI_FLOAT
#endif
  
  MPI_Sendrecv(send_buf,send_bufsize,mpi_type,send_nodeid,tx_tag,
	       rec_buf,rec_bufsize,mpi_type,rec_nodeid,rx_tag,
	       MPI_COMM_WORLD,&stat);
#undef mpi_type
#endif
}
		 
void gmx_wait(int left,int right)
{
#ifndef USE_MPI
  MYFATAL("gmx_wait");
#else
  gmx_tx_wait(left);
  gmx_rx_wait(right);
#endif
}

void gmx_sync_ring(int nodeid,int nnodes,int left,int right)
{
#ifndef USE_MPI
  MYFATAL("gmx_sync_ring");
#else
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
#endif
}

void gmx_stat(FILE *fp,char *msg)
{
  fprintf(fp,"def_stat: %s (from %s, %d)\n",msg,__FILE__,__LINE__);
}

void gmx_reset_idle(void)
{
  ;
}

void gmx_abort(int nodeid,int nnodes,int errorno)
{
#ifndef USE_MPI
  MYFATAL("gmx_abort");
#else
  MPI_Abort(MPI_COMM_WORLD,errorno);
#endif
}

void gmx_sumd(int nr,double r[],const t_commrec *cr)
{
#ifndef USE_MPI
  MYFATAL("gmx_sumd");
#else
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
#endif
}

void gmx_sumf(int nr,float r[],const t_commrec *cr)
{
#ifndef USE_MPI
  MYFATAL("gmx_sumf");
#else
  float *buf[2];
  int  NR,bufs,j,i,cur=0;
#define next (1-cur)

  bufs=nr*sizeof(float);
  NR=nr;

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
#endif
}

void gmx_sumi(int nr,int r[],const t_commrec *cr)
{
#ifndef USE_MPI
  MYFATAL("gmx_sumi");
#else
  int *buf[2];
  int  NR,bufs,j,i,cur=0;
#define next (1-cur)

  bufs=nr*sizeof(int);
  NR=nr;

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
#endif
}

void gmx_finalize(t_commrec *cr)
{
#ifndef USE_MPI
  MYFATAL("gmx_finalize");
#else
#ifdef MPICH_NAME
#ifdef DEBUG
  fprintf(stdlog,"In gmx_finalize. Will try to synchronize the ring\n");
#endif
  gmx_sync_ring(cr->nodeid,cr->nnodes,cr->left,cr->right);
#ifdef DEBUG
  fprintf(stdlog,"Succesfully did so! Exiting now.\n");
#endif
  thanx(stdlog);
  exit(0);
#else
#ifdef DEBUG
  fprintf(stdlog,"Will call MPI_Finalize now\n");
  fflush(stdlog);
#endif
  ret = MPI_Finalize();
#ifdef DEBUG
  fprintf(stdlog,"Return code from MPI_Finalize = %d\n",ret);
  fflush(stdlog);
#endif
#endif
#endif
}
