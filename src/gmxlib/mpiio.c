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
 * Great Red Owns Many ACres of Sand 
 */
static char *SRCID_mpiio_c = "$Id$";

#include <string.h>
#include "mpiio.h"
#include "fatal.h"
#include "main.h"
#include "smalloc.h"

static int  mpi_num_procs=0;
static int  mpi_my_rank=-1;
static char mpi_hostname[MPI_MAX_PROCESSOR_NAME];

static MPI_Request mpi_req_tx=MPI_REQUEST_NULL,mpi_req_rx;

/*#define DEBUG*/

#ifdef _SGI_
#define MPI_TEST
#endif

void mpiio_tx(int pid,void *buf,int bufsize)
{
  int        tag,flag;
  MPI_Status status;
  
#ifdef DEBUG
  fprintf(stderr,"mpiio_tx: pid=%d, buf=%x, bufsize=%d\n",pid,buf,bufsize);
#endif
#ifdef MPI_TEST
  /* workaround for crashes encountered with MPI on IRIX 6.5 */
  if (mpi_req_tx != MPI_REQUEST_NULL) {
    MPI_Test(&mpi_req_tx,&flag,&status);
    if (flag==FALSE) {
      fprintf(stdlog,"mpiio_tx called before previous send was complete: pid=%d, buf=%x, bufsize=%d\n",
	      pid,buf,bufsize);
      mpiio_tx_wait(pid);
    }
  }
#endif
  tag = 0;
  if (MPI_Isend(buf,bufsize,MPI_BYTE,pid,tag,MPI_COMM_WORLD,&mpi_req_tx) != 0)
    fatal_error(0,"MPI_Isend Failed !");
}

void mpiio_tx_wait(int pid)
{
  MPI_Status  status;
  int mpi_result;
  
  if ((mpi_result=MPI_Wait(&mpi_req_tx,&status)) != 0)
    fatal_error(0,"MPI_Wait: result=%d",mpi_result);
}

void mpiio_txs(int pid,void *buf,int bufsize)
{
  int tag;

#ifdef DEBUG
  fprintf(stderr,"mpiio_txs: pid=%d, buf=%x, bufsize=%d\n",pid,buf,bufsize);
#endif
  tag = 0;
  if (MPI_Send(buf,bufsize,MPI_BYTE,pid,tag,MPI_COMM_WORLD) != 0)
    fatal_error(0,"MPI_Send Failed !");
  /* mpiio_tx(pid,buf,bufsize);
     mpiio_tx_wait(pid);*/
}

void mpiio_rx(int pid,void *buf,int bufsize)
{
  int        tag;

#ifdef DEBUG
  fprintf(stderr,"mpiio_rx: pid=%d, buf=%x, bufsize=%d\n",pid,buf,bufsize);
#endif
  tag = 0;
  if (MPI_Irecv( buf, bufsize, MPI_BYTE, pid, tag, MPI_COMM_WORLD, &mpi_req_rx) != 0 )
    fatal_error(0,"MPI_Recv Failed !");
}

void mpiio_rx_wait(int pid)
{
  MPI_Status  status;
  int mpi_result;
  
  if ((mpi_result=MPI_Wait(&mpi_req_rx,&status)) != 0)
    fatal_error(0,"MPI_Wait: result=%d",mpi_result);
}

void mpiio_rxs(int pid,void *buf,int bufsize)
{
  MPI_Status stat;
  int        tag;

#ifdef DEBUG
  fprintf(stderr,"mpiio_rxs: pid=%d, buf=%x, bufsize=%d\n",pid,buf,bufsize);
#endif
  tag = 0;
  if (MPI_Recv( buf, bufsize, MPI_BYTE, pid, tag, MPI_COMM_WORLD, &stat) != 0 )
    fatal_error(0,"MPI_Recv Failed !");
  /* mpiio_rx(pid,buf,bufsize);
     mpiio_rx_wait(pid);*/
}

int mpiio_setup(int *argc,char **argv,int *nprocs)
{
  char   buf[256];
  int    resultlen;               /* actual length of processor name      */
  int    i,flag;

  /* Call the MPI routines */
  (void) MPI_Init(argc,&argv);
  (void) MPI_Comm_size( MPI_COMM_WORLD, &mpi_num_procs );
  (void) MPI_Comm_rank( MPI_COMM_WORLD, &mpi_my_rank );
  (void) MPI_Get_processor_name( mpi_hostname, &resultlen );

  fprintf(stderr,"NPROCS=%d, MYRANK=%d, HOSTNAME=%s\n",
	  mpi_num_procs,mpi_my_rank,mpi_hostname);
  
  *nprocs=mpi_num_procs;
  
  return mpi_my_rank;
}

/* If this is the master, spawn some kids. Return pid */

void mpiio_stat(FILE *fp,char *msg)
{
  ;
}

int  mpinodecount(void)
{
  return mpi_num_procs;
}

int mpinodenumber(void)
{
  return mpi_my_rank;
}

int mpi_idle_send(void)
{
  return 0;
}

int mpi_idle_rec(void)
{
  return 0;
}

void mpi_left_right(int nprocs,int pid,int *left,int *right)
{
  *left  = (nprocs+pid-1) % nprocs;
  *right = (pid+1) % nprocs;
}

void mpiio_tx_rx(int send_pid,void *send_buf,int send_bufsize,
		 int rec_pid,void *rec_buf,int rec_bufsize)
{
  fatal_error(0,"mpiio_tx_rx not implemented!");
}
		 
void mpiio_wait(int left,int right)
{
  mpiio_tx_wait(left);
  mpiio_rx_wait(right);
}

void mpiio_sync_ring(int pid,int nprocs,int left,int right)
{
  fatal_error(0,"mpiio_sync_ring called");;
}

void mpi_reset_idle(void)
{
  ;
}

void mpi_abort(int pid,int nprocs,int errorno)
{
  MPI_Abort(MPI_COMM_WORLD,errorno);
}

#ifdef DEBUG
#undef DEBUG
#endif
