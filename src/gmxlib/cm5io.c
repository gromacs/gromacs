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
 * S  C  A  M  O  R  G
 */
static char *SRCID_cm5io_c = "$Id$";

#include <stdio.h>
#include <cm/cmmd.h>
/*#include <cm/cmmd-io.h>*/
#include <network.h>

#define TAG   CMMD_DEFAULT_TAG

static CMMD_mcb mcb_send;
static CMMD_mcb mcb_rec;
static unsigned long idle_send=0,idle_rec=0;

void cm5io_tx(int pid,void *buf,int bufsize)
{
  mcb_send=CMMD_send_async(pid,TAG,buf,bufsize,NULL,NULL);
}

void cm5io_tx_wait(int pid)
{
  while (!CMMD_msg_done(mcb_send))
#ifdef PROFILING
    idle_send++
#endif
      ;
  CMMD_free_mcb(mcb_send);
}

void cm5io_txs(int pid,void *buf,int bufsize)
{
  cm5io_tx(pid,buf,bufsize);
  cm5io_tx_wait(pid);
}

void cm5io_rx(int pid,void *buf,int bufsize)
{
  mcb_rec=CMMD_receive_async(pid,TAG,buf,bufsize,NULL,NULL);
}

void cm5io_rx_wait(int pid)
{
  while (!CMMD_msg_done(mcb_rec))
#ifdef PROFILING
    idle_rec++
#endif
      ;
  CMMD_free_mcb(mcb_rec);
}

void cm5io_wait(int send,int receive)
{
  int bSend=0,bRec=0;

  do {
    bSend = bSend || CMMD_msg_done(mcb_send);
    bRec  = bRec  || CMMD_msg_done(mcb_rec);
    if (!bSend) idle_send++;
    if (!bRec)  idle_rec++;
  } while (!(bSend && bRec));
  CMMD_free_mcb(mcb_rec);
  CMMD_free_mcb(mcb_send);
}

void cm5io_rxs(int pid,void *buf,int bufsize)
{
  cm5io_rx(pid,buf,bufsize);
  cm5io_rx_wait(pid);
}

void cm5io_tx_rx(int send_pid,void *send_buf,int send_bufsize,
		 int rec_pid,void *rec_buf,int rec_bufsize)
{
  int ret;

  ret=CMMD_send_and_receive(rec_pid,TAG,rec_buf,rec_bufsize,
			    send_pid,TAG,send_buf,send_bufsize);
  if (ret != TRUE) {
    fprintf(stderr,"Only %d bytes sent and received...\n",ret);
    exit(1);
  }
}

void cm5io_init(int pid,int nprocs)
{
  CMMD_fset_io_mode(stdin,CMMD_independent);
  CMMD_fset_io_mode(stdout,CMMD_independent);
  CMMD_fset_io_mode(stderr,CMMD_independent);
  CMMD_sync_with_nodes();
  /* CMMD_node_timer_clear(0); */
}

void cm5io_stat(FILE *fp,char *msg)
{
  fprintf(fp,"cm5io_stat message: %s\n",msg);
  fprintf(fp,"Idle Send:    %d\n",idle_send);
  fprintf(fp,"Idle Receive: %d\n",idle_rec);
}

int cm5_cpu_id()
{ 
  return(CMMD_self_address());
}

int cm5_cpu_num()
{
  return (CMMD_partition_size());
}

int cm5_idle_send()
{
  return idle_send;
}

int cm5_idle_rec()
{
  return idle_rec;
}

void cm5_left_right(int nprocs,int pid,int *left,int *right)
{
  *left=(pid-1+nprocs) % nprocs;
  *right=(pid+1) % nprocs;
}

void reset_idle()
{
  idle_send=0,idle_rec=0;
}

void cm5_sync_ring(int pid,int nprocs,int left,int right)
{
  CMMD_sync_with_nodes();
}
