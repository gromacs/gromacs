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
static char *SRCID_cm5io_c = "$Id$";

#include <stdio.h>
#include <cm/cmmd.h>
/*#include <cm/cmmd-io.h>*/
#include <network.h>

#define TAG   CMMD_DEFAULT_TAG

static CMMD_mcb mcb_send;
static CMMD_mcb mcb_rec;
static unsigned long idle_send=0,idle_rec=0;

void cm5io_tx(int nodeid,void *buf,int bufsize)
{
  mcb_send=CMMD_send_async(nodeid,TAG,buf,bufsize,NULL,NULL);
}

void cm5io_tx_wait(int nodeid)
{
  while (!CMMD_msg_done(mcb_send))
#ifdef PROFILING
    idle_send++
#endif
      ;
  CMMD_free_mcb(mcb_send);
}

void cm5io_txs(int nodeid,void *buf,int bufsize)
{
  cm5io_tx(nodeid,buf,bufsize);
  cm5io_tx_wait(nodeid);
}

void cm5io_rx(int nodeid,void *buf,int bufsize)
{
  mcb_rec=CMMD_receive_async(nodeid,TAG,buf,bufsize,NULL,NULL);
}

void cm5io_rx_wait(int nodeid)
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

void cm5io_rxs(int nodeid,void *buf,int bufsize)
{
  cm5io_rx(nodeid,buf,bufsize);
  cm5io_rx_wait(nodeid);
}

void cm5io_tx_rx(int send_nodeid,void *send_buf,int send_bufsize,
		 int rec_nodeid,void *rec_buf,int rec_bufsize)
{
  int ret;

  ret=CMMD_send_and_receive(rec_nodeid,TAG,rec_buf,rec_bufsize,
			    send_nodeid,TAG,send_buf,send_bufsize);
  if (ret != TRUE)
    fatal_error(0,"Only %d bytes sent and received...",ret);
}

void cm5io_init(int nodeid,int nnodes)
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

int cm5_node_id()
{ 
  return(CMMD_self_address());
}

int cm5_node_num()
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

void cm5_left_right(int nnodes,int nodeid,int *left,int *right)
{
  *left=(nodeid-1+nnodes) % nnodes;
  *right=(nodeid+1) % nnodes;
}

void reset_idle()
{
  idle_send=0,idle_rec=0;
}

void cm5_sync_ring(int nodeid,int nnodes,int left,int right)
{
  CMMD_sync_with_nodes();
}
