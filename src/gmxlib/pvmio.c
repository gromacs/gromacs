/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 1.6
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
 * Great Red Oystrich Makes All Chemists Sane
 */
static char *SRCID_pvmio_c = "$Id$";

#include <stdio.h>
#include "typedefs.h"
#include "network.h"
#include "macros.h"
#include "main.h"

#ifdef USE_PVM3

#include "pvm3.h"

#define DEFAULT_STRIDE  1
#define ANY_TID        -1

typedef struct { 
  int errnok; 
  char *msgk;
} t_pvmerr;

static int   tids[MAXPROC];
static int   mypid=0,mynprocs=1;
static ulong idle_send=0,idle_rec=0;

void pvmio_tx(int dest,void *buf,int bufsize)
{
  int bufid,info;
#ifdef DEBUG
  fprintf(stderr,"dest: %d, type: %d \n",dest,PVM_BYTE);
#endif
  bufid = pvm_initsend(PvmDataRaw);
  info  = pvm_pkbyte(buf,bufsize,DEFAULT_STRIDE);
  info  = pvm_send(dest,PVM_BYTE);
}

char *pvm_error(int errorno)
{
  static t_pvmerr pvmerrs[] = {
    -2,  "PvmBadParam (giving an invalid argument value)",
    -6,  "PvmNoHost (specified host is not on the virtual machine)",
    -7,  "PvmNoFile (Can not find executable)",
    -10, "PvmNoMem (not enough memory)",
    -14, "PvmSysErr (pvmd not responding)",
    -27, "PvmOutOfRes (out of resources)"
    };
  static char *deferr="Unknown Pvm Error";
  int i;
  
  for(i=0; (i<asize(pvmerrs)); i++)
    if (errorno == pvmerrs[i].errnok)
      return pvmerrs[i].msgk;
      
  return deferr;
}

void pvmio_tx_wait(int dest)
{
#ifdef PROFILING
  idle_send++;
#endif
}

void pvmio_txs(int dest,void *buf,int bufsize)
{
  pvmio_tx(dest,buf,bufsize);
  pvmio_tx_wait(dest);
}

void pvmio_rx(int src,void *buf,int bufsize)
{
  int bufid;
  
  bufid = pvm_recv(src,PVM_BYTE);
  if (bufid) 
    pvm_upkbyte(buf,bufsize,DEFAULT_STRIDE);
  else
    fatal_error(0,"Did not receive %d bytes from src %x",bufsize,src);
}

void pvmio_rx_wait(int src)
{
#ifdef PROFILING
  idle_rec++;
#endif
}

void pvmio_rxs(int src,void *buf,int bufsize)
{
  pvmio_rx(src,buf,bufsize);
  pvmio_rx_wait(src);
}


int pvmio_setup(char *argv[],int nprocs)
{
  int i,bufid,info;
  int pid,numt;

  if (pvm_parent() == PvmNoParent) {
    /* I am the MASTER! Obey me... */
    pid=0;
    tids[0]=pvm_mytid();
    fprintf(stderr,"PVM_MD: Spawning %d tasks\n",nprocs-1);
    numt=pvm_spawn(argv[0],&(argv[1]),PvmTaskDefault,
		   "",nprocs-1,tids+1);
    if (numt != nprocs-1) {
      fprintf(stderr,"Could only spawn %d tasks!\n",numt);
      for(i=0; (i<nprocs); i++) {
	if (tids[i] < 0)
	  fprintf(stderr,"PVM_MD: tids[%d] = %d --> %s\n",
		  i,tids[i],pvm_error(tids[i]));
	else
	  fprintf(stderr,"PVM_MD: tids[%d] = %d\n",i,tids[i]);
      }
      fatal_error(0,"Could not spawn %d processes, numt=%d",nprocs-1,numt);
    }
    for(i=0; (i<nprocs); i++)
      fprintf(stderr,"tids[%d]=%x\n",i,tids[i]);
    
    for (i=1; (i<nprocs); i++) {
      fprintf(stderr,"pid %d: sending tids to pid %d = tid %x \n",
	      pid,i,tids[i]); 
      /* Send target id's to everyone */
      pvmio_tx(tids[i],tids,nprocs*sizeof(tids[0]));
      
      /* Send index to other CPUS */
      pvmio_tx(tids[i],&i,sizeof(i));
    } 
  }
  else {
    /* I am a slave, get some tids */
    bufid = pvm_recv(ANY_TID,PVM_BYTE);
    info  = pvm_upkbyte((char *)tids,nprocs*sizeof(tids[0]),DEFAULT_STRIDE);
    fprintf(stderr,"Got tids...\n");
    
    /* Get my own number */
    pvmio_rx(tids[0],&pid,sizeof(pid));
  }
  mypid    = pid;
  mynprocs = nprocs;
  
  fprintf(stderr,"my pid is: %d\n",pid);
  
  return pid;
}

void pvmio_stat(FILE *fp,char *msg)
{
  fprintf(fp,"pvmio_stat message: %s\n",msg);
  fprintf(fp,"Idle Send:    %d\n",idle_send);
  fprintf(fp,"Idle Receive: %d\n",idle_rec);
}

int pvmnodenumber(void)
{ 
  return mypid;
}

int pvmnodecount(void)
{ 
  return mynprocs;
}

int pvm_idle_send(void)
{
  return idle_send;
}

int pvm_idle_rec(void) 
{
  return idle_rec;
}

void pvm_left_right(int nprocs,int pid,int *left,int *right)
{
  *left=tids[(pid-1+nprocs) % nprocs];
  *right=tids[(pid+1) % nprocs];
  fprintf(stderr,"PVM: left=%d, right=%d\n",*left,*right); fflush(stderr);
}

void pvm_reset_idle()
{
  idle_send=0,idle_rec=0;
}

void pvmio_tx_rx(int send_pid,void *send_buf,int send_bufsize,
		 int rec_pid,void *rec_buf,int rec_bufsize)
{
  fatal_error(0,"pvmio_tx_rx called");
}

void pvmio_wait(int send,int receive)
{
  pvmio_tx_wait(send);
  pvmio_rx_wait(receive);
}

void pvmio_sync_ring(int pid,int nprocs,int left,int right)
{
  fatal_error(0,"pvmio_sync_ring called");
}

#endif
/* USE_PVM3 */
