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

static int   tids[MAXNODES];
static int   mynodeid=0,mynnodes=1;
static unsigned long idle_send=0,idle_rec=0;

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


int pvmio_setup(char *argv[],int nnodes)
{
  int i,bufid,info;
  int nodeid,numt;

  if (pvm_parent() == PvmNoParent) {
    /* I am the MASTER! Obey me... */
    nodeid=0;
    tids[0]=pvm_mytid();
    fprintf(stderr,"PVM_MD: Spawning %d tasks\n",nnodes-1);
    numt=pvm_spawn(argv[0],&(argv[1]),PvmTaskDefault,
		   "",nnodes-1,tids+1);
    if (numt != nnodes-1) {
      fprintf(stderr,"Could only spawn %d tasks!\n",numt);
      for(i=0; (i<nnodes); i++) {
	if (tids[i] < 0)
	  fprintf(stderr,"PVM_MD: tids[%d] = %d --> %s\n",
		  i,tids[i],pvm_error(tids[i]));
	else
	  fprintf(stderr,"PVM_MD: tids[%d] = %d\n",i,tids[i]);
      }
      fatal_error(0,"Could not spawn %d nodes, numt=%d",nnodes-1,numt);
    }
    for(i=0; (i<nnodes); i++)
      fprintf(stderr,"tids[%d]=%x\n",i,tids[i]);
    
    for (i=1; (i<nnodes); i++) {
      fprintf(stderr,"nodeid %d: sending tids to nodeid %d = tid %x \n",
	      nodeid,i,tids[i]); 
      /* Send target id's to everyone */
      pvmio_tx(tids[i],tids,nnodes*sizeof(tids[0]));
      
      /* Send index to other CPUS */
      pvmio_tx(tids[i],&i,sizeof(i));
    } 
  }
  else {
    /* I am a slave, get some tids */
    bufid = pvm_recv(ANY_TID,PVM_BYTE);
    info  = pvm_upkbyte((char *)tids,nnodes*sizeof(tids[0]),DEFAULT_STRIDE);
    fprintf(stderr,"Got tids...\n");
    
    /* Get my own number */
    pvmio_rx(tids[0],&nodeid,sizeof(nodeid));
  }
  mynodeid    = nodeid;
  mynnodes = nnodes;
  
  fprintf(stderr,"my nodeid is: %d\n",nodeid);
  
  return nodeid;
}

void pvmio_stat(FILE *fp,char *msg)
{
  fprintf(fp,"pvmio_stat message: %s\n",msg);
  fprintf(fp,"Idle Send:    %d\n",idle_send);
  fprintf(fp,"Idle Receive: %d\n",idle_rec);
}

int pvmnodenumber(void)
{ 
  return mynodeid;
}

int pvmnodecount(void)
{ 
  return mynnodes;
}

int pvm_idle_send(void)
{
  return idle_send;
}

int pvm_idle_rec(void) 
{
  return idle_rec;
}

void pvm_left_right(int nnodes,int nodeid,int *left,int *right)
{
  *left=tids[(nodeid-1+nnodes) % nnodes];
  *right=tids[(nodeid+1) % nnodes];
  fprintf(stderr,"PVM: left=%d, right=%d\n",*left,*right); fflush(stderr);
}

void pvm_reset_idle()
{
  idle_send=0,idle_rec=0;
}

void pvmio_tx_rx(int send_nodeid,void *send_buf,int send_bufsize,
		 int rec_nodeid,void *rec_buf,int rec_bufsize)
{
  fatal_error(0,"pvmio_tx_rx called");
}

void pvmio_wait(int send,int receive)
{
  pvmio_tx_wait(send);
  pvmio_rx_wait(receive);
}

void pvmio_sync_ring(int nodeid,int nnodes,int left,int right)
{
  fatal_error(0,"pvmio_sync_ring called");
}

void pvm_abort(int nodeid,int nnodes,int errno)
{
  int i;

  for(i=0; i<nnodes; i++)
    if (i != nodeid)
      pvm_sendsig(tids[i],SIGKILL);
  
  exit(errno);
}

#endif
/* USE_PVM3 */
