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
 * Glycine aRginine prOline Methionine Alanine Cystine Serine
 */
static char *SRCID_main_c = "$Id$";
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include "smalloc.h"
#include "fatal.h"
#include "network.h"
#include "main.h"
#include "macros.h"
#include "futil.h"
#include "filenm.h"

#define BUFSIZE	1024

FILE *stdlog=NULL;
int  gmx_parallel=0;

static void mem_init(void)
{
#ifdef _860_
  void *p;
  unsigned long avail;

  avail=maxavail();
  if ((p=malloc(avail))==NULL)
    fatal_error(errno,"mem init");
  else
    {
      memset(p,0,avail);
      free(p);
    }
#endif
}

static int get_nodeid(FILE *log,int left,int right,int *nodeid,int *nnodes)
     /*
      * The ring of nodes is only defined by the interconnection
      * via the supplied communication channels (left and right). Thus
      * it is not defined what the (hardware) node id's are in the
      * ring. To be independent of the node id assignment (to allow
      * node switching without modifying a node id) this routine
      * determines the node id and the number of nodes. On entry
      * nodeid needs to be set to an unique value in the system and nnodes
      * needs to be set to the maximum number of nodes in the system.
      * The lowest nodeid in the ring will then get the value 0 assigned. The
      * rest is then assigned by incrementing node id's to the right
      * until the ring is closed. The function returns 1 in case it succeeded
      * in determining the values for nodeid and nnodes, else it returns 0. If
      * the hardware does not implement a ring structure this will hang the
      * system!
      */
{
  int i,nodeids[MAXNODES],min_index,min_nodeid,send_nodeid,receive_nodeid;

  *nnodes=0;
  send_nodeid=*nodeid;
  min_nodeid=send_nodeid;
  min_index=*nnodes;
  do {
#ifdef DEBUGPAR
    fprintf(log,"Sending: %d\n",send_nodeid);
#endif
    gmx_tx(left,&send_nodeid,sizeof(send_nodeid));
    gmx_rx(right,&receive_nodeid,sizeof(receive_nodeid));
    gmx_tx_wait(left);
    gmx_rx_wait(right);
#ifdef DEBUGPAR
    fprintf(log,"Received: %d\n",receive_nodeid);
#endif
    if (send_nodeid<min_nodeid) {
      min_nodeid=send_nodeid;
      min_index=*nnodes;
    }
    nodeids[(*nnodes)++]=send_nodeid;
    send_nodeid=receive_nodeid;
  } while (receive_nodeid!=*nodeid);
  
#ifdef DEBUGPAR  
  fprintf(log,"min_index=%d\n",min_index);
  fprintf(log,"nnodes   =%d\n",*nnodes);
  fprintf(log,"nodeid      =%d\n",*nodeid);
#endif

  for (i=min_index; (*nodeid)!=nodeids[i%(*nnodes)]; i++)
    ;
  (*nodeid)=(i-min_index+(*nnodes))%(*nnodes);
#ifdef DEBUGPAR
  fprintf(log,"min_index=%d\n",min_index);
  fprintf(log,"nnodes   =%d\n",*nnodes);
  fprintf(log,"nodeid      =%d\n",*nodeid);
  for (i=0; i<(*nnodes); i++) {
    fprintf(log,"%d translated %d --> %d",
	    i,nodeids[i],(i-min_index+(*nnodes))%(*nnodes));
    if (nodeids[i]==(*nodeid)) 
      fprintf(log," *");
    fprintf(log,"\n");
  }
#endif
  return 1;
}

char *par_fn(char *base,int ftp,t_commrec *cr)
{
  static char buf[256];
  
  /* Copy to buf, and strip extension */
  strcpy(buf,base);
  buf[strlen(base)-4] = '\0';
  
  /* Add node info */
  if (PAR(cr)) 
    sprintf(buf+strlen(buf),"%d",cr->nodeid);
  strcat(buf,".");
  
  /* Add extension again */
  strcat(buf,(ftp == efTPX) ? "tpr" : (ftp == efENX) ? "edr" : ftp2ext(ftp));
  
  return buf;
}

void open_log(char *lognm,t_commrec *cr)
{
  int  len,testlen,pid;
  char *buf,*host;
  
  debug_gmx();
  
  /* Communicate the filename for logfile */
  if (cr->nnodes > 1) {
    if (MASTER(cr)) {
      len = strlen(lognm)+1;
      gmx_txs(cr->right,&len,sizeof(len));
      gmx_rxs(cr->left,&testlen,sizeof(testlen));
      
      debug_gmx();
      
      gmx_txs(cr->right,lognm,len);
      gmx_rxs(cr->left,lognm,len);
      if (len != testlen)
	fatal_error(0,"Communication error on NODE 0!");
      
    }
    else {
      gmx_rxs(cr->left,&len,sizeof(len));
      debug_gmx();
      
      gmx_txs(cr->right,&len,sizeof(len));
      snew(lognm,len);
      gmx_rxs(cr->left,lognm,len);
      gmx_txs(cr->right,lognm,len);
    }
  }
  
  debug_gmx();

  /* Since log always ends with '.log' let's use this info */
  buf    = par_fn(lognm,efLOG,cr);
  stdlog = ffopen(buf,"w");
  
  /* Get some machine parameters */
  host = getenv("HOST");
#ifndef NO_GETPID
  pid = getpid();
#else
  pid = 0;
#endif
  fprintf(stdlog,"Log file opened: nodeid %d, nnodes = %d, host = %s, process = %d\n",
	  cr->nodeid,cr->nnodes,host ? host : "unknown",pid);
  fflush(stdlog);
  debug_gmx();
}

static void comm_args(t_commrec *cr,int *argc,char ***argv)
{
  int i,len;
  char **argv_tmp=NULL,*buf;
  
  if (!MASTER(cr))
    *argc=0;
  gmx_sumi(1,argc,cr);
  
  if (!MASTER(cr))
    snew(argv_tmp,*argc+1);
  fprintf(stderr,"NODEID=%d argc=%d\n",cr->nodeid,*argc);
  for(i=0; (i<*argc); i++) {
    if (MASTER(cr)) {
      len = strlen((*argv)[i])+1;
      gmx_txs(cr->right,&len,sizeof(len));
      gmx_rxs(cr->left,&len,sizeof(len));
      gmx_txs(cr->right,(*argv)[i],len);
      snew(buf,len);
      gmx_rxs(cr->left,buf,len);
      if (strcmp(buf,(*argv)[i]) != 0)
	fatal_error(0,"Communicating argv[%d]=%s\n",i,(*argv)[i]);
      sfree(buf);
    }
    else {
      gmx_rxs(cr->left,&len,sizeof(len));
      gmx_txs(cr->right,&len,sizeof(len));
      snew(argv_tmp[i],len);
      gmx_rxs(cr->left,argv_tmp[i],len);
      gmx_txs(cr->right,argv_tmp[i],len);
    }
  }
  if (!MASTER(cr)) {
    argv_tmp[*argc] = NULL;
    *argv = argv_tmp;
  }
  debug_gmx();
}

t_commrec *init_par(int *argc,char ***argv_ptr)
{
  t_commrec *cr;
  char      **argv;
  int       i;
  
  argv = *argv_ptr;
  snew(cr,1);
  
  cr->nnodes=1;
  /* Get the number of nodes.
   * This is useless for newer MPI versions.
   */
  for(i=0; (argv[i] != NULL); i++) {
    if (strcmp(argv[i],"-np")==0)
      if (argv[i+1]!=NULL)
	cr->nnodes=atoi(argv[i+1]);
  }
  
#ifdef USE_MPI
  gmx_parallel = 1;
#ifdef CHECK_MPI_ENV
  /* Do not use MPI calls when env.var. CHECK_MPI_ENV is not set */
  if (getenv(CHECK_MPI_ENV) == NULL)
    gmx_parallel = 0;
#endif
  if (gmx_parallel)
    cr->nodeid=gmx_setup(argc,argv,&cr->nnodes);
  else
    cr->nodeid=0;
#else
  cr->nodeid=0;
  if (cr->nnodes > 1) {
    gmx_parallel = 1; 
#ifdef USE_PVM3
    cr->nodeid=pvmio_setup(argv,cr->nnodes);
#else
    cr->nodeid=gmx_node_id();
#endif
  }
#endif
  
#ifdef _amb_
  pageMode(1);
#endif
  mem_init();

  if (!PAR(cr) && (cr->nodeid != 0))
    fatal_error(0,"(!PAR(cr) && (cr->nodeid != 0))");
  
  if (PAR(cr)) {
    gmx_left_right(cr->nnodes,cr->nodeid,&cr->left,&cr->right);
#ifdef DEBUGPAR
    fprintf(stderr,"Going to initialise network\n");
#endif

#ifndef USE_PVM3
#ifdef DEBUGPAR
    fprintf(stderr,"Initialised network\n");
    fprintf(stderr,"Getting new node id's\n");
#endif
    if (get_nodeid(stderr,cr->left,cr->right,&cr->nodeid,&cr->nnodes)==0)
      fatal_error(0,"could not get nodeid & nnodes from ring topology");
#ifdef DEBUGPAR
    fprintf(stderr,"Got new node id's\n");
    fprintf(stderr,"nnodes=%d, nodeid=%d\n",cr->nnodes,cr->nodeid);
#endif

#endif
  }

  /* Communicate arguments if parallel */
  if (PAR(cr))
    comm_args(cr,argc,argv_ptr);

  return cr;
}

