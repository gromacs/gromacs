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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <time.h>
#include "smalloc.h"
#include "gmx_fatal.h"
#include "network.h"
#include "main.h"
#include "macros.h"
#include "futil.h"
#include "filenm.h"

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#define BUFSIZE	1024






FILE *stdlog=NULL;
int  gmx_parallel_env=0;

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

static void par_fn(char *base,int ftp,const t_commrec *cr,
		   char buf[],int bufsize)
{
  int n;
  
  if(bufsize<(strlen(base)+4))
     gmx_mem("Character buffer too small!");

  /* Copy to buf, and strip extension */
  strcpy(buf,base);
  buf[strlen(base) - strlen(ftp2ext(fn2ftp(base))) - 1] = '\0';

  /* Add node info */
  if (PAR(cr)) 
    sprintf(buf+strlen(buf),"%d",cr->nodeid);
  strcat(buf,".");
  
  /* Add extension again */
  strcat(buf,(ftp == efTPX) ? "tpr" : (ftp == efENX) ? "edr" : ftp2ext(ftp));
}

void check_multi_int(FILE *log,const t_commrec *mcr,int val,char *name)
{
  int  *ibuf,p;
  bool bCompatible;

  fprintf(log,"Multi-checking %s ... ",name);
  
  snew(ibuf,mcr->nnodes);
  ibuf[mcr->nodeid] = val;
  if((mcr) && PAR(mcr))
    gmx_sumi(mcr->nnodes,ibuf,mcr);
  
  bCompatible = TRUE;
  for(p=1; p<mcr->nnodes; p++)
    bCompatible = bCompatible && (ibuf[p-1] == ibuf[p]);
  
  if (bCompatible) 
    fprintf(log,"OK\n");
  else {
    fprintf(log,"\n%s is not equal for all subsystems\n",name);
    for(p=0; p<mcr->nnodes; p++)
      fprintf(log,"  subsystem %d: %d\n",p,ibuf[p]);
    gmx_fatal(FARGS,"The %d subsystems are not compatible\n",mcr->nnodes);
  }
  
  sfree(ibuf);
}

void open_log(char *lognm,const t_commrec *cr)
{
  int  len,testlen,pid;
  char buf[256],host[256];
  time_t t;

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
	gmx_comm("Communication error on NODE 0!");
      
    }
    else {
      gmx_rxs(cr->left,&len,sizeof(len));
      debug_gmx();
      
      gmx_txs(cr->right,&len,sizeof(len));
      snew(lognm,len+8);
      gmx_rxs(cr->left,lognm,len);
      gmx_txs(cr->right,lognm,len);
    }
  }
  
  debug_gmx();

  /* Since log always ends with '.log' let's use this info */
  par_fn(lognm,efLOG,cr,buf,255);
  stdlog = ffopen(buf,"w");
  
  /* Get some machine parameters */
#ifdef HAVE_UNISTD_H
  if( gethostname(host,255) != 0)
    sprintf(host,"unknown");
#else
  sprintf(host,"unknown");
#endif  

  time(&t);

#ifndef NO_GETPID
  pid = getpid();
#else
  pid = 0;
#endif

  fprintf(stdlog,
	  "Log file opened on %s"
	  "Host: %s  pid: %d  nodeid: %d  nnodes:  %d\n",
	  ctime(&t),host,pid,cr->nodeid,cr->nnodes);

#if (defined BUILD_MACHINE && defined BUILD_TIME && defined BUILD_USER) 
  fprintf(stdlog,
	  "The Gromacs distribution was built %s by\n"
	  "%s (%s)\n\n\n",BUILD_TIME,BUILD_USER,BUILD_MACHINE);
#endif

  fflush(stdlog);
  debug_gmx();
}

static void comm_args(const t_commrec *cr,int *argc,char ***argv)
{
  int i,len;
  char **argv_tmp=NULL,*buf;
  
  if (!MASTER(cr))
    *argc=0;

  if((cr) && PAR(cr))
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
	gmx_fatal(FARGS,"Communicating argv[%d]=%s\n",i,(*argv)[i]);
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

t_commrec *init_multisystem(t_commrec *cr,int nfile,t_filenm fnm[],bool bParFn)
{
  t_commrec *mcr;
  int  i,ftp;
  char buf[256];
  
  snew(mcr,1);

  mcr->nodeid = cr->nodeid;
  mcr->nnodes = cr->nnodes;
  mcr->left   = cr->left;
  mcr->right  = cr->right;
  cr->nodeid  = 0;
  cr->nnodes  = 1;
  
  if (bParFn) {
    /* Patch output and tpx file names (except log which has been done already)
     */
    for(i=0; (i<nfile); i++) {
      /* Because of possible multiple extensions per type we must look 
       * at the actual file name 
       */
      if ((is_output(&fnm[i]) || fnm[i].ftp == efTPX) && fnm[i].ftp != efLOG) {
	ftp = fn2ftp(fnm[i].fns[0]);
	par_fn(fnm[i].fns[0],ftp,mcr,buf,255);
	sfree(fnm[i].fns[0]);
	fnm[i].fns[0] = strdup(buf);
      }
    }
  }

  return mcr;
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
  
#ifdef GMX_MPI
  gmx_parallel_env = 1;
#ifdef GMX_CHECK_MPI_ENV
  /* Do not use MPI calls when env.var. GMX_CHECK_MPI_ENV is not set */
  if (getenv(GMX_CHECK_MPI_ENV) == NULL)
    gmx_parallel_env = 0;
#endif
  if (gmx_parallel_env)
    cr->nodeid = gmx_setup(argc,argv,&cr->nnodes);
  else
    cr->nodeid = 0;
#else
  cr->nodeid   = 0;
  cr->nnodes   = 1;
  gmx_parallel_env = 0; 
#endif
  
  if (!PAR(cr) && (cr->nodeid != 0))
    gmx_comm("(!PAR(cr) && (cr->nodeid != 0))");
  
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
      gmx_comm("could not get nodeid & nnodes from ring topology");
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

