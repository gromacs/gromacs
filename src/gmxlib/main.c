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
static char *SRCID_main_c = "$Id$";

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include "smalloc.h"
#include "fatal.h"
#include "led.h"
#include "network.h"
#include "main.h"
#include "macros.h"
#include "futil.h"
#include "filenm.h"

/*#define DEBUGPAR */
#define BUFSIZE	1024

FILE *stdlog=NULL;

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

static int get_pid(FILE *log,int left,int right,int *pid,int *nprocs)
     /*
      * The ring of processors is only defined by the interconnection
      * via the supplied communication channels (left and right). Thus
      * it is not defined what the (hardware) processor id's are in the
      * ring. To be independent of the processor id assignment (to allow
      * procesor switching without modifying a processor id) this routine
      * determines the processor id and the number of processors. On entry
      * pid needs to be set to an unique value in the system and nprocs
      * needs to be set to the maximum number of processors in the system.
      * The lowest pid in the ring will then get the value 0 assigned. The
      * rest is then assigned by incrementing processor id's to the right
      * until the ring is closed. The function returns 1 in case it succeeded
      * in determining the values for pid and nprocs, else it returns 0. If
      * the hardware does not implement a ring structure this will hang the
      * system!
      */
{
  int i,pids[MAXPROC],min_index,min_pid,send_pid,receive_pid;

  *nprocs=0;
  send_pid=*pid;
  min_pid=send_pid;
  min_index=*nprocs;
  do {
#ifdef DEBUGPAR
    fprintf(log,"Sending: %d\n",send_pid);
#endif
    gmx_tx(left,record(send_pid));
    gmx_rx(right,record(receive_pid));
    gmx_tx_wait(left);
    gmx_rx_wait(right);
#ifdef DEBUGPAR
    fprintf(log,"Received: %d\n",receive_pid);
#endif
    if (send_pid<min_pid) {
      min_pid=send_pid;
      min_index=*nprocs;
    }
    pids[(*nprocs)++]=send_pid;
    send_pid=receive_pid;
  } while (receive_pid!=*pid);
  
#ifdef DEBUGPAR  
  fprintf(log,"min_index=%d\n",min_index);
  fprintf(log,"nprocs   =%d\n",*nprocs);
  fprintf(log,"pid      =%d\n",*pid);
#endif

  for (i=min_index; (*pid)!=pids[i%(*nprocs)]; i++)
    ;
  (*pid)=(i-min_index+(*nprocs))%(*nprocs);
#ifdef DEBUGPAR
  fprintf(log,"min_index=%d\n",min_index);
  fprintf(log,"nprocs   =%d\n",*nprocs);
  fprintf(log,"pid      =%d\n",*pid);
  for (i=0; i<(*nprocs); i++) {
    fprintf(log,"%d translated %d --> %d",
	    i,pids[i],(i-min_index+(*nprocs))%(*nprocs));
    if (pids[i]==(*pid)) 
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
  
  /* Add processor info */
  if (PAR(cr)) 
    sprintf(buf+strlen(buf),"%d",cr->pid);
  strcat(buf,".");
  
  /* Add extension again */
  strcat(buf,(ftp == efTPX) ? "tpr" : (ftp == efENX) ? "edr" : ftp2ext(ftp));
  
  return buf;
}

void open_log(char *lognm,t_commrec *cr)
{
  int  len,testlen,pid;
  char *buf,*host;
  
#ifdef DEBUGPAR
  debug_par();
  fprintf(stderr,"OPEN_LOG: cr->pid = %d\n",cr->pid);
#endif
  /* Communicate the filename for logfile */
  if (cr->nprocs > 1) {
    if (MASTER(cr)) {
      len = strlen(lognm)+1;
      gmx_txs(cr->right,record(len));
      gmx_rxs(cr->left,record(testlen));
#ifdef DEBUGPAR
      fprintf(stderr,"cr->pid = %d, len = %d\n",cr->pid,testlen);
#endif
      gmx_txs(cr->right,lognm,len);
      gmx_rxs(cr->left,lognm,len);
      if (len != testlen)
	fatal_error(0,"Communication error on PROC 0!");
      
    }
    else {
      gmx_rxs(cr->left,record(len));
#ifdef DEBUGPAR
      fprintf(stderr,"cr->pid = %d, len = %d\n",cr->pid,len);
#endif
      gmx_txs(cr->right,record(len));
      snew(lognm,len);
      gmx_rxs(cr->left,lognm,len);
      gmx_txs(cr->right,lognm,len);
    }
  }
  
#ifdef DEBUGPAR
  debug_par();
#endif
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
  fprintf(stdlog,"Log file opened: pid %d, nprocs = %d, host = %s, process = %d\n",
	  cr->pid,cr->nprocs,host ? host : "unknown",pid);
  fflush(stdlog);
#ifdef DEBUGPAR
  debug_par();
#endif
}

static void comm_args(t_commrec *cr,int *argc,char ***argv)
{
  int i,len;
  char **argv_tmp,*buf;
  
  if (!MASTER(cr))
    *argc=0;
  gmx_sumi(1,argc,cr);
  
  if (!MASTER(cr))
    snew(argv_tmp,*argc+1);
  fprintf(stderr,"PID=%d argc=%d\n",cr->pid,*argc);
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
  debug_par();
}

t_commrec *init_par(int *argc,char ***argv_ptr)
{
  t_commrec *cr;
  char      **argv;
  int       i;
  
  argv = *argv_ptr;
  snew(cr,1);
  
  cr->nprocs=1;
  /* Get the number of processors.
   * This is useless for newer MPI versions.
   */
  for(i=0; (argv[i] != NULL); i++) {
    if (strcmp(argv[i],"-np")==0)
      if (argv[i+1]!=NULL)
	cr->nprocs=atoi(argv[i+1]);
  }

#ifdef USE_MPI
  cr->pid=mpiio_setup(argc,argv,&cr->nprocs);
#else
  cr->pid=0;
  if (cr->nprocs > 1) {
#ifdef USE_PVM3
    cr->pid=pvmio_setup(argv,cr->nprocs);
#else
    cr->pid=gmx_cpu_id();
#endif
  }
#endif
  
#ifdef _amb_
  pageMode(1);
#endif
  mem_init();

  if (!PAR(cr) && (cr->pid != 0))
    exit(0);
  
  if (PAR(cr)) {
    gmx_left_right(cr->nprocs,cr->pid,&cr->left,&cr->right);
#ifdef DEBUGPAR
    fprintf(stderr,"Going to initialise network\n");
#endif

#ifndef USE_PVM3
#ifdef DEBUGPAR
    fprintf(stderr,"Initialised network\n");
    fprintf(stderr,"Getting new processor id's\n");
#endif
    if (get_pid(stderr,cr->left,cr->right,&cr->pid,&cr->nprocs)==0)
      fatal_error(0,"could not get pid & nprocs from ring topology");
#ifdef DEBUGPAR
    fprintf(stderr,"Got new processor id's\n");
    fprintf(stderr,"nprocs=%d, pid=%d\n",cr->nprocs,cr->pid);
#endif

#endif
  }
  
  /* Communicate arguments if parallel */
  if (PAR(cr))
    comm_args(cr,argc,argv_ptr);
    
  return cr;
}

