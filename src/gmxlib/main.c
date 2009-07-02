/*
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
#include "mdrun.h"
#include "gmxfio.h"

/* The source code in this file should be thread-safe. 
         Please keep it that way. */


#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#if ((defined WIN32 || defined _WIN32 || defined WIN64 || defined _WIN64) && !defined __CYGWIN__ && !defined __CYGWIN32__)
#include <process.h>
#endif


#define BUFSIZE	1024

/* this is not strictly thread-safe, but it's only written to at the beginning
   of the simulation, once by each thread with the same value. We assume
   that writing to an int is atomic.*/
int  gmx_parallel_env=0;

static void par_fn(char *base,int ftp,const t_commrec *cr,
		   bool bUnderScore,
		   char buf[],int bufsize)
{
  int n;
  
  if(bufsize<(strlen(base)+4))
     gmx_mem("Character buffer too small!");

  /* Copy to buf, and strip extension */
  strcpy(buf,base);
  buf[strlen(base) - strlen(ftp2ext(fn2ftp(base))) - 1] = '\0';

  /* Add node info */
  if (bUnderScore)
    strcat(buf,"_");
  if (MULTISIM(cr) && !bUnderScore) {
    sprintf(buf+strlen(buf),"%d",cr->ms->sim);
  } else if (PAR(cr)) {
    sprintf(buf+strlen(buf),"%d",cr->nodeid);
  }
  strcat(buf,".");
  
  /* Add extension again */
  strcat(buf,(ftp == efTPX) ? "tpr" : (ftp == efEDR) ? "edr" : ftp2ext(ftp));
}

void check_multi_int(FILE *log,const gmx_multisim_t *ms,int val,
                     const char *name)
{
  int  *ibuf,p;
  bool bCompatible;

  fprintf(log,"Multi-checking %s ... ",name);
  
  if (ms == NULL)
    gmx_fatal(FARGS,
	      "check_multi_int called with a NULL communication pointer");

  snew(ibuf,ms->nsim);
  ibuf[ms->sim] = val;
  gmx_sumi_sim(ms->nsim,ibuf,ms);
  
  bCompatible = TRUE;
  for(p=1; p<ms->nsim; p++)
    bCompatible = bCompatible && (ibuf[p-1] == ibuf[p]);
  
  if (bCompatible) 
    fprintf(log,"OK\n");
  else {
    fprintf(log,"\n%s is not equal for all subsystems\n",name);
    for(p=0; p<ms->nsim; p++)
      fprintf(log,"  subsystem %d: %d\n",p,ibuf[p]);
    gmx_fatal(FARGS,"The %d subsystems are not compatible\n",ms->nsim);
  }
  
  sfree(ibuf);
}

FILE *gmx_log_open(char *lognm,const t_commrec *cr,bool bMasterOnly, 
                   unsigned long Flags)
{
  int  len,testlen,pid;
  char buf[256],host[256];
  time_t t;
  FILE *fp;

  bool bAppend = Flags & MD_APPENDFILES;	
  
  debug_gmx();
  
  /* Communicate the filename for logfile */
  if (cr->nnodes > 1 && !bMasterOnly) {
    if (MASTER(cr))
      len = strlen(lognm)+1;
    gmx_bcast(sizeof(len),&len,cr);
    if (!MASTER(cr))
      snew(lognm,len+8);
    gmx_bcast(len*sizeof(*lognm),lognm,cr);
  }
  
  debug_gmx();

  if (PAR(cr) && !bMasterOnly) {
    /* Since log always ends with '.log' let's use this info */
    par_fn(lognm,efLOG,cr,cr->ms!=NULL,buf,255);
	  fp = gmx_fio_fopen(buf, bAppend ? "a" : "w" );
  } else {
	  fp = gmx_fio_fopen(lognm, bAppend ? "a" : "w" );
  }

  gmx_fatal_set_log_file(fp);
  
  /* Get some machine parameters */
#ifdef HAVE_UNISTD_H
  if( gethostname(host,255) != 0)
    sprintf(host,"unknown");
#else
  sprintf(host,"unknown");
#endif  

  time(&t);

#ifndef NO_GETPID
#   if ((defined WIN32 || defined _WIN32 || defined WIN64 || defined _WIN64) && !defined __CYGWIN__ && !defined __CYGWIN32__)
	  pid = _getpid();
#   else
	  pid = getpid();
#   endif
#else
	pid = 0;
#endif

  if(bAppend)
  {
	  fprintf(fp,
			  "\n\n"
			  "-----------------------------------------------------------\n"
			  "Restarting from checkpoint, appending to previous log file.\n\n"
			  );
  }
	
  fprintf(fp,
	  "Log file opened on %s"
	  "Host: %s  pid: %d  nodeid: %d  nnodes:  %d\n",
	  ctime(&t),host,pid,cr->nodeid,cr->nnodes);

#if (defined BUILD_MACHINE && defined BUILD_TIME && defined BUILD_USER) 
  fprintf(fp,
	  "The Gromacs distribution was built %s by\n"
	  "%s (%s)\n\n\n",BUILD_TIME,BUILD_USER,BUILD_MACHINE);
#endif

  fflush(fp);
  debug_gmx();

  return fp;
}

void gmx_log_close(FILE *fp)
{
  if (fp) {
    gmx_fatal_set_log_file(NULL);
    gmx_fio_fclose(fp);
  }
}

static void comm_args(const t_commrec *cr,int *argc,char ***argv)
{
  int i,len;
  
  if ((cr) && PAR(cr))
    gmx_bcast(sizeof(*argc),argc,cr);
  
  if (!MASTER(cr))
    snew(*argv,*argc+1);
  fprintf(stderr,"NODEID=%d argc=%d\n",cr->nodeid,*argc);
  for(i=0; (i<*argc); i++) {
    if (MASTER(cr))
      len = strlen((*argv)[i])+1;
    gmx_bcast(sizeof(len),&len,cr);
    if (!MASTER(cr))
      snew((*argv)[i],len);
    /*gmx_bcast(len*sizeof((*argv)[i][0]),(*argv)[i],cr);*/
    gmx_bcast(len*sizeof(char),(*argv)[i],cr);
  }
  debug_gmx();
}

void init_multisystem(t_commrec *cr,int nsim,
		      int nfile,t_filenm fnm[],bool bParFn)
{
  gmx_multisim_t *ms;
  int  nnodes,nnodpersim,sim,i,ftp;
  char buf[256];
#ifdef GMX_MPI
  MPI_Group mpi_group_world;
#endif  
  int *rank;

  nnodes  = cr->nnodes;
  if (nnodes % nsim != 0)
    gmx_fatal(FARGS,"The number of nodes (%d) is not a multiple of the number of simulations (%d)",nnodes,nsim);

  nnodpersim = nnodes/nsim;
  sim = cr->nodeid/nnodpersim;

  if (debug)
    fprintf(debug,"We have %d simulations, %d nodes per simulation, local simulation is %d\n",nsim,nnodpersim,sim);

  snew(ms,1);
  cr->ms = ms;
  ms->nsim = nsim;
  ms->sim  = sim;
#ifdef GMX_MPI
  /* Create a communicator for the master nodes */
  snew(rank,ms->nsim);
  for(i=0; i<ms->nsim; i++)
    rank[i] = i*nnodpersim;
  MPI_Comm_group(MPI_COMM_WORLD,&mpi_group_world);
  MPI_Group_incl(mpi_group_world,nsim,rank,&ms->mpi_group_masters);
  sfree(rank);
  MPI_Comm_create(MPI_COMM_WORLD,ms->mpi_group_masters,
		  &ms->mpi_comm_masters);
#endif

  /* Reduce the intra-simulation communication */
  cr->sim_nodeid = cr->nodeid % nnodpersim;
  cr->nnodes = nnodpersim;
#ifdef GMX_MPI
  MPI_Comm_split(MPI_COMM_WORLD,sim,cr->sim_nodeid,&cr->mpi_comm_mysim);
  cr->mpi_comm_mygroup = cr->mpi_comm_mysim;
  cr->nodeid = cr->sim_nodeid;
#endif

  if (debug) {
    fprintf(debug,"This is simulation %d",cr->ms->sim);
    if (PAR(cr))
      fprintf(debug,", local number of nodes %d, local nodeid %d",
	      cr->nnodes,cr->sim_nodeid);
    fprintf(debug,"\n\n");
  }

  if (bParFn) {
    /* Patch output and tpx file names (except log which has been done already)
     */
    for(i=0; (i<nfile); i++) {
      /* Because of possible multiple extensions per type we must look 
       * at the actual file name 
       */
      if (is_output(&fnm[i]) ||
	  fnm[i].ftp == efTPX || fnm[i].ftp == efCPT ||
	  strcmp(fnm[i].opt,"-rerun") == 0) {
	ftp = fn2ftp(fnm[i].fns[0]);
	par_fn(fnm[i].fns[0],ftp,cr,FALSE,buf,255);
	sfree(fnm[i].fns[0]);
	fnm[i].fns[0] = strdup(buf);
      }
    }
  }
}

t_commrec *init_par(int *argc,char ***argv_ptr)
{
  t_commrec *cr;
  char      **argv;
  int       i;
  
  snew(cr,1);

  argv = *argv_ptr;

#ifdef GMX_MPI
#ifdef GMX_THREAD_MPI
  if (tMPI_Get_N(argc, argv_ptr)>1)
    gmx_parallel_env=1;
  else
    gmx_parallel_env=0;
#endif
#ifdef GMX_LIB_MPI
  gmx_parallel_env = 1;
#ifdef GMX_CHECK_MPI_ENV
  /* Do not use MPI calls when env.var. GMX_CHECK_MPI_ENV is not set */
  if (getenv(GMX_CHECK_MPI_ENV) == NULL)
    gmx_parallel_env = 0;
#endif
#endif
  if (gmx_parallel_env) {
    cr->sim_nodeid = gmx_setup(argc,argv,&cr->nnodes);
  } else {
    cr->nnodes     = 1;
    cr->sim_nodeid = 0;
  }
#else
  gmx_parallel_env = 0; 
  cr->sim_nodeid   = 0;
  cr->nnodes       = 1;
#endif

  if (!PAR(cr) && (cr->sim_nodeid != 0))
    gmx_comm("(!PAR(cr) && (cr->sim_nodeid != 0))");
  
  if (PAR(cr)) {
#ifdef GMX_MPI
    cr->mpi_comm_mysim = MPI_COMM_WORLD;
    cr->mpi_comm_mygroup = cr->mpi_comm_mysim;
#endif
  }
  cr->nodeid = cr->sim_nodeid;

  cr->duty = (DUTY_PP | DUTY_PME);

  /* Communicate arguments if parallel */
  if (PAR(cr))
    comm_args(cr,argc,argv_ptr);

  return cr;
}

t_commrec *init_cr_nopar(void)
{
  t_commrec *cr;

  snew(cr,1);

  cr->nnodes     = 1; 
  cr->sim_nodeid = 0;
  cr->nodeid     = 0;
  cr->nthreads   = 1;
  cr->threadid   = 0;
  cr->duty       = (DUTY_PP | DUTY_PME);

  return cr;
}
