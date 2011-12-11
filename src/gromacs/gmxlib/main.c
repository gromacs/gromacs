/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
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

#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
#endif

#ifdef HAVE_DIRECT_H
/* windows-specific include for _chdir() */
#include <direct.h>
#endif


#include "smalloc.h"
#include "gmx_fatal.h"
#include "network.h"
#include "main.h"
#include "macros.h"
#include "futil.h"
#include "filenm.h"
#include "mdrun.h"
#include "gmxfio.h"
#include "string2.h"

#ifdef GMX_THREADS
#include "thread_mpi.h"
#endif

/* The source code in this file should be thread-safe. 
         Please keep it that way. */


#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#if ((defined WIN32 || defined _WIN32 || defined WIN64 || defined _WIN64) && !defined __CYGWIN__ && !defined __CYGWIN32__)
#include <process.h>
#endif


/* Portable version of ctime_r implemented in src/gmxlib/string2.c, but we do not want it declared in public installed headers */
char *
gmx_ctime_r(const time_t *clock,char *buf, int n);


#define BUFSIZE	1024

/* this is not strictly thread-safe, but it's only written to at the beginning
   of the simulation, once by each thread with the same value. We assume
   that writing to an int is atomic.*/
static gmx_bool parallel_env_val;
#ifdef GMX_THREADS
tMPI_Thread_mutex_t parallel_env_mutex=TMPI_THREAD_MUTEX_INITIALIZER;
#endif


/* returns 1 when running in a parallel environment, so could also be 1 if
   mdrun was started with: mpirun -np 1.
     
   Use this function only to check whether a parallel environment has   
   been initialized, for example when checking whether gmx_finalize()   
   needs to be called. Use PAR(cr) to check whether the simulation actually
   has more than one node/thread.  */
gmx_bool gmx_parallel_env_initialized(void)
{
    gmx_bool ret;
#ifdef GMX_THREADS
    tMPI_Thread_mutex_lock(&parallel_env_mutex);
#endif
    ret=parallel_env_val;
#ifdef GMX_THREADS
    tMPI_Thread_mutex_unlock(&parallel_env_mutex);
#endif
    return ret;
}

static void set_parallel_env(gmx_bool val)
{
#ifdef GMX_THREADS
    tMPI_Thread_mutex_lock(&parallel_env_mutex);
#endif
    if (!parallel_env_val)
    {
        /* we only allow it to be set, not unset */
        parallel_env_val=val;
    }
#ifdef GMX_THREADS
    tMPI_Thread_mutex_unlock(&parallel_env_mutex);
#endif
}


static void par_fn(char *base,int ftp,const t_commrec *cr,
		   gmx_bool bAppendSimId,gmx_bool bAppendNodeId,
		   char buf[],int bufsize)
{
  int n;
  
  if((size_t)bufsize<(strlen(base)+10))
     gmx_mem("Character buffer too small!");

  /* Copy to buf, and strip extension */
  strcpy(buf,base);
  buf[strlen(base) - strlen(ftp2ext(fn2ftp(base))) - 1] = '\0';

  if (bAppendSimId) {
    sprintf(buf+strlen(buf),"%d",cr->ms->sim);
  }
  if (bAppendNodeId) {
    strcat(buf,"_node");
    sprintf(buf+strlen(buf),"%d",cr->nodeid);
  }
  strcat(buf,".");
  
  /* Add extension again */
  strcat(buf,(ftp == efTPX) ? "tpr" : (ftp == efEDR) ? "edr" : ftp2ext(ftp));
  if (cr->nodeid == 0) {
    printf("node %d par_fn '%s'\n",cr->nodeid,buf);
    if (fn2ftp(buf) == efLOG) {
      printf("log\n");
    }
  }
}

void check_multi_int(FILE *log,const gmx_multisim_t *ms,int val,
                     const char *name)
{
  int  *ibuf,p;
  gmx_bool bCompatible;

  if (NULL != log)
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
  {
      if (NULL != log)
          fprintf(log,"OK\n");
  }
  else 
  {
      if (NULL != log)
      {
          fprintf(log,"\n%s is not equal for all subsystems\n",name);
          for(p=0; p<ms->nsim; p++)
              fprintf(log,"  subsystem %d: %d\n",p,ibuf[p]);
      }
      gmx_fatal(FARGS,"The %d subsystems are not compatible\n",ms->nsim);
  }
  
  sfree(ibuf);
}

void check_multi_large_int(FILE *log,const gmx_multisim_t *ms,
                           gmx_large_int_t val, const char *name)
{
  gmx_large_int_t  *ibuf;
  int p;
  gmx_bool bCompatible;

  if (NULL != log)
      fprintf(log,"Multi-checking %s ... ",name);
  
  if (ms == NULL)
    gmx_fatal(FARGS,
	      "check_multi_int called with a NULL communication pointer");

  snew(ibuf,ms->nsim);
  ibuf[ms->sim] = val;
  gmx_sumli_sim(ms->nsim,ibuf,ms);
  
  bCompatible = TRUE;
  for(p=1; p<ms->nsim; p++)
    bCompatible = bCompatible && (ibuf[p-1] == ibuf[p]);
  
  if (bCompatible) 
  {
      if (NULL != log)
          fprintf(log,"OK\n");
  }
  else 
  {
      if (NULL != log)
      {
          fprintf(log,"\n%s is not equal for all subsystems\n",name);
          for(p=0; p<ms->nsim; p++)
          {
              char strbuf[255];
              /* first make the format string */
              snprintf(strbuf, 255, "  subsystem %%d: %s\n", 
                       gmx_large_int_pfmt);
              fprintf(log,strbuf,p,ibuf[p]);
          }
      }
      gmx_fatal(FARGS,"The %d subsystems are not compatible\n",ms->nsim);
  }
  
  sfree(ibuf);
}


void gmx_log_open(const char *lognm,const t_commrec *cr,gmx_bool bMasterOnly, 
                   unsigned long Flags, FILE** fplog)
{
    int  len,testlen,pid;
    char buf[256],host[256];
    time_t t;
    char timebuf[STRLEN];
    FILE *fp=*fplog;
    char *tmpnm;

    gmx_bool bAppend = Flags & MD_APPENDFILES;	
  
    debug_gmx();
  
    /* Communicate the filename for logfile */
    if (cr->nnodes > 1 && !bMasterOnly
#ifdef GMX_THREADS
        /* With thread MPI the non-master log files are opened later
         * when the files names are already known on all nodes.
         */
        && FALSE
#endif
        )
    {
        if (MASTER(cr))
        {
            len = strlen(lognm) + 1;
        }
        gmx_bcast(sizeof(len),&len,cr);
        if (!MASTER(cr))
        {
            snew(tmpnm,len+8);
        }
        else
        {
            tmpnm=gmx_strdup(lognm);
        }
        gmx_bcast(len*sizeof(*tmpnm),tmpnm,cr);
    }
    else
    {
        tmpnm=gmx_strdup(lognm);
    }
  
    debug_gmx();

    if (!bMasterOnly && !MASTER(cr))
    {
        /* Since log always ends with '.log' let's use this info */
        par_fn(tmpnm,efLOG,cr,FALSE,!bMasterOnly,buf,255);
        fp = gmx_fio_fopen(buf, bAppend ? "a+" : "w+" );
    }
    else if (!bAppend)
    {
        fp = gmx_fio_fopen(tmpnm, bAppend ? "a+" : "w+" );
    }

    sfree(tmpnm);

    gmx_fatal_set_log_file(fp);
  
    /* Get some machine parameters */
#ifdef HAVE_UNISTD_H
    if (gethostname(host,255) != 0)
    {
        sprintf(host,"unknown");
    }
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

    if (bAppend)
    {
        fprintf(fp,
                "\n"
                "\n"
                "-----------------------------------------------------------\n"
                "Restarting from checkpoint, appending to previous log file.\n"
                "\n"
            );
    }
	
    gmx_ctime_r(&t,timebuf,STRLEN);

    fprintf(fp,
            "Log file opened on %s"
            "Host: %s  pid: %d  nodeid: %d  nnodes:  %d\n",
            timebuf,host,pid,cr->nodeid,cr->nnodes);

#if (defined BUILD_MACHINE && defined BUILD_TIME && defined BUILD_USER) 
    fprintf(fp,
            "The Gromacs distribution was built %s by\n"
            "%s (%s)\n\n\n",BUILD_TIME,BUILD_USER,BUILD_MACHINE);
#endif

    fflush(fp);
    debug_gmx();

    *fplog = fp;
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

void init_multisystem(t_commrec *cr,int nsim, char **multidirs,
                      int nfile, const t_filenm fnm[],gmx_bool bParFn)
{
    gmx_multisim_t *ms;
    int  nnodes,nnodpersim,sim,i,ftp;
    char buf[256];
#ifdef GMX_MPI
    MPI_Group mpi_group_world;
#endif  
    int *rank;

#ifndef GMX_MPI
    if (nsim > 1)
    {
        gmx_fatal(FARGS,"This binary is compiled without MPI support, can not do multiple simulations.");
    }
#endif

    nnodes  = cr->nnodes;
    if (nnodes % nsim != 0)
    {
        gmx_fatal(FARGS,"The number of nodes (%d) is not a multiple of the number of simulations (%d)",nnodes,nsim);
    }

    nnodpersim = nnodes/nsim;
    sim = cr->nodeid/nnodpersim;

    if (debug)
    {
        fprintf(debug,"We have %d simulations, %d nodes per simulation, local simulation is %d\n",nsim,nnodpersim,sim);
    }

    snew(ms,1);
    cr->ms = ms;
    ms->nsim = nsim;
    ms->sim  = sim;
#ifdef GMX_MPI
    /* Create a communicator for the master nodes */
    snew(rank,ms->nsim);
    for(i=0; i<ms->nsim; i++)
    {
        rank[i] = i*nnodpersim;
    }
    MPI_Comm_group(MPI_COMM_WORLD,&mpi_group_world);
    MPI_Group_incl(mpi_group_world,nsim,rank,&ms->mpi_group_masters);
    sfree(rank);
    MPI_Comm_create(MPI_COMM_WORLD,ms->mpi_group_masters,
                    &ms->mpi_comm_masters);

#if !defined(GMX_THREADS) && !defined(MPI_IN_PLACE_EXISTS)
    /* initialize the MPI_IN_PLACE replacement buffers */
    snew(ms->mpb, 1);
    ms->mpb->ibuf=NULL;
    ms->mpb->libuf=NULL;
    ms->mpb->fbuf=NULL;
    ms->mpb->dbuf=NULL;
    ms->mpb->ibuf_alloc=0;
    ms->mpb->libuf_alloc=0;
    ms->mpb->fbuf_alloc=0;
    ms->mpb->dbuf_alloc=0;
#endif

#endif

    /* Reduce the intra-simulation communication */
    cr->sim_nodeid = cr->nodeid % nnodpersim;
    cr->nnodes = nnodpersim;
#ifdef GMX_MPI
    MPI_Comm_split(MPI_COMM_WORLD,sim,cr->sim_nodeid,&cr->mpi_comm_mysim);
    cr->mpi_comm_mygroup = cr->mpi_comm_mysim;
    cr->nodeid = cr->sim_nodeid;
#endif

    if (debug)
    {
        fprintf(debug,"This is simulation %d",cr->ms->sim);
        if (PAR(cr))
        {
            fprintf(debug,", local number of nodes %d, local nodeid %d",
                    cr->nnodes,cr->sim_nodeid);
        }
        fprintf(debug,"\n\n");
    }

    if (multidirs)
    {
        int ret;
        if (debug)
        {
            fprintf(debug,"Changing to directory %s\n",multidirs[cr->ms->sim]);
        }
        if (chdir(multidirs[cr->ms->sim]) != 0)
        {
            gmx_fatal(FARGS, "Couldn't change directory to %s: %s",
                      multidirs[cr->ms->sim],
                      strerror(errno));
        }
    }
    else if (bParFn)
    {
        /* Patch output and tpx, cpt and rerun input file names */
        for(i=0; (i<nfile); i++)
        {
            /* Because of possible multiple extensions per type we must look 
             * at the actual file name 
             */
            if (is_output(&fnm[i]) ||
                fnm[i].ftp == efTPX || fnm[i].ftp == efCPT ||
                strcmp(fnm[i].opt,"-rerun") == 0)
            {
                ftp = fn2ftp(fnm[i].fns[0]);
                par_fn(fnm[i].fns[0],ftp,cr,TRUE,FALSE,buf,255);
                sfree(fnm[i].fns[0]);
                fnm[i].fns[0] = gmx_strdup(buf);
            }
        }
    }
}

t_commrec *init_par(int *argc,char ***argv_ptr)
{
    t_commrec *cr;
    char      **argv;
    int       i;
    gmx_bool      pe=FALSE;

    snew(cr,1);

    argv = *argv_ptr;

#ifdef GMX_MPI
#ifdef GMX_LIB_MPI
    pe = TRUE;
#ifdef GMX_CHECK_MPI_ENV
    /* Do not use MPI calls when env.var. GMX_CHECK_MPI_ENV is not set */
    if (getenv(GMX_CHECK_MPI_ENV) == NULL)
        pe = FALSE;
#endif /* GMX_CHECK_MPI_ENV */
#endif /* GMX_LIB_MPI  */
    set_parallel_env(pe);
    if (pe) {
        cr->sim_nodeid = gmx_setup(argc,argv,&cr->nnodes);
    } else {
        cr->nnodes     = 1;
        cr->sim_nodeid = 0;
    }
#else /* GMX_MPI */
    pe=FALSE;
    set_parallel_env(pe);
    cr->sim_nodeid   = 0;
    cr->nnodes       = 1;
#endif /* GMX_MPI */

    if (!PAR(cr) && (cr->sim_nodeid != 0))
        gmx_comm("(!PAR(cr) && (cr->sim_nodeid != 0))");

    if (PAR(cr)) 
    {
#ifdef GMX_MPI
        cr->mpi_comm_mysim = MPI_COMM_WORLD;
        cr->mpi_comm_mygroup = cr->mpi_comm_mysim;
#endif /* GMX_MPI */
    }
    cr->nodeid = cr->sim_nodeid;

    cr->duty = (DUTY_PP | DUTY_PME);

    /* Communicate arguments if parallel */
#ifndef GMX_THREADS
    if (PAR(cr))
        comm_args(cr,argc,argv_ptr);
#endif /* GMX_THREADS */

#ifdef GMX_MPI
#if !defined(GMX_THREADS) && !defined(MPI_IN_PLACE_EXISTS)
  /* initialize the MPI_IN_PLACE replacement buffers */
  snew(cr->mpb, 1);
  cr->mpb->ibuf=NULL;
  cr->mpb->libuf=NULL;
  cr->mpb->fbuf=NULL;
  cr->mpb->dbuf=NULL;
  cr->mpb->ibuf_alloc=0;
  cr->mpb->libuf_alloc=0;
  cr->mpb->fbuf_alloc=0;
  cr->mpb->dbuf_alloc=0;
#endif
#endif

    return cr;
}

t_commrec *init_par_threads(const t_commrec *cro)
{
#ifdef GMX_THREADS
    int initialized;
    t_commrec *cr;

    /* make a thread-specific commrec */
    snew(cr,1);
    /* now copy the whole thing, so settings like the number of PME nodes
       get propagated. */
    *cr=*cro;

    /* and we start setting our own thread-specific values for things */
    MPI_Initialized(&initialized);
    if (!initialized)
        gmx_comm("Initializing threads without comm");
    set_parallel_env(TRUE);
    /* once threads will be used together with MPI, we'll
       fill the cr structure with distinct data here. This might even work: */
    cr->sim_nodeid = gmx_setup(0,NULL, &cr->nnodes);

    cr->mpi_comm_mysim = MPI_COMM_WORLD;
    cr->mpi_comm_mygroup = cr->mpi_comm_mysim;
    cr->nodeid = cr->sim_nodeid;
    cr->duty = (DUTY_PP | DUTY_PME);

    return cr;
#else
    return NULL;
#endif
}


t_commrec *init_cr_nopar(void)
{
    t_commrec *cr;

    snew(cr,1);

    cr->nnodes     = 1; 
    /* cr->nthreads   = 1; */
    cr->sim_nodeid = 0;
    cr->nodeid     = 0;
    /* cr->threadid   = 0; */
    cr->duty       = (DUTY_PP | DUTY_PME);

    return cr;
}
