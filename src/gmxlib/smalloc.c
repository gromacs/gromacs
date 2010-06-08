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

/* This file is completely threadsafe - keep it that way! */

#ifdef GMX_THREADS
#include "thread_mpi/threads.h"
#endif 


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "gmx_fatal.h"
#include "smalloc.h"
#include "main.h"
#ifdef WITH_DMALLOC
#include "dmalloc.h"
#endif

#ifdef DEBUG
static void log_action(int bMal,const char *what,const char *file,int line,
                       int nelem,int size,void *ptr)
{
  static int btot=0;
  char *NN = "NULL";
  int        bytes;
  
  bytes=size*nelem;
  if (!bMal)
    bytes=-bytes;
  
#ifdef GMX_THREADS
  tMPI_Thread_mutex_lock(&gmx_logfile_mtx);
#endif

  /* This total memory count is not correct, since with realloc
   * it adds the whole size again, not just the increment.
   */
  /* This static variable is protected by the mutex too... */
  btot+=bytes;
    
  bytes/=1024;
  if (debug && (bytes != 0)) {
    fprintf(debug,"%s:%d kB (%7d kB) [%s, line %d, nelem %d, size %d]\n",
	    what ? what : NN,bytes,btot/1024,
	    file ? file : NN,line,nelem,size);
  }
  /* Print to stderr for things larger than 1 MB */
  if (bytes >= 1024 || bytes <= -1024) {
    char *fname=NULL;
    if (file) {
      fname = strrchr(file,'/');
      if (fname) {
	fname++;
      } else {
	fname = file;
      }
    }
    printf("%s: %.1f MB [%s, line %d, nelem %d, size %d]\n",
	   what ? what  : NN,bytes/1024.0,
	   file ? fname : NN,line,nelem,size);
  }
#ifdef GMX_THREADS
  tMPI_Thread_mutex_unlock(&gmx_logfile_mtx);
#endif
}
#endif

void *save_malloc(const char *name,const char *file,int line,size_t size)
{
  void *p;
  
  p=NULL;
  if (size==0)
    p=NULL;
  else
    {
      if ((p=malloc(size))==NULL) 
        gmx_fatal(errno,__FILE__,__LINE__,
		  "Not enough memory. Failed to malloc %s bytes for %" gmx_large_int_fmt "\n(called from file %s, line %d)",
		  name,(gmx_large_int_t)size,file,line);
      (void) memset(p,0,size);
    }
#ifdef DEBUG
  log_action(1,name,file,line,1,size,p);
#endif
  return p;
}

void *save_calloc(const char *name,const char *file,int line,
                  size_t nelem,size_t elsize)
{
  void *p;
  
  p=NULL;
  if ((nelem==0)||(elsize==0))
    p=NULL;
  else
    {
#ifdef PRINT_ALLOC_KB
      int rank=0;
      if (nelem*elsize >= PRINT_ALLOC_KB*1024) {
#ifdef GMX_MPI
#include <mpi.h>
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#endif
	printf("Allocating %.1f MB for %s (called from file %s, line %d on %d)\n",
	       nelem*elsize/1048576.0,name,file,line,rank);
      }
#endif
#ifdef GMX_BROKEN_CALLOC
      /* emulate calloc(3) with malloc/memset on machines with 
         a broken calloc, e.g. in -lgmalloc on cray xt3. */
      if ((p=malloc((size_t)nelem*(size_t)elsize))==NULL) 
        gmx_fatal(errno,__FILE__,__LINE__,
		  "Not enough memory. Failed to calloc %"gmx_large_int_fmt
		  " elements of size %"gmx_large_int_fmt
		  " for %s\n(called from file %s, line %d)",
		  (gmx_large_int_t)nelem,(gmx_large_int_t)elsize,
		  name,file,line);
      memset(p, 0,(size_t) (nelem * elsize));
#else
      if ((p=calloc((size_t)nelem,(size_t)elsize))==NULL) 
        gmx_fatal(errno,__FILE__,__LINE__,
		  "Not enough memory. Failed to calloc %"gmx_large_int_fmt
		  " elements of size %"gmx_large_int_fmt
		  " for %s\n(called from file %s, line %d)",
		  (gmx_large_int_t)nelem,(gmx_large_int_t)elsize,name,file,line);
#endif
    }
#ifdef DEBUG
  log_action(1,name,file,line,nelem,elsize,p);
#endif
  return p;
}

void *save_realloc(const char *name,const char *file,int line,void *ptr,
		   size_t nelem,size_t elsize)
{
  void *p;
  size_t size = nelem*elsize;
  
  p=NULL;
  if (size==0)
    {
      save_free(name, file, line, ptr);
    }
  else
    {
#ifdef PRINT_ALLOC_KB
      int rank=0;
      if (size >= PRINT_ALLOC_KB*1024) {
#ifdef GMX_MPI
#include <mpi.h>
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#endif
	printf("Reallocating %.1f MB for %s (called from file %s, line %d on %d)\n",
	       size/1048576.0,name,file,line,rank);
      }
#endif
      if (ptr==NULL) 
	p=malloc((size_t)size); 
      else 
	p=realloc(ptr,(size_t)size);
      if (p == NULL) 
        gmx_fatal(errno,__FILE__,__LINE__,
		  "Not enough memory. Failed to realloc %"gmx_large_int_fmt
		  " bytes for %s, %s=%x\n(called from file %s, line %d)",
		  (gmx_large_int_t)size,name,name,ptr,file,line);
#ifdef DEBUG
      log_action(1,name,file,line,1,size,p);
#endif
    }
  return p;
}

void save_free(const char *name,const char *file,int line, void *ptr)
{
#ifdef DEBUG
  log_action(0,name,file,line,0,0,ptr);
#endif
  if (ptr != NULL)
    free(ptr);
}

size_t maxavail(void)
{
  char *ptr;
  size_t low,high,size;
  
  low=0;
  high=256e6;
  while ((high-low) > 4) {
    size=(high+low)/2;
    if ((ptr=(char *)malloc((size_t)size))==NULL)
      high=size;
    else {
      free(ptr);
      low=size;
    }
  }
  return low;
}

size_t memavail(void)
{
  char *ptr;
  size_t size;
  
  size = maxavail(); 
  if (size != 0) { 
    if ((ptr=(char *)malloc((size_t)size)) != NULL) {
      size += memavail();
      free(ptr);
    }
  }
  return size;
}
