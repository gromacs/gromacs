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

/* If we don't have useful routines for allocating aligned memory,
 * then we have to use the old-style GROMACS approach bitwise-ANDing
 * pointers to ensure alignment. Freeing such a pointer requires
 * keeping track of the original pointer, so we set up an array
 * to store the original and aligned pointers. */

#if (!defined HAVE_POSIX_MEMALIGN && !defined HAVE_MEMALIGN)

#define SAVED_POINTERS_REALLOC 32

typedef struct
{
    void *freeable_ptr;
    void *aligned_ptr;
} saved_ptr_t;

static saved_ptr_t *saved_ptrs = NULL;
static int num_saved_ptrs = 0;

#endif

/* Pointers allocated with this routine should only be freed
 * with save_free_aligned, however this will only matter
 * on systems that lack posix_memalign() and memalign() when 
 * freeing memory that needed to be adjusted to achieve
 * the necessary alignment. */
void *save_calloc_aligned(char *name,char *file,int line,unsigned nelem,
                          unsigned elsize,unsigned alignment)
{
    void *p0,*p;
    bool allocate_fail;

    if (alignment == 0)
    {
        gmx_fatal(errno,__FILE__,__LINE__,
                  "Cannot allocate aligned memory with alignment of zero!\n(called from file %s, line %d)",file,line);
    }

#if (!defined HAVE_POSIX_MEMALIGN && !defined HAVE_MEMALIGN)
    if (0 == num_saved_ptrs % SAVED_POINTERS_REALLOC) {
#ifdef DEBUG
        log_action(0,name,file,line,0,0,ptr);
#endif
        srealloc(saved_ptrs, num_saved_ptrs + SAVED_POINTERS_REALLOC);
    }
#endif
    
    p0 = NULL;
    if ((nelem==0)||(elsize==0))
    {
        p  = NULL;
    }
    else
    {
#ifdef PRINT_ALLOC_KB
        if (nelem*elsize >= PRINT_ALLOC_KB*1024)
        {
            printf("Allocating %.1f MB for %s\n",
                   nelem*elsize/(PRINT_ALLOC_KB*1024.0),name);
        }
#endif

        allocate_fail = FALSE; /* stop compiler warnings */
#ifdef HAVE_POSIX_MEMALIGN
        allocate_fail = (0 != posix_memalign(&p, alignment, (size_t)nelem*(size_t)elsize));
#elif HAVE_MEMALIGN
        allocate_fail = ((p = memalign(alignment, (size_t)nelem*(size_t)elsize)) == NULL);
#else
        allocate_fail = ((p0 = malloc((size_t)nelem*(size_t)elsize+(size_t)alignment))==NULL);
#endif
        if (allocate_fail)
        {
            gmx_fatal(errno,__FILE__,__LINE__,
                      "Not enough memory. Failed to allocate %u aligned elements of size %u for %s\n(called from file %s, line %d)",nelem,elsize,name,file,line);
        }
  
#if (!defined HAVE_POSIX_MEMALIGN && !defined HAVE_MEMALIGN)
	/* Make the aligned pointer p, and save the underlying pointer that
	 * we're allowed to free(). */
	p = (void *) (((size_t) p0 + alignment - 1) & (~((size_t) (alignment-1))));
	saved_ptrs[num_saved_ptrs].freeable_ptr = p0;
	saved_ptrs[num_saved_ptrs].aligned_ptr = p;
	num_saved_ptrs++;
#endif

	memset(p, 0,(size_t) (nelem * elsize));
    }
    return p;
}

/* This routine can be called with any pointer */
void save_free_aligned(char *name,char *file,int line,void *ptr)
{
    int i, j;
    if (NULL != ptr)
    {
#if (!defined HAVE_POSIX_MEMALIGN && !defined HAVE_MEMALIGN)
        /* Manage the saved-pointers data structure. */
        for (i = num_saved_ptrs-1; i >= 0; i--)
        {
            if ((size_t) ptr == (size_t) saved_ptrs[i].aligned_ptr)
            {
	        ptr = saved_ptrs[i].freeable_ptr;
                /* Now remove the record of this saved pointer, replacing
                 * it with the one at the end of the array. */
		saved_ptrs[i] = saved_ptrs[num_saved_ptrs-1];
                num_saved_ptrs--;
                break;
            }
        }
#endif
        /* (Now) we're allowed to use a normal free() on this pointer. */
        save_free(name,file,line,ptr);
    }
}

