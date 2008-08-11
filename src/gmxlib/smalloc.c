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

/* This file is completely threadsafe - keep it that way! */
#include <gmx_thread.h>


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "gmx_fatal.h"
#include "smalloc.h"
#include "main.h"


#ifdef DEBUG
static void log_action(int bMal,char *what,char *file,int line,
                       int nelem,int size,void *ptr)
{
  static int btot=0;
  char *NN = "NULL";
  int        bytes;
  
  bytes=size*nelem;
  if (!bMal)
    bytes=-bytes;
  
#ifdef GMX_THREAD_PTHREAD
  pthread_mutex_lock(&gmx_logfile_mtx);
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
#ifdef GMX_THREAD_PTHREAD
  pthread_mutex_unlock(&gmx_logfile_mtx);
#endif
}
#endif

void *save_malloc(char *name,char *file,int line,unsigned size)
{
  void *p;
  
  p=NULL;
  if (size==0)
    p=NULL;
  else
    {
      if ((p=malloc(size))==NULL) 
        gmx_fatal(errno,__FILE__,__LINE__,
		  "Not enough memory. Failed to malloc %u bytes for %s\n(called from file %s, line %d)",
		  name,size,file,line);
      (void) memset(p,0,size);
    }
#ifdef DEBUG
  log_action(1,name,file,line,1,size,p);
#endif
  return p;
}

void *save_calloc(char *name,char *file,int line,
                  unsigned nelem,unsigned elsize)
{
  void *p;
  
  p=NULL;
  if ((nelem==0)||(elsize==0))
    p=NULL;
  else
    {
#ifdef PRINT_ALLOC_KB
      if (nelem*elsize >= PRINT_ALLOC_KB*1024) {
	printf("Allocating %.1f MB for %s\n",
	       nelem*elsize/(PRINT_ALLOC_KB*1024.0),name);
      }
#endif
#ifdef GMX_BROKEN_CALLOC
      /* emulate calloc(3) with malloc/memset on machines with 
         a broken calloc, e.g. in -lgmalloc on cray xt3. */
      if ((p=malloc((size_t)nelem*(size_t)elsize))==NULL) 
        gmx_fatal(errno,__FILE__,__LINE__,
		  "Not enough memory. Failed to calloc %u elements of size %u for %s\n(called from file %s, line %d)",nelem,elsize,name,file,line);
      memset(p, 0,(size_t) (nelem * elsize));
#else
      if ((p=calloc((size_t)nelem,(size_t)elsize))==NULL) 
        gmx_fatal(errno,__FILE__,__LINE__,
		  "Not enough memory. Failed to calloc %u elements of size %u for %s\n(called from file %s, line %d)",nelem,elsize,name,file,line);
#endif
    }
#ifdef DEBUG
  log_action(1,name,file,line,nelem,elsize,p);
#endif
  return p;
}

void *save_realloc(char *name,char *file,int line,void *ptr,
		   unsigned nelem,unsigned elsize)
{
  void *p;
  unsigned long size = nelem*elsize;
  
  p=NULL;
  if (size==0)
    p=NULL;
  else
    {
#ifdef PRINT_ALLOC_KB
      if (size >= PRINT_ALLOC_KB*1024) {
	printf("Reallocating %.1f MB for %s\n",
	       size/(PRINT_ALLOC_KB*1024.0),name);
      }
#endif
      if (ptr==NULL) 
	p=malloc((size_t)size); 
      else 
	p=realloc(ptr,(size_t)size);
      if (p == NULL) 
        gmx_fatal(errno,__FILE__,__LINE__,
		  "Not enough memory. Failed to realloc %u bytes for %s, %s=%x\n(called from file %s, line %d)",
		  size,name,name,ptr,file,line);
    }
#ifdef DEBUG
  log_action(1,name,file,line,1,size,p);
#endif
  return p;
}

void save_free(char *name,char *file,int line, void *ptr)
{
#ifdef DEBUG
  log_action(0,name,file,line,0,0,ptr);
#endif
  if (ptr != NULL)
    free(ptr);
}

unsigned maxavail(void)
{
  char *ptr;
  unsigned low,high,size;
  
  low=0;
  high=256e6;
  while ((high-low) > 4) {
    size=(high+low)/2;
    if ((ptr=malloc((size_t)size))==NULL)
      high=size;
    else {
      free(ptr);
      low=size;
    }
  }
  return low;
}

unsigned memavail(void)
{
  char *ptr;
  unsigned size;
  
  size = maxavail(); 
  if (size != 0) { 
    if ((ptr=malloc((size_t)size)) != NULL) {
      size += memavail();
      free(ptr);
    }
  }
  return size;
}
