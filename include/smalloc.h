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
 * Gromacs Runs On Most of All Computer Systems
 */

#ifndef _smalloc_h
#define _smalloc_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>

/*
 * Memory allocation routines in gromacs:
 *
 * If an allocation fails, the program is halted by means of the
 * fatal_error routine, which outputs source file and line number
 * and the name of the variable involved.
 *
 * Macro's which can be used:
 *
 * snew(ptr,nelem)
 *    Allocates memory for nelem elements and returns this in ptr.
 *    The allocated memory is initialized to zeros.
 *
 * srenew(ptr,nelem)
 *    Reallocates memory for nelem elements and returns this in ptr.
 *
 * smalloc(ptr,size)
 *    Allocates memory for size bytes and returns this in ptr.
 *
 * scalloc(ptr,nelem,elsize)
 *    Allocates memory for nelem elements of size elsize and returns 
 *    this in ptr.
 *
 * srealloc(ptr,size)
 *    Reallocates memory for size bytes and returns this in ptr.
 *
 * sfree(ptr)
 *    Frees memory referenced by ptr.
 *
 ****************************************************************************
 *
 * Functions which are used by the macro's:
 *
 * extern void *save_malloc(char *name,char *file,int line,int size);
 *    Like alloc, returns a pointer to the allocated space, uses name, file
 *    and line to generate an error message when allocation failed.
 *
 * extern void *save_calloc(char *name,char *file,int line, 
 *                          size_t nelem,size_t elsize);
 *    Like calloc, returns a pointer to the allocated space, uses name, file
 *    and line to generate an error message when allocation failed.
 *
 * extern void *save_realloc(char *name,char *file,int line,
 *                           void *ptr,size_t size);
 *    Like realloc, returns a pointer to the allocated space, uses name, file
 *    and line to generate an error message when allocation failed.
 *    If ptr equals NULL, malloc is called in stead of realloc, in this way
 *    it is possible to combine first and later allocations.
 *
 * extern void save_free(char *name,char *file,int line, void *ptr);
 *    Like free, uses name, file and line to generate an error message when 
 *    the free failed.
 *
 * extern size_t maxavail();
 *    Returns the maximum available allocation unit, by applying a binary
 *    search on the largest block of memory available. After allocation
 *    it invokes free to restore the original state. So it is important
 *    that free can undo the effect of a malloc.
 * 
 * extern size_t memavail();
 *    Returns the total of available allocation unit, by applying maxavail
 *    until no space is left, it then frees all allocated space and returns
 *    the sum of the previously allocated space. As mentioned with maxavail,
 *    it is important that free can undo the effect of a malloc.
 * 
 */

#ifdef __cplusplus
extern "C" { 
#endif

void *save_malloc(const char *name,const char *file,int line,size_t size); 
void *save_calloc(const char *name,const char *file,int line,
		  size_t nelem,size_t elsize); 
void *save_realloc(const char *name,const char *file,int line,
		   void *ptr,size_t nelem,size_t elsize);
void save_free(const char *name,const char *file,int line, void *ptr);
size_t maxavail(void);
size_t memavail(void);

#ifdef __cplusplus
}

template <typename T>
void _snew(const char *name, T *&ptr, int nelem)
{
    ptr = (T *)save_calloc(name, __FILE__, __LINE__, nelem, sizeof(T));
}
template <typename T>
void _srenew(const char *name, T *&ptr, int nelem)
{
    ptr = (T *)save_realloc(name, __FILE__, __LINE__, ptr, nelem, sizeof(T));
}
template <typename T>
void _smalloc(const char *name, T *&ptr, size_t size)
{
    ptr = (T *)save_malloc(name, __FILE__, __LINE__, size);
}
template <typename T>
void _srealloc(const char *name, T *&ptr, size_t size)
{
    ptr = (T *)save_realloc(name, __FILE__, __LINE__, ptr, size, sizeof(char));
}



#define snew(ptr,nelem) _snew(#ptr,(ptr),(nelem))
#define srenew(ptr,nelem) _srenew(#ptr,(ptr),(nelem))
#define smalloc(ptr, size) _smalloc(#ptr,(ptr),(size))
#define srealloc(ptr, size) _srealloc(#ptr,(ptr),(size))

#else

/* These macros work in C, not in C++ */
#define snew(ptr,nelem) (ptr)=save_calloc(#ptr,__FILE__,__LINE__,\
			(nelem),sizeof(*(ptr)))
#define srenew(ptr,nelem) (ptr)=save_realloc(#ptr,__FILE__,__LINE__,\
			(ptr),(nelem),sizeof(*(ptr)))
#define smalloc(ptr,size) (ptr)=save_malloc(#ptr,__FILE__,__LINE__,size)
#define scalloc(ptr,nelem,elsize)\
		(ptr)=save_calloc(#ptr,__FILE__,__LINE__,nelem,elsize)
#define srealloc(ptr,size) (ptr)=save_realloc(#ptr,__FILE__,__LINE__,\
			(ptr),size,1)
#endif

#define sfree(ptr) save_free(#ptr,__FILE__,__LINE__,(ptr))

#endif	/* _smalloc_h */
