/*
 *       @(#) copyrgt.c 1.12 9/30/97
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0b
 * 
 * Copyright (c) 1990-1997,
 * BIOSON Research Institute, Dept. of Biophysical Chemistry,
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
 * Gromacs Runs On Most of All Computer Systems
 */

#ifndef	_memtab_h
#define	_memtab_h

#ifdef HAVE_IDENT
#ident	"@(#) memtab.h 1.12 12/16/92"
#endif /* HAVE_IDENT */

/*
 * This module is intended for alloc(at)ing memory in one contiguous
 * block, using only local references. All references within this block
 * are made relative to the start of the block. This means that before
 * using a pointer in this block, they must be reallocated, simply by
 * adding the base address to the pointer value.
 */

#define	NOENTRY		-1	/* Denotes a NULL pointer reference	*/

typedef struct
{
  void *ref; 		/* The "physical" address of an entry		*/
  void *handle;		/* The associated handle			*/
} t_mref;

typedef struct
{
  int nref;		/* The number of inserted references		*/
  t_mref *refs;		/* The inserted references and their handles	*/
  int msize;		/* The total size of the memory allocated sofar	*/
  char *mtab;		/* The allocated memory (one contiguous block)	*/
} t_memtab;

extern void init_memtab(t_memtab *mtab);
     /*
      * Initialises the struct, should be called before any other action
      * The function init_memtab() will initialise a struct which can be
      * extended by successive invokations of put_memtab().
      */

extern void dispose_memtab(t_memtab *mtab);
     /*
      * Disposes all the memory allocated for the struct. After this
      * any reference within the block may become invalid (depends on
      * the definition of the free() call).
      */

extern void *put_memtab(t_memtab *mtab,int size,void *ref);
     /*
      * The function put_memtab() returns a handle to the memory block
      * specified by pointer ref and size bytes. If the address was inserted
      * before in the struct, this (previous) handle will be returned else
      * the struct will be extended, a new handle will be created and the
      * data will be copied into the struct. Note that the returned handle
      * is actually the offset of the copied block within the struct.
      * NULL pointers and null sizes will return a handle that will convert 
      * to NULL after invokation of get_memtab(). Note that after every
      * invokation absolute addresses, got before from get_memtab() might
      * have changed due to a move within the reallocation.
      */

extern void *get_memtab(t_memtab *mtab,void *handle);
     /*
      * Returns the (physical) address corresponding to the specified handle.
      * Handle should be a value returned by put_memtab(). When invoked
      * with NULL for handle, get_memtab() returns the base address of
      * the allocated block.
      */

long wr_memtab(FILE *fp,t_memtab *mtab);
     /*
      * Writes the referenced contents of the struct to the file specified
      * by fp. The function wr_memtab() returnes the number of bytes written.
      */

void *rd_memtab(FILE *fp,int size,t_memtab *mtab);
     /*
      * Extends the struct by reading size bytes from the file specified
      * by fp. The function rd_memtab() will return a handle to the read
      * block.
      */

#endif	/* _memtab_h */
