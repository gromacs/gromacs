/*
 * Copyright (c) 1997 Massachusetts Institute of Technology
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to use, copy, modify, and distribute the Software without
 * restriction, provided the Software, including any modified copies made
 * under this license, is not distributed for a fee, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE MASSACHUSETTS INSTITUTE OF TECHNOLOGY BE LIABLE
 * FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
 * CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 * Except as contained in this notice, the name of the Massachusetts
 * Institute of Technology shall not be used in advertising or otherwise
 * to promote the sale, use or other dealings in this Software without
 * prior written authorization from the Massachusetts Institute of
 * Technology.
 *
 */

/*
 * malloc.c -- memory allocation related functions
 */

/* $Id$ */
#ifdef FFTW_USING_CILK
#include <cilk.h>
#include <cilk-compat.h>
#endif

#include <fftw.h>
#include <stdio.h>
#include <stdlib.h>

int fftw_malloc_cnt = 0;
void *(*fftw_malloc_hook) (size_t n) = (void *(*)(size_t n)) 0;
void (*fftw_free_hook) (void *p) = (void (*)(void *p)) 0;

#define FFTW_MALLOC_DEBUG 0
/* sorry for this debugging hack ... */
#define COMMA ,

#if FFTW_MALLOC_DEBUG
#define WHEN_DEBUG(a) a

/*
 * debugging malloc/free.  Initialize every malloced and freed area to
 * random values, just to make sure we are not using uninitialized
 * pointers.  Also check for writes past the ends of allocated blocks,
 * and a couple of other things.
 *
 * This code is a quick and dirty hack -- use at your own risk.
 */

int fftw_malloc_total = 0;

#define MAGIC 0xABadCafe
#define PAD_FACTOR 2
#define TWOINTS (2 * sizeof(int))

#define VERBOSE_ALLOCATION 0

#if VERBOSE_ALLOCATION
#define WHEN_VERBOSE(a) a
#else
#define WHEN_VERBOSE(a) 
#endif

void *fftw_malloc(size_t n)
{
     char *p;
     int i;

     WHEN_VERBOSE({
	  printf("FFTW_MALLOC %d\n",n);
	  fflush(stdout);
     })

     if (n == 0)
	  fftw_die("Tried to allocate a block of zero size!\n");

     fftw_malloc_total += n;

     p = (char *) malloc(PAD_FACTOR*n + TWOINTS);
     if (!p)
	  fftw_die("fftw_malloc: out of memory\n");

     /* store the size in a known position */
     ((int *) p)[0] = n;
     ((int *) p)[1] = MAGIC;
     for (i = 0; i < PAD_FACTOR*n; ++i)
	  p[i + TWOINTS] = (char) (i ^ 0xDEADBEEF);

     ++fftw_malloc_cnt;

     /* skip the size we stored previously */
     return (void *) (p + TWOINTS);
}

void fftw_free(void *p)
{
     char *q = ((char *) p) - TWOINTS;

     if (!p)
	  fftw_die("fftw_free: tried to free NULL pointer!\n");

     if (!q)
	  fftw_die("fftw_free: tried to free NULL+TWOINTS pointer!\n");

     {
	  int n = ((int *) q)[0];
	  int magic = ((int *) q)[1];
	  int i;
	  
	  WHEN_VERBOSE({
	       printf("FFTW_FREE %d\n",n);
	       fflush(stdout);
	  })
	  
	  if (n == 0)
	       fftw_die("Tried to free a freed pointer!\n");
	  *((int *) q) = 0; /* set to zero to detect duplicate free's */
	  
	  if (magic != MAGIC)
	       fftw_die("Wrong magic in fftw_free()!\n");	       
	  ((int *) q)[1] = ~MAGIC;

	  if (n < 0)
	       fftw_die("Tried to free block with corrupt size descriptor!\n");
	  
	  fftw_malloc_total -= n;
	  
	  if (fftw_malloc_total < 0)
	       fftw_die("fftw_malloc_total went negative!\n");
	  
	  /* check for writing past end of array: */
	  for (i = n; i < PAD_FACTOR*n; ++i)
	       if (q[i+TWOINTS] != (char) (i ^ 0xDEADBEEF)) {
		    fprintf(stderr, "Byte %d past end of array has changed!\n",
			    i - n + 1);
		    fftw_die("Array bounds overwritten!\n");
	       }
	  
	  for (i = 0; i < PAD_FACTOR*n; ++i)
	       q[i + TWOINTS] = (char) (i ^ 0xBEEFDEAD);

	  --fftw_malloc_cnt;
	  free(q);
     }
}

#else				/* production version, no hacks */
#define WHEN_DEBUG(a) 

void *fftw_malloc(size_t n)
{
     void *p;

     if (fftw_malloc_hook)
	  return fftw_malloc_hook(n);

     if (n == 0)
	  n = 1;

     p = malloc(n);

     if (!p)
	  fftw_die("fftw_malloc: out of memory\n");

     return p;
}

void fftw_free(void *p)
{
     if (p) {
	  if (fftw_free_hook) {
	       fftw_free_hook(p);
	       return;
	  }
	  free(p);
     }
}

#endif

/* die when fatal errors occur */
void fftw_die(char *s)
{
     fprintf(stderr, "%s", s);
     exit(1);
}

/* check for memory leaks when debugging */
void fftw_check_memory_leaks(void)
{
     extern int fftw_node_cnt, fftw_plan_cnt, fftw_twiddle_size;

     if (WHEN_DEBUG(fftw_malloc_cnt ||)
	 WHEN_DEBUG(fftw_malloc_total ||)
	 fftw_node_cnt || fftw_plan_cnt || fftw_twiddle_size) {
	  fprintf(stderr,
		  "MEMORY LEAK!!!\n"
		  WHEN_DEBUG("fftw_malloc = %d")
		  " node=%d plan=%d twiddle=%d\n"
		  WHEN_DEBUG("fftw_malloc_total = %d\n"), 
		  WHEN_DEBUG(fftw_malloc_cnt COMMA)
		  fftw_node_cnt, fftw_plan_cnt, fftw_twiddle_size
		  WHEN_DEBUG(COMMA fftw_malloc_total));
	  exit(1);
     }
}
