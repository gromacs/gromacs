/*
 * Copyright (c) 1997-1999 Massachusetts Institute of Technology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
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

#include <fftw-int.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif

fftw_malloc_type_function fftw_malloc_hook = 0;
fftw_free_type_function fftw_free_hook = 0;
fftw_die_type_function fftw_die_hook = 0;

/**********************************************************
 *   DEBUGGING CODE
 **********************************************************/
#ifdef FFTW_DEBUG
static int fftw_malloc_cnt = 0;

/*
 * debugging malloc/free.  Initialize every malloced and freed area to
 * random values, just to make sure we are not using uninitialized
 * pointers.  Also check for writes past the ends of allocated blocks,
 * and a couple of other things.
 *
 * This code is a quick and dirty hack -- use at your own risk.
 */

static int fftw_malloc_total = 0, fftw_malloc_max = 0, fftw_malloc_cnt_max = 0;

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

     fftw_malloc_total += n;

     if (fftw_malloc_total > fftw_malloc_max)
	  fftw_malloc_max = fftw_malloc_total;

     p = (char *) malloc(PAD_FACTOR * n + TWOINTS);
     if (!p)
	  fftw_die("fftw_malloc: out of memory\n");

     /* store the size in a known position */
     ((int *) p)[0] = n;
     ((int *) p)[1] = MAGIC;
     for (i = 0; i < PAD_FACTOR * n; ++i)
	  p[i + TWOINTS] = (char) (i ^ 0xDEADBEEF);

     ++fftw_malloc_cnt;

     if (fftw_malloc_cnt > fftw_malloc_cnt_max)
	  fftw_malloc_cnt_max = fftw_malloc_cnt;

     /* skip the size we stored previously */
     return (void *) (p + TWOINTS);
}

void fftw_free(void *p)
{
     char *q;

     if (!p)
	  return;

     q = ((char *) p) - TWOINTS;
     if (!q)
	  fftw_die("fftw_free: tried to free NULL+TWOINTS pointer!\n");

     {
	  int n = ((int *) q)[0];
	  int magic = ((int *) q)[1];
	  int i;

	  WHEN_VERBOSE( {
		       printf("FFTW_FREE %d\n", n);
		       fflush(stdout);
		       })

	  *((int *) q) = 0;	/* set to zero to detect duplicate free's */

	  if (magic != MAGIC)
	       fftw_die("Wrong magic in fftw_free()!\n");
	  ((int *) q)[1] = ~MAGIC;

	  if (n < 0)
	       fftw_die("Tried to free block with corrupt size descriptor!\n");

	  fftw_malloc_total -= n;

	  if (fftw_malloc_total < 0)
	       fftw_die("fftw_malloc_total went negative!\n");

	  /* check for writing past end of array: */
	  for (i = n; i < PAD_FACTOR * n; ++i)
	       if (q[i + TWOINTS] != (char) (i ^ 0xDEADBEEF)) {
		    fflush(stdout);
		    fprintf(stderr, "Byte %d past end of array has changed!\n",
			    i - n + 1);
		    fftw_die("Array bounds overwritten!\n");
	       }
	  for (i = 0; i < PAD_FACTOR * n; ++i)
	       q[i + TWOINTS] = (char) (i ^ 0xBEEFDEAD);

	  --fftw_malloc_cnt;

	  if (fftw_malloc_cnt < 0)
	       fftw_die("fftw_malloc_cnt went negative!\n");

	  if (fftw_malloc_cnt == 0 && fftw_malloc_total > 0 ||
	      fftw_malloc_cnt > 0 && fftw_malloc_total == 0)
	       fftw_die("fftw_malloc_cnt/total not zero at the same time!\n");

	  free(q);
     }
}

#else
/**********************************************************
 *   NON DEBUGGING CODE
 **********************************************************/
/* production version, no hacks */

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
void fftw_die(const char *s)
{
     if (fftw_die_hook)
	  fftw_die_hook(s);

     fflush(stdout);
     fprintf(stderr, "fftw: %s", s);
     exit(EXIT_FAILURE);
}

/* check for memory leaks when debugging */
void fftw_check_memory_leaks(void)
{
     extern int fftw_node_cnt, fftw_plan_cnt, fftw_twiddle_size;

#ifdef FFTW_DEBUG
     if (fftw_malloc_cnt || fftw_malloc_total ||
	 fftw_node_cnt || fftw_plan_cnt || fftw_twiddle_size) {
	  fflush(stdout);
	  fprintf(stderr,
		  "MEMORY LEAK!!!\n"
		  "fftw_malloc = %d"
		  " node=%d plan=%d twiddle=%d\n"
		  "fftw_malloc_total = %d\n",
		  fftw_malloc_cnt,
		  fftw_node_cnt, fftw_plan_cnt, fftw_twiddle_size,
		  fftw_malloc_total);
	  exit(EXIT_FAILURE);
     }
#else
     if (fftw_node_cnt || fftw_plan_cnt || fftw_twiddle_size) {
	  fflush(stdout);
	  fprintf(stderr,
		  "MEMORY LEAK!!!\n"
		  " node=%d plan=%d twiddle=%d\n",
		  fftw_node_cnt, fftw_plan_cnt, fftw_twiddle_size);
	  exit(EXIT_FAILURE);
     }
#endif
}

void fftw_print_max_memory_usage(void)
{
#ifdef FFTW_DEBUG
     printf("\nMaximum number of blocks allocated = %d\n"
	    "Maximum number of bytes allocated  = %0.3f kB\n",
	    fftw_malloc_cnt_max, fftw_malloc_max / 1024.0);
#endif
}
