
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
 * twiddle.c -- compute twiddle factors
 * These are the twiddle factors for *direct* fft.  Flip sign to get
 * the inverse
 */

/* $Id$ */
#ifdef FFTW_USING_CILK
#include <cilk.h>
#include <cilk-compat.h>
#endif

#include <fftw.h>
#include <math.h>
#include <stdlib.h>

#define FFTW_K2PI 6.2831853071795864769252867665590057683943387987502

/*
 * compute the W coefficients (that is, powers of the root of 1)
 * and store them into an array.
 */
static void fftw_compute_twiddle(int n, int r, int m, FFTW_COMPLEX *W)
{
     double twoPiOverN;
     int i, j;

     twoPiOverN = FFTW_K2PI / (double) n;
     for (i = 0; i < m; ++i)
	  for (j = 1; j < r; ++j) {
	       int k = i * (r - 1) + (j - 1);
	       c_re(W[k]) = cos(twoPiOverN * (double) i * (double) j);
	       c_im(W[k]) = -sin(twoPiOverN * (double) i * (double) j);
	  }
}

/*
 * these routines implement a simple reference-count-based 
 * management of twiddle structures
 */
static fftw_twiddle *twlist = (fftw_twiddle *) 0;
int fftw_twiddle_size = 0;	/* total allocated size, for debugging */

fftw_twiddle *fftw_create_twiddle(int n, int r, int m)
{
     fftw_twiddle *tw;
     FFTW_COMPLEX *W;

     /* lookup for this n in the twiddle list */
     for (tw = twlist; tw; tw = tw->next)
	  if (tw->n == n && tw->r == r && tw->m == m) {
	       ++tw->refcnt;
	       return tw;
	  }
     /* not found --- allocate a new struct twiddle */
     tw = (fftw_twiddle *) fftw_malloc(sizeof(fftw_twiddle));
     W = (FFTW_COMPLEX *) fftw_malloc(m * (r - 1) * sizeof(FFTW_COMPLEX));
     fftw_twiddle_size += n;

     tw->n = n;
     tw->r = r;
     tw->m = m;
     tw->twarray = W;
     tw->refcnt = 1;
     fftw_compute_twiddle(n, r, m, W);

     /* enqueue the new struct */
     tw->next = twlist;
     twlist = tw;

     return tw;
}

void fftw_destroy_twiddle(fftw_twiddle * tw)
{
     fftw_twiddle **p;
     --tw->refcnt;

     if (tw->refcnt == 0) {
	  /* remove from the list of known twiddle factors */
	  for (p = &twlist; p; p = &((*p)->next))
	       if (*p == tw) {
		    *p = tw->next;
		    fftw_twiddle_size -= tw->n;
		    fftw_free(tw->twarray);
		    fftw_free(tw);
		    return;
	       }
	  fftw_die("BUG in fftw_destroy_twiddle\n");
     }
}
