/*
 * Copyright (c) 1997,1998 Massachusetts Institute of Technology
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
 * twiddle.c -- compute twiddle factors
 * These are the twiddle factors for *direct* fft.  Flip sign to get
 * the inverse
 */

/* $Id$ */
#ifdef FFTW_USING_CILK
#include <cilk.h>
#include <cilk-compat.h>
#endif

#include <fftw-int.h>
#include <math.h>
#include <stdlib.h>

#ifndef TRUE
#define TRUE (1 == 1)
#endif

#ifndef FALSE
#define FALSE (1 == 0)
#endif

static fftw_complex *fftw_compute_rader_twiddle(int n, int r, int g)
{
     FFTW_TRIG_REAL twoPiOverN;
     int m = n / r;
     int i, j, gpower;
     fftw_complex *W;

     twoPiOverN = FFTW_K2PI / (FFTW_TRIG_REAL) n;
     W = (fftw_complex *) fftw_malloc((r - 1) * m * sizeof(fftw_complex));
     for (i = 0; i < m; ++i)
	  for (gpower = 1, j = 0; j < r - 1; ++j, gpower = (gpower * g) % r) {
	       int k = i * (r - 1) + j;
	       FFTW_TRIG_REAL
		   ij = (FFTW_TRIG_REAL) (i * gpower);
	       c_re(W[k]) = FFTW_TRIG_COS(twoPiOverN * ij);
	       c_im(W[k]) = FFTW_FORWARD * FFTW_TRIG_SIN(twoPiOverN * ij);
	  }

     return W;
}

/*
 * compute the W coefficients (that is, powers of the root of 1)
 * and store them into an array.
 */
static fftw_complex *fftw_compute_twiddle(int n, const fftw_codelet_desc *d)
{
     FFTW_TRIG_REAL twoPiOverN;
     int i, j;
     fftw_complex *W;

     twoPiOverN = FFTW_K2PI / (FFTW_TRIG_REAL) n;

     if (!d) {
	  /* generic codelet, needs all twiddles in order */
	  W = (fftw_complex *) fftw_malloc(n * sizeof(fftw_complex));
	  for (i = 0; i < n; ++i) {
	       c_re(W[i]) = FFTW_TRIG_COS(twoPiOverN * (FFTW_TRIG_REAL) i);
	       c_im(W[i]) = FFTW_FORWARD * FFTW_TRIG_SIN(twoPiOverN * (FFTW_TRIG_REAL) i);
	  }
     } else if (d->type == FFTW_RADER)
	  W = fftw_compute_rader_twiddle(n, d->size, d->signature);
     else {
	  int r = d->size;
	  int m = n / r, m_alloc;
	  int r1 = d->ntwiddle;
	  int istart;

	  if (d->type == FFTW_TWIDDLE) {
	       istart = 0;
	       m_alloc = m;
	  } else if (d->type == FFTW_HC2HC) {
	       m = m / 2 + 1;
	       m_alloc = m - 1;
	       istart = 1;
	  } else {
	       fftw_die("compute_twiddle: invalid argument\n");
	       /* paranoia for gcc */
	       m_alloc = 0;
	       istart = 0;
	  }

	  W = (fftw_complex *) fftw_malloc(r1 * m_alloc * sizeof(fftw_complex));
	  for (i = istart; i < m; ++i)
	       for (j = 0; j < r1; ++j) {
		    int k = (i - istart) * r1 + j;
		    FFTW_TRIG_REAL
			ij = (FFTW_TRIG_REAL) (i * d->twiddle_order[j]);
		    c_re(W[k]) = FFTW_TRIG_COS(twoPiOverN * ij);
		    c_im(W[k]) = FFTW_FORWARD * FFTW_TRIG_SIN(twoPiOverN * ij);
	       }
     }

     return W;
}

/*
 * these routines implement a simple reference-count-based 
 * management of twiddle structures
 */
static fftw_twiddle *twlist = (fftw_twiddle *) 0;
int fftw_twiddle_size = 0;	/* total allocated size, for debugging */

/* true if the two codelets can share the same twiddle factors */
static int compatible(const fftw_codelet_desc *d1, const fftw_codelet_desc *d2)
{
     int i;

     /* true if they are the same codelet */
     if (d1 == d2)
	  return TRUE;

     /* false if one is null and the other is not */
     if (!d1 || !d2)
	  return FALSE;

     /* false if size is different */
     if (d1->size != d2->size)
	  return FALSE;

     /* false if different types (FFTW_TWIDDLE/FFTW_HC2HC/FFTW_RADER) */
     if (d1->type != d2->type)
	  return FALSE;

     /* false if they need different # of twiddles */
     if (d1->ntwiddle != d2->ntwiddle)
	  return FALSE;

     /* false if the twiddle orders are different */
     for (i = 0; i < d1->ntwiddle; ++i)
	  if (d1->twiddle_order[i] != d2->twiddle_order[i])
	       return FALSE;

     return TRUE;
}

fftw_twiddle *fftw_create_twiddle(int n, const fftw_codelet_desc *d)
{
     fftw_twiddle *tw;

     /* lookup this n in the twiddle list */
     for (tw = twlist; tw; tw = tw->next)
	  if (n == tw->n && compatible(d, tw->cdesc)) {
	       ++tw->refcnt;
	       return tw;
	  }
     /* not found --- allocate a new struct twiddle */
     tw = (fftw_twiddle *) fftw_malloc(sizeof(fftw_twiddle));
     fftw_twiddle_size += n;

     tw->n = n;
     tw->cdesc = d;
     tw->twarray = fftw_compute_twiddle(n, d);
     tw->refcnt = 1;

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
