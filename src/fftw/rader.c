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
 * Compute transforms of prime sizes using Rader's trick: turn them
 * into convolutions of size n - 1, which you then perform via a pair
 * of FFTs. 
 */

#include <stdlib.h>
#include <math.h>

#include <fftw-int.h>

#ifdef FFTW_USING_CILK
#include <cilk.h>
#include <cilk-compat.h>
#endif

#ifdef FFTW_DEBUG
#define WHEN_DEBUG(a) a
#else
#define WHEN_DEBUG(a)
#endif

/* compute n^m mod p, where m >= 0 and p > 0. */
static int power_mod(int n, int m, int p)
{
     if (m == 0)
	  return 1;
     else if (m % 2 == 0) {
	  int x = power_mod(n, m / 2, p);
	  return ((x * x) % p);
     } else
	  return ((n * power_mod(n, m - 1, p)) % p);
}

/*
 * Find the period of n in the multiplicative group mod p (p prime).
 * That is, return the smallest m such that n^m == 1 mod p.
 */
static int period(int n, int p)
{
     int prod = n, period = 1;

     while (prod != 1) {
	  prod = (prod * n) % p;
	  ++period;
	  if (prod == 0)
	       fftw_die("non-prime order in Rader\n");
     }
     return period;
}

/* find a generator for the multiplicative group mod p, where p is prime */
static int find_generator(int p)
{
     int g;

     for (g = 1; g < p; ++g)
	  if (period(g, p) == p - 1)
	       break;
     if (g == p)
	  fftw_die("couldn't find generator for Rader\n");
     return g;
}

/***************************************************************************/

static fftw_rader_data *create_rader_aux(int p, int flags)
{
     fftw_complex *omega, *work;
     int g, ginv, gpower;
     int i;
     FFTW_TRIG_REAL twoPiOverN;
     fftw_real scale = 1.0 / (p - 1);	/* for convolution */
     fftw_plan plan;
     fftw_rader_data *d;

     if (p < 2)
	  fftw_die("non-prime order in Rader\n");

     flags &= ~FFTW_IN_PLACE;

     d = (fftw_rader_data *) fftw_malloc(sizeof(fftw_rader_data));

     g = find_generator(p);
     ginv = power_mod(g, p - 2, p);

     omega = (fftw_complex *) fftw_malloc((p - 1) * sizeof(fftw_complex));

     plan = fftw_create_plan(p - 1, FFTW_FORWARD, flags);

     work = (fftw_complex *) fftw_malloc((p - 1) * sizeof(fftw_complex));

     twoPiOverN = FFTW_K2PI / (FFTW_TRIG_REAL) p;
     gpower = 1;
     for (i = 0; i < p - 1; ++i) {
	  c_re(work[i]) = scale * FFTW_TRIG_COS(twoPiOverN * gpower);
	  c_im(work[i]) = FFTW_FORWARD * scale * FFTW_TRIG_SIN(twoPiOverN * gpower);
	  gpower = (gpower * ginv) % p;
     }

     /* fft permuted roots of unity */
     fftw_executor_simple(p - 1, work, omega, plan->root, 1, 1);

     fftw_free(work);

     d->plan = plan;
     d->omega = omega;
     d->g = g;
     d->ginv = ginv;
     d->p = p;
     d->flags = flags;
     d->refcount = 1;
     d->next = NULL;

     d->cdesc = (fftw_codelet_desc *) fftw_malloc(sizeof(fftw_codelet_desc));
     d->cdesc->name = NULL;
     d->cdesc->codelet = NULL;
     d->cdesc->size = p;
     d->cdesc->dir = FFTW_FORWARD;
     d->cdesc->type = FFTW_RADER;
     d->cdesc->signature = g;
     d->cdesc->ntwiddle = 0;
     d->cdesc->twiddle_order = NULL;
     return d;
}

/***************************************************************************/

static fftw_rader_data *fftw_create_rader(int p, int flags)
{
     fftw_rader_data *d = fftw_rader_top;

     flags &= ~FFTW_IN_PLACE;
     while (d && (d->p != p || d->flags != flags))
	  d = d->next;
     if (d) {
	  d->refcount++;
	  return d;
     }
     d = create_rader_aux(p, flags);
     d->next = fftw_rader_top;
     fftw_rader_top = d;
     return d;
}

/***************************************************************************/

/* Compute the prime FFTs, premultiplied by twiddle factors.  Below, we
 * extensively use the identity that fft(x*)* = ifft(x) in order to
 * share data between forward and backward transforms and to obviate
 * the necessity of having separate forward and backward plans. */

void fftw_twiddle_rader(fftw_complex *A, const fftw_complex *W,
			int m, int r, int stride,
			fftw_rader_data * d)
{
     fftw_complex *tmp = (fftw_complex *)
     fftw_malloc((r - 1) * sizeof(fftw_complex));
     int i, k, gpower = 1, g = d->g, ginv = d->ginv;
     fftw_real a0r, a0i;
     fftw_complex *omega = d->omega;

     for (i = 0; i < m; ++i, A += stride, W += r - 1) {
	  /* 
	   * Here, we fft W[k-1] * A[k*(m*stride)], using Rader.
	   * (Actually, W is pre-permuted to match the permutation that we 
	   * will do on A.) 
	   */

	  /* First, permute the input and multiply by W, storing in tmp: */
	  /* gpower == g^k mod r in the following loop */
	  for (k = 0; k < r - 1; ++k, gpower = (gpower * g) % r) {
	       fftw_real rA, iA, rW, iW;
	       rW = c_re(W[k]);
	       iW = c_im(W[k]);
	       rA = c_re(A[gpower * (m * stride)]);
	       iA = c_im(A[gpower * (m * stride)]);
	       c_re(tmp[k]) = rW * rA - iW * iA;
	       c_im(tmp[k]) = rW * iA + iW * rA;
	  }

	  WHEN_DEBUG( {
		     if (gpower != 1)
		     fftw_die("incorrect generator in Rader\n");
		     }
	  );

	  /* FFT tmp to A: */
	  fftw_executor_simple(r - 1, tmp, A + (m * stride),
			       d->plan->root, 1, m * stride);

	  /* set output DC component: */
	  a0r = c_re(A[0]);
	  a0i = c_im(A[0]);
	  c_re(A[0]) += c_re(A[(m * stride)]);
	  c_im(A[0]) += c_im(A[(m * stride)]);

	  /* now, multiply by omega: */
	  for (k = 0; k < r - 1; ++k) {
	       fftw_real rA, iA, rW, iW;
	       rW = c_re(omega[k]);
	       iW = c_im(omega[k]);
	       rA = c_re(A[(k + 1) * (m * stride)]);
	       iA = c_im(A[(k + 1) * (m * stride)]);
	       c_re(A[(k + 1) * (m * stride)]) = rW * rA - iW * iA;
	       c_im(A[(k + 1) * (m * stride)]) = -(rW * iA + iW * rA);
	  }

	  /* this will add A[0] to all of the outputs after the ifft */
	  c_re(A[(m * stride)]) += a0r;
	  c_im(A[(m * stride)]) -= a0i;

	  /* inverse FFT: */
	  fftw_executor_simple(r - 1, A + (m * stride), tmp,
			       d->plan->root, m * stride, 1);

	  /* finally, do inverse permutation to unshuffle the output: */
	  for (k = 0; k < r - 1; ++k, gpower = (gpower * ginv) % r) {
	       c_re(A[gpower * (m * stride)]) = c_re(tmp[k]);
	       c_im(A[gpower * (m * stride)]) = -c_im(tmp[k]);
	  }

	  WHEN_DEBUG( {
		     if (gpower != 1)
		     fftw_die("incorrect generator in Rader\n");
		     }
	  );

     }

     fftw_free(tmp);
}

void fftwi_twiddle_rader(fftw_complex *A, const fftw_complex *W,
			 int m, int r, int stride,
			 fftw_rader_data * d)
{
     fftw_complex *tmp = (fftw_complex *)
     fftw_malloc((r - 1) * sizeof(fftw_complex));
     int i, k, gpower = 1, g = d->g, ginv = d->ginv;
     fftw_real a0r, a0i;
     fftw_complex *omega = d->omega;

     for (i = 0; i < m; ++i, A += stride, W += r - 1) {
	  /* 
	   * Here, we fft W[k-1]* * A[k*(m*stride)], using Rader. 
	   * (Actually, W is pre-permuted to match the permutation that
	   * we will do on A.) 
	   */

	  /* First, permute the input and multiply by W*, storing in tmp: */
	  /* gpower == g^k mod r in the following loop */
	  for (k = 0; k < r - 1; ++k, gpower = (gpower * g) % r) {
	       fftw_real rA, iA, rW, iW;
	       rW = c_re(W[k]);
	       iW = c_im(W[k]);
	       rA = c_re(A[gpower * (m * stride)]);
	       iA = c_im(A[gpower * (m * stride)]);
	       c_re(tmp[k]) = rW * rA + iW * iA;
	       c_im(tmp[k]) = iW * rA - rW * iA;
	  }

	  WHEN_DEBUG( {
		     if (gpower != 1)
		     fftw_die("incorrect generator in Rader\n");
		     }
	  );

	  /* FFT tmp to A: */
	  fftw_executor_simple(r - 1, tmp, A + (m * stride),
			       d->plan->root, 1, m * stride);

	  /* set output DC component: */
	  a0r = c_re(A[0]);
	  a0i = c_im(A[0]);
	  c_re(A[0]) += c_re(A[(m * stride)]);
	  c_im(A[0]) -= c_im(A[(m * stride)]);

	  /* now, multiply by omega: */
	  for (k = 0; k < r - 1; ++k) {
	       fftw_real rA, iA, rW, iW;
	       rW = c_re(omega[k]);
	       iW = c_im(omega[k]);
	       rA = c_re(A[(k + 1) * (m * stride)]);
	       iA = c_im(A[(k + 1) * (m * stride)]);
	       c_re(A[(k + 1) * (m * stride)]) = rW * rA - iW * iA;
	       c_im(A[(k + 1) * (m * stride)]) = -(rW * iA + iW * rA);
	  }

	  /* this will add A[0] to all of the outputs after the ifft */
	  c_re(A[(m * stride)]) += a0r;
	  c_im(A[(m * stride)]) += a0i;

	  /* inverse FFT: */
	  fftw_executor_simple(r - 1, A + (m * stride), tmp,
			       d->plan->root, m * stride, 1);

	  /* finally, do inverse permutation to unshuffle the output: */
	  for (k = 0; k < r - 1; ++k, gpower = (gpower * ginv) % r) {
	       A[gpower * (m * stride)] = tmp[k];
	  }

	  WHEN_DEBUG( {
		     if (gpower != 1)
		     fftw_die("incorrect generator in Rader\n");
		     }
	  );
     }

     fftw_free(tmp);
}

/***************************************************************************/

/*
 * Make an FFTW_RADER plan node.  Note that this function must go
 * here, rather than in putils.c, because it indirectly calls the
 * fftw_planner.  If we included it in putils.c, which is also used
 * by rfftw, then any program using rfftw would be linked with all
 * of the FFTW codelets, even if they were not needed.   I wish that the
 * darn linkers operated on a function rather than a file granularity. 
 */
fftw_plan_node *fftw_make_node_rader(int n, int size, fftw_direction dir,
				     fftw_plan_node *recurse,
				     int flags)
{
     fftw_plan_node *p = fftw_make_node();

     p->type = FFTW_RADER;
     p->nodeu.rader.size = size;
     p->nodeu.rader.codelet = dir == FFTW_FORWARD ?
	 fftw_twiddle_rader : fftwi_twiddle_rader;
     p->nodeu.rader.rader_data = fftw_create_rader(size, flags);
     p->nodeu.rader.recurse = recurse;
     fftw_use_node(recurse);

     if (flags & FFTW_MEASURE)
	  p->nodeu.rader.tw =
	      fftw_create_twiddle(n, p->nodeu.rader.rader_data->cdesc);
     else
	  p->nodeu.rader.tw = 0;
     return p;
}
