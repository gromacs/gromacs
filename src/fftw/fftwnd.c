
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

/* $Id$ */

#include <stdlib.h>

#include <fftw.h>

/* Prototypes for functions used internally in this file: */

static void fftw2d_out_of_place_aux(fftwnd_plan p, int howmany,
				FFTW_COMPLEX *in, int istride, int idist,
			      FFTW_COMPLEX *out, int ostride, int odist);
static void fftw3d_out_of_place_aux(fftwnd_plan p, int howmany,
				FFTW_COMPLEX *in, int istride, int idist,
			      FFTW_COMPLEX *out, int ostride, int odist);
static void fftwnd_out_of_place_aux(fftwnd_plan p, int howmany,
				FFTW_COMPLEX *in, int istride, int idist,
			      FFTW_COMPLEX *out, int ostride, int odist);

static void fftw2d_in_place_aux(fftwnd_plan p, int howmany,
			   FFTW_COMPLEX *in_out, int istride, int idist);
static void fftw3d_in_place_aux(fftwnd_plan p, int howmany,
			   FFTW_COMPLEX *in_out, int istride, int idist);
static void fftwnd_in_place_aux(fftwnd_plan p, int howmany,
			   FFTW_COMPLEX *in_out, int istride, int idist);

/*********** Initializing the FFTWND Auxiliary Data **********/

fftwnd_plan fftw2d_create_plan(int nx, int ny, fftw_direction dir, int flags)
{
     int n[2];

     n[0] = nx;
     n[1] = ny;

     return fftwnd_create_plan(2, n, dir, flags);
}

fftwnd_plan fftw3d_create_plan(int nx, int ny, int nz, fftw_direction dir,
			       int flags)
{
     int n[3];

     n[0] = nx;
     n[1] = ny;
     n[2] = nz;

     return fftwnd_create_plan(3, n, dir, flags);
}

fftwnd_plan fftwnd_create_plan(int rank, const int *n, 
			       fftw_direction dir, int flags)
{
     int i, j, max_dim = 0;
     fftwnd_plan p;
     int cur_flags;

     if (rank < 0)
	  return 0;

     for (i = 0; i < rank; ++i)
	  if (n[i] <= 0)
	       return 0;

     p = (fftwnd_plan) fftw_malloc(sizeof(fftwnd_aux_data));
     p->n = 0;
     p->n_before = 0;
     p->n_after = 0;
     p->plans = 0;
     p->work = 0;

     p->rank = rank;
     p->is_in_place = flags & FFTW_IN_PLACE;

     if (rank == 0)
	  return 0;

     p->n = (int *) fftw_malloc(sizeof(int) * rank);
     p->n_before = (int *) fftw_malloc(sizeof(int) * rank);
     p->n_after = (int *) fftw_malloc(sizeof(int) * rank);
     p->plans = (fftw_plan *) fftw_malloc(rank * sizeof(fftw_plan));
     p->n_before[0] = 1;
     p->n_after[rank - 1] = 1;

     for (i = 0; i < rank; ++i) {
	  p->n[i] = n[i];

	  if (i) {
	       p->n_before[i] = p->n_before[i - 1] * n[i - 1];
	       p->n_after[rank - 1 - i] = p->n_after[rank - i] * n[rank - i];
	  }
	  if (i < rank - 1 || (flags & FFTW_IN_PLACE)) {
	       /* fft's except the last dimension are always in-place */
	       cur_flags = flags | FFTW_IN_PLACE;
	       for (j = i - 1; j >= 0 && n[i] != n[j]; --j);

	       if (n[i] > max_dim)
		    max_dim = n[i];
	  } else {
	       cur_flags = flags;
	       /* we must create a separate plan for the last dimension */
	       j = -1;
	  }

	  if (j >= 0) {
	       /* 
	        * If a plan already exists for this size
	        * array, reuse it: 
	        */
	       p->plans[i] = p->plans[j];
	  } else {
	       /* generate a new plan: */
	       p->plans[i] = fftw_create_plan(n[i], dir, cur_flags);
	       if (!p->plans[i]) {
		    fftwnd_destroy_plan(p);
		    return 0;
	       }
	  }
     }

     /* Create work array for in-place FFTs: */
     if (max_dim > 0)
	  p->work = (FFTW_COMPLEX *)
	      fftw_malloc(sizeof(FFTW_COMPLEX) * max_dim);

     return p;
}

/************* Freeing the FFTWND Auxiliary Data *************/

void fftwnd_destroy_plan(fftwnd_plan plan)
{
     if (plan) {
	  if (plan->plans) {
	       int i, j;

	       for (i = 0; i < plan->rank; ++i) {
		    for (j = i - 1;
			 j >= 0 && plan->plans[i] != plan->plans[j];
			 --j);
		    if (j < 0 && plan->plans[i])
			 fftw_destroy_plan(plan->plans[i]);
	       }
	       fftw_free(plan->plans);
	  }
	  if (plan->n)
	       fftw_free(plan->n);

	  if (plan->n_before)
	       fftw_free(plan->n_before);

	  if (plan->n_after)
	       fftw_free(plan->n_after);

	  if (plan->work)
	       fftw_free(plan->work);

	  fftw_free(plan);
     }
}

/************** Computing the N-Dimensional FFT **************/

void fftwnd(fftwnd_plan plan, int howmany,
	    FFTW_COMPLEX *in, int istride, int idist,
	    FFTW_COMPLEX *out, int ostride, int odist)
{
     if (plan->is_in_place)	/* fft is in-place */
	  switch (plan->rank) {
	      case 0:
		   break;
	      case 1:
		   fftw(plan->plans[0], howmany, in, istride, idist,
			plan->work, 1, 0);
		   break;
	      case 2:
		   fftw2d_in_place_aux(plan, howmany, in, istride, idist);
		   break;
	      case 3:
		   fftw3d_in_place_aux(plan, howmany, in, istride, idist);
		   break;
	      default:
		   fftwnd_in_place_aux(plan, howmany, in, istride, idist);
     } else {
	  if (in == out || out == 0)
	       fftw_die("Illegal attempt to perform in-place FFT!\n");
	  switch (plan->rank) {
	      case 0:
		   break;
	      case 1:
		   fftw(plan->plans[0], howmany, in, istride, idist,
			out, ostride, odist);
		   break;
	      case 2:
		   fftw2d_out_of_place_aux(plan, howmany, in, istride,
					   idist, out, ostride, odist);
		   break;
	      case 3:
		   fftw3d_out_of_place_aux(plan, howmany, in, istride,
					   idist, out, ostride, odist);
		   break;
	      default:
		   fftwnd_out_of_place_aux(plan, howmany, in, istride,
					   idist, out, ostride, odist);
	  }
     }
}

static void fftw2d_out_of_place_aux(fftwnd_plan p, int howmany,
				FFTW_COMPLEX *in, int istride, int idist,
			       FFTW_COMPLEX *out, int ostride, int odist)
{
     int fft_iter;
     fftw_plan p0, p1;
     int n0, n1;

     p0 = p->plans[0];
     p1 = p->plans[1];
     n0 = p->n[0];
     n1 = p->n[1];

     for (fft_iter = 0; fft_iter < howmany; ++fft_iter) {
	  /* FFT y dimension (out-of-place): */
	  fftw(p1, n0,
	       in + fft_iter * idist, istride, n1 * istride,
	       out + fft_iter * odist, ostride, n1 * ostride);
	  /* FFT x dimension (in-place): */
	  fftw(p0, n1,
	       out + fft_iter * odist, n1 * ostride, ostride,
	       p->work, 1, 1);
     }
}

static void fftw3d_out_of_place_aux(fftwnd_plan p, int howmany,
				FFTW_COMPLEX *in, int istride, int idist,
			       FFTW_COMPLEX *out, int ostride, int odist)
{
     int fft_iter;
     int i;
     fftw_plan p0, p1, p2;
     int n0, n1, n2;

     p0 = p->plans[0];
     p1 = p->plans[1];
     p2 = p->plans[2];
     n0 = p->n[0];
     n1 = p->n[1];
     n2 = p->n[2];

     for (fft_iter = 0; fft_iter < howmany; ++fft_iter) {
	  /* FFT z dimension (out-of-place): */
	  fftw(p2, n0 * n1,
	       in + fft_iter * idist, istride, n2 * istride,
	       out + fft_iter * odist, ostride, n2 * ostride);
	  /* FFT y dimension (in-place): */
	  for (i = 0; i < n0; ++i)
	       fftw(p1, n2,
		    out + fft_iter * odist + i * n1 * n2 * ostride,
		    n2 * ostride, ostride, p->work, 1, 0);
	  /* FFT x dimension (in-place): */
	  fftw(p0, n1 * n2,
	       out + fft_iter * odist, n1 * n2 * ostride, ostride,
	       p->work, 1, 0);
     }
}

static void fftwnd_out_of_place_aux(fftwnd_plan p, int howmany,
				FFTW_COMPLEX *in, int istride, int idist,
			       FFTW_COMPLEX *out, int ostride, int odist)
{
     int fft_iter;
     int j, i;

     /* Do FFT for rank > 3: */

     for (fft_iter = 0; fft_iter < howmany; ++fft_iter) {
	  /* do last dimension (out-of-place): */
	  fftw(p->plans[p->rank - 1], p->n_before[p->rank - 1],
	     in + fft_iter * idist, istride, p->n[p->rank - 1] * istride,
	   out + fft_iter * odist, ostride, p->n[p->rank - 1] * ostride);

	  /* do first dimension (in-place): */
	  fftw(p->plans[0], p->n_after[0],
	       out + fft_iter * odist, p->n_after[0] * ostride, ostride,
	       p->work, 1, 0);

	  /* do other dimensions (in-place): */
	  for (j = 1; j < p->rank - 1; ++j)
	       for (i = 0; i < p->n_before[j]; ++i)
		    fftw(p->plans[j], p->n_after[j],
			 out + fft_iter * odist + i * ostride * p->n[j] *
			 p->n_after[j], p->n_after[j] * ostride,
			 ostride, p->work, 1, 0);
     }
}

static void fftw2d_in_place_aux(fftwnd_plan p, int howmany,
			    FFTW_COMPLEX *in_out, int istride, int idist)
{
     int fft_iter;
     fftw_plan p0, p1;
     int n0, n1;

     p0 = p->plans[0];
     p1 = p->plans[1];
     n0 = p->n[0];
     n1 = p->n[1];

     for (fft_iter = 0; fft_iter < howmany; ++fft_iter) {
	  /* FFT y dimension: */
	  fftw(p1, n0,
	       in_out + fft_iter * idist, istride, istride * n1,
	       p->work, 1, 0);
	  /* FFT x dimension: */
	  fftw(p0, n1,
	       in_out + fft_iter * idist, istride * n1, istride,
	       p->work, 1, 0);
     }
}

static void fftw3d_in_place_aux(fftwnd_plan p, int howmany,
			    FFTW_COMPLEX *in_out, int istride, int idist)
{
     int i;
     int fft_iter;
     fftw_plan p0, p1, p2;
     int n0, n1, n2;

     p0 = p->plans[0];
     p1 = p->plans[1];
     p2 = p->plans[2];
     n0 = p->n[0];
     n1 = p->n[1];
     n2 = p->n[2];

     for (fft_iter = 0; fft_iter < howmany; ++fft_iter) {
	  /* FFT z dimension: */
	  fftw(p2, n0 * n1,
	       in_out + fft_iter * idist, istride, n2 * istride,
	       p->work, 1, 0);
	  /* FFT y dimension: */
	  for (i = 0; i < n0; ++i)
	       fftw(p1, n2,
		    in_out + fft_iter * idist + i * n1 *
		    n2 * istride, n2 * istride, istride, p->work, 1, 0);
	  /* FFT x dimension: */
	  fftw(p0, n1 * n2,
	       in_out + fft_iter * idist, n1 * n2 * istride, istride,
	       p->work, 1, 0);
     }
}

static void fftwnd_in_place_aux(fftwnd_plan p, int howmany,
			    FFTW_COMPLEX *in_out, int istride, int idist)
/* Do FFT for rank > 3: */
{
     int fft_iter;
     int j, i;

     for (fft_iter = 0; fft_iter < howmany; ++fft_iter) {
	  /* do last dimension: */
	  fftw(p->plans[p->rank - 1], p->n_before[p->rank - 1],
	  in_out + fft_iter * idist, istride, p->n[p->rank - 1] * istride,
	       p->work, 1, 0);

	  /* do first dimension: */
	  fftw(p->plans[0], p->n_after[0],
	     in_out + fft_iter * idist, p->n_after[0] * istride, istride,
	       p->work, 1, 0);

	  /* do other dimensions: */
	  for (j = 1; j < p->rank - 1; ++j)
	       for (i = 0; i < p->n_before[j]; ++i)
		    fftw(p->plans[j], p->n_after[j],
		      in_out + fft_iter * idist + i * istride * p->n[j] *
			 p->n_after[j], p->n_after[j] * istride, istride,
			 p->work, 1, 0);
     }
}
