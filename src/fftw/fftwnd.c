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

/* $Id$ */

#include <fftw-int.h>

/* the number of buffers to use for buffered transforms: */
#define FFTWND_NBUFFERS 8

/* the default number of buffers to use: */
#define FFTWND_DEFAULT_NBUFFERS 0

/* the number of "padding" elements between consecutive buffer lines */
#define FFTWND_BUFFER_PADDING 8

static void destroy_plan_array(int rank, fftw_plan *plans);

static void init_test_array(fftw_complex *arr, int stride, int n)
{
     int j;

     for (j = 0; j < n; ++j) {
	  c_re(arr[stride * j]) = 0.0;
	  c_im(arr[stride * j]) = 0.0;
     }
}

/*
 * Same as fftw_measure_runtime, except for fftwnd plan.
 */
double fftwnd_measure_runtime(fftwnd_plan plan,
			      fftw_complex *in, int istride,
			      fftw_complex *out, int ostride)
{
     fftw_time begin, end, start;
     double t, tmax, tmin;
     int i, iter;
     int n;
     int repeat;

     if (plan->rank == 0)
	  return 0.0;

     n = 1;
     for (i = 0; i < plan->rank; ++i)
	  n *= plan->n[i];

     iter = 1;

     for (;;) {
	  tmin = 1.0E10;
	  tmax = -1.0E10;
	  init_test_array(in, istride, n);

	  start = fftw_get_time();
	  /* repeat the measurement FFTW_TIME_REPEAT times */
	  for (repeat = 0; repeat < FFTW_TIME_REPEAT; ++repeat) {
	       begin = fftw_get_time();
	       for (i = 0; i < iter; ++i) {
		    fftwnd(plan, 1, in, istride, 0, out, ostride, 0);
	       }
	       end = fftw_get_time();

	       t = fftw_time_to_sec(fftw_time_diff(end, begin));
	       if (t < tmin)
		    tmin = t;
	       if (t > tmax)
		    tmax = t;

	       /* do not run for too long */
	       t = fftw_time_to_sec(fftw_time_diff(end, start));
	       if (t > FFTW_TIME_LIMIT)
		    break;
	  }

	  if (tmin >= FFTW_TIME_MIN)
	       break;

	  iter *= 2;
     }

     tmin /= (double) iter;
     tmax /= (double) iter;

     return tmin;
}

/********************** Initializing the FFTWND Plan ***********************/

/* Initialize everything except for the 1D plans and the work array: */
fftwnd_plan fftwnd_create_plan_aux(int rank, const int *n,
				   fftw_direction dir, int flags)
{
     int i;
     fftwnd_plan p;

     if (rank < 0)
	  return 0;

     for (i = 0; i < rank; ++i)
	  if (n[i] <= 0)
	       return 0;

     p = (fftwnd_plan) fftw_malloc(sizeof(fftwnd_data));
     p->n = 0;
     p->n_before = 0;
     p->n_after = 0;
     p->plans = 0;
     p->work = 0;
     p->dir = dir;

     p->rank = rank;
     p->is_in_place = flags & FFTW_IN_PLACE;

     p->nwork = 0;
     p->nbuffers = 0;

     if (rank == 0)
	  return 0;

     p->n = (int *) fftw_malloc(sizeof(int) * rank);
     p->n_before = (int *) fftw_malloc(sizeof(int) * rank);
     p->n_after = (int *) fftw_malloc(sizeof(int) * rank);
     p->n_before[0] = 1;
     p->n_after[rank - 1] = 1;

     for (i = 0; i < rank; ++i) {
	  p->n[i] = n[i];

	  if (i) {
	       p->n_before[i] = p->n_before[i - 1] * n[i - 1];
	       p->n_after[rank - 1 - i] = p->n_after[rank - i] * n[rank - i];
	  }
     }

     return p;
}

/* create an empty new array of rank 1d plans */
fftw_plan *fftwnd_new_plan_array(int rank)
{
     fftw_plan *plans;
     int i;

     plans = (fftw_plan *) fftw_malloc(rank * sizeof(fftw_plan));
     if (!plans)
	  return 0;
     for (i = 0; i < rank; ++i)
	  plans[i] = 0;
     return plans;
}

/* 
 * create an array of plans using the ordinary 1d fftw_create_plan,
 * which allocates its own array and creates plans optimized for
 * contiguous data. 
 */
fftw_plan *fftwnd_create_plans_generic(fftw_plan *plans,
				       int rank, const int *n,
				       fftw_direction dir, int flags)
{
     if (rank <= 0)
	  return 0;

     if (plans) {
	  int i, j;
	  int cur_flags;

	  for (i = 0; i < rank; ++i) {
	       if (i < rank - 1 || (flags & FFTW_IN_PLACE)) {
		    /* 
		     * fft's except the last dimension are always in-place 
		     */
		    cur_flags = flags | FFTW_IN_PLACE;
		    for (j = i - 1; j >= 0 && n[i] != n[j]; --j);
	       } else {
		    cur_flags = flags;
		    /* 
		     * we must create a separate plan for the last
		     * dimension 
		     */
		    j = -1;
	       }

	       if (j >= 0) {
		    /* 
		     * If a plan already exists for this size
		     * array, reuse it: 
		     */
		    plans[i] = plans[j];
	       } else {
		    /* generate a new plan: */
		    plans[i] = fftw_create_plan(n[i], dir, cur_flags);
		    if (!plans[i]) {
			 destroy_plan_array(rank, plans);
			 return 0;
		    }
	       }
	  }
     }
     return plans;
}

static int get_maxdim(int rank, const int *n, int flags)
{
     int i;
     int maxdim = 0;

     for (i = 0; i < rank - 1; ++i)
	  if (n[i] > maxdim)
	       maxdim = n[i];
     if (rank > 0 && flags & FFTW_IN_PLACE && n[rank - 1] > maxdim)
	  maxdim = n[rank - 1];

     return maxdim;
}

/* compute number of elements required for work array (has to
   be big enough to hold ncopies of the largest dimension in
   n that will need an in-place transform. */
int fftwnd_work_size(int rank, const int *n, int flags, int ncopies)
{
     return (ncopies * get_maxdim(rank, n, flags)
	     + (ncopies - 1) * FFTWND_BUFFER_PADDING);
}

/*
 * create plans using the fftw_create_plan_specific planner, which
 * allows us to create plans for each dimension that are specialized
 * for the strides that we are going to use. 
 */
fftw_plan *fftwnd_create_plans_specific(fftw_plan *plans,
					int rank, const int *n,
					const int *n_after,
					fftw_direction dir, int flags,
					fftw_complex *in, int istride,
					fftw_complex *out, int ostride)
{
     if (rank <= 0)
	  return 0;

     if (plans) {
	  int i, stride, cur_flags;
	  fftw_complex *work = 0;
	  int nwork;

	  nwork = fftwnd_work_size(rank, n, flags, 1);
	  if (nwork)
	       work = (fftw_complex*)fftw_malloc(nwork * sizeof(fftw_complex));

	  for (i = 0; i < rank; ++i) {
	       /* fft's except the last dimension are always in-place */
	       if (i < rank - 1)
		    cur_flags = flags | FFTW_IN_PLACE;
	       else
		    cur_flags = flags;

	       /* stride for transforming ith dimension */
	       stride = n_after[i];

	       if (cur_flags & FFTW_IN_PLACE)
		    plans[i] = fftw_create_plan_specific(n[i], dir, cur_flags,
						    in, istride * stride,
							 work, 1);
	       else
		    plans[i] = fftw_create_plan_specific(n[i], dir, cur_flags,
						    in, istride * stride,
						  out, ostride * stride);
	       if (!plans[i]) {
		    destroy_plan_array(rank, plans);
		    fftw_free(work);
		    return 0;
	       }
	  }

	  if (work)
	       fftw_free(work);
     }
     return plans;
}

/*
 * Create an fftwnd_plan specialized for specific arrays.  (These
 * arrays are ignored, however, if they are NULL or if the flags do
 * not include FFTW_MEASURE.)  The main advantage of being provided
 * arrays like this is that we can do runtime timing measurements of
 * our options, without worrying about allocating excessive scratch
 * space.
 */
fftwnd_plan fftwnd_create_plan_specific(int rank, const int *n,
					fftw_direction dir, int flags,
					fftw_complex *in, int istride,
					fftw_complex *out, int ostride)
{
     fftwnd_plan p;

     if (!(p = fftwnd_create_plan_aux(rank, n, dir, flags)))
	  return 0;

     if (!(flags & FFTW_MEASURE) || in == 0
	 || (!p->is_in_place && out == 0)) {

/**** use default plan ****/

	  p->plans = fftwnd_create_plans_generic(fftwnd_new_plan_array(rank),
						 rank, n, dir, flags);
	  if (!p->plans) {
	       fftwnd_destroy_plan(p);
	       return 0;
	  }
	  if (flags & FFTWND_FORCE_BUFFERED)
	       p->nbuffers = FFTWND_NBUFFERS;
	  else
	       p->nbuffers = FFTWND_DEFAULT_NBUFFERS;

	  p->nwork = fftwnd_work_size(rank, n, flags, p->nbuffers + 1);
	  if (p->nwork && !(flags & FFTW_THREADSAFE)) {
	       p->work = (fftw_complex*) fftw_malloc(p->nwork 
						     * sizeof(fftw_complex));
	       if (!p->work) {
		    fftwnd_destroy_plan(p);
		    return 0;
	       }
	  }
     } else {
/**** use runtime measurements to pick plan ****/

	  fftw_plan *plans_buf, *plans_nobuf;
	  double t_buf, t_nobuf;

	  p->nwork = fftwnd_work_size(rank, n, flags, FFTWND_NBUFFERS + 1);
	  if (p->nwork && !(flags & FFTW_THREADSAFE)) {
	       p->work = (fftw_complex*) fftw_malloc(p->nwork 
						     * sizeof(fftw_complex));
	       if (!p->work) {
		    fftwnd_destroy_plan(p);
		    return 0;
	       }
	  }
	  else
	       p->work = (fftw_complex*) NULL;

	  /* two possible sets of 1D plans: */
	  plans_buf = fftwnd_create_plans_generic(fftwnd_new_plan_array(rank),
						  rank, n, dir, flags);
	  plans_nobuf = 
	       fftwnd_create_plans_specific(fftwnd_new_plan_array(rank),
					    rank, n, p->n_after, dir,
					    flags, in, istride,
					    out, ostride);
	  if (!plans_buf || !plans_nobuf) {
	       destroy_plan_array(rank, plans_nobuf);
	       destroy_plan_array(rank, plans_buf);
	       fftwnd_destroy_plan(p);
	       return 0;
	  }
	  /* time the two possible plans */
	  p->plans = plans_nobuf;
	  p->nbuffers = 0;
	  p->nwork = fftwnd_work_size(rank, n, flags, p->nbuffers + 1);
	  t_nobuf = fftwnd_measure_runtime(p, in, istride, out, ostride);
	  p->plans = plans_buf;
	  p->nbuffers = FFTWND_NBUFFERS;
	  p->nwork = fftwnd_work_size(rank, n, flags, p->nbuffers + 1);
	  t_buf = fftwnd_measure_runtime(p, in, istride, out, ostride);

	  /* pick the better one: */
	  if (t_nobuf < t_buf) {	/* use unbuffered transform */
	       p->plans = plans_nobuf;
	       p->nbuffers = 0;

	       /* work array is unnecessarily large */
	       if (p->work)
		    fftw_free(p->work);
	       p->work = 0;

	       destroy_plan_array(rank, plans_buf);

	       /* allocate a work array of the correct size: */
	       p->nwork = fftwnd_work_size(rank, n, flags, p->nbuffers + 1);
	       if (p->nwork && !(flags & FFTW_THREADSAFE)) {
		    p->work = (fftw_complex*) fftw_malloc(p->nwork 
						       * sizeof(fftw_complex));
		    if (!p->work) {
			 fftwnd_destroy_plan(p);
			 return 0;
		    }
	       }
	  } else {		/* use buffered transform */
	       destroy_plan_array(rank, plans_nobuf);
	  }
     }

     return p;
}

fftwnd_plan fftw2d_create_plan_specific(int nx, int ny,
					fftw_direction dir, int flags,
					fftw_complex *in, int istride,
					fftw_complex *out, int ostride)
{
     int n[2];

     n[0] = nx;
     n[1] = ny;

     return fftwnd_create_plan_specific(2, n, dir, flags,
					in, istride, out, ostride);
}

fftwnd_plan fftw3d_create_plan_specific(int nx, int ny, int nz,
					fftw_direction dir, int flags,
					fftw_complex *in, int istride,
					fftw_complex *out, int ostride)
{
     int n[3];

     n[0] = nx;
     n[1] = ny;
     n[2] = nz;

     return fftwnd_create_plan_specific(3, n, dir, flags,
					in, istride, out, ostride);
}

/* Create a generic fftwnd plan: */

fftwnd_plan fftwnd_create_plan(int rank, const int *n,
			       fftw_direction dir, int flags)
{
     return fftwnd_create_plan_specific(rank, n, dir, flags, 0, 1, 0, 1);
}

fftwnd_plan fftw2d_create_plan(int nx, int ny,
			       fftw_direction dir, int flags)
{
     return fftw2d_create_plan_specific(nx, ny, dir, flags, 0, 1, 0, 1);
}

fftwnd_plan fftw3d_create_plan(int nx, int ny, int nz,
			       fftw_direction dir, int flags)
{
     return fftw3d_create_plan_specific(nx, ny, nz, dir, flags, 0, 1, 0, 1);
}

/************************ Freeing the FFTWND Plan ************************/

static void destroy_plan_array(int rank, fftw_plan *plans)
{
     if (plans) {
	  int i, j;

	  for (i = 0; i < rank; ++i) {
	       for (j = i - 1;
		    j >= 0 && plans[i] != plans[j];
		    --j);
	       if (j < 0 && plans[i])
		    fftw_destroy_plan(plans[i]);
	  }
	  fftw_free(plans);
     }
}

void fftwnd_destroy_plan(fftwnd_plan plan)
{
     if (plan) {
	  destroy_plan_array(plan->rank, plan->plans);

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

/************************ Printing the FFTWND Plan ************************/

void fftwnd_fprint_plan(FILE *f, fftwnd_plan plan)
{
     if (plan) {
	  int i, j;

	  if (plan->rank == 0) {
	       fprintf(f, "plan for rank 0 (null) transform.\n");
	       return;
	  }
	  fprintf(f, "plan for ");
	  for (i = 0; i < plan->rank; ++i)
	       fprintf(f, "%s%d", i ? "x" : "", plan->n[i]);
	  fprintf(f, " transform:\n");

	  if (plan->nbuffers > 0)
	       fprintf(f, "  -- using buffered transforms (%d buffers)\n",
		       plan->nbuffers);
	  else
	       fprintf(f, "  -- using unbuffered transform\n");

	  for (i = 0; i < plan->rank; ++i) {
	       fprintf(f, "* dimension %d (size %d) ", i, plan->n[i]);

	       for (j = i - 1; j >= 0; --j)
		    if (plan->plans[j] == plan->plans[i])
			 break;

	       if (j < 0)
		    fftw_fprint_plan(f, plan->plans[i]);
	       else
		    fprintf(f, "plan is same as dimension %d plan.\n", j);
	  }
     }
}

void fftwnd_print_plan(fftwnd_plan plan)
{
     fftwnd_fprint_plan(stdout, plan);
}

/********************* Buffered FFTW (in-place) *********************/

void fftw_buffered(fftw_plan p, int howmany,
		   fftw_complex *in, int istride, int idist,
		   fftw_complex *work,
		   int nbuffers, fftw_complex *buffers)
{
     int i = 0, n, nb;

     n = p->n;
     nb = n + FFTWND_BUFFER_PADDING;

     do {
	  for (; i <= howmany - nbuffers; i += nbuffers) {
	       fftw_complex *cur_in = in + i * idist;
	       int j, buf;

	       /* 
	        * First, copy nbuffers strided arrays to the
	        * contiguous buffer arrays (reading consecutive
	        * locations, assuming that idist is 1):
	        */
	       for (j = 0; j < n; ++j) {
		    fftw_complex *cur_in2 = cur_in + j * istride;
		    fftw_complex *cur_buffers = buffers + j;

		    for (buf = 0; buf <= nbuffers - 4; buf += 4) {
			 *cur_buffers = *cur_in2;
			 *(cur_buffers += nb) = *(cur_in2 += idist);
			 *(cur_buffers += nb) = *(cur_in2 += idist);
			 *(cur_buffers += nb) = *(cur_in2 += idist);
			 cur_buffers += nb;
			 cur_in2 += idist;
		    }
		    for (; buf < nbuffers; ++buf) {
			 *cur_buffers = *cur_in2;
			 cur_buffers += nb;
			 cur_in2 += idist;
		    }
	       }

	       /* 
	        * Now, compute the FFTs in the buffers (in-place
	        * using work): 
	        */
	       fftw(p, nbuffers, buffers, 1, nb, work, 1, 0);

	       /* 
	        * Finally, copy the results back from the contiguous
	        * buffers to the strided arrays (writing consecutive
	        * locations):
	        */
	       for (j = 0; j < n; ++j) {
		    fftw_complex *cur_in2 = cur_in + j * istride;
		    fftw_complex *cur_buffers = buffers + j;

		    for (buf = 0; buf <= nbuffers - 4; buf += 4) {
			 *cur_in2 = *cur_buffers;
			 *(cur_in2 += idist) = *(cur_buffers += nb);
			 *(cur_in2 += idist) = *(cur_buffers += nb);
			 *(cur_in2 += idist) = *(cur_buffers += nb);
			 cur_buffers += nb;
			 cur_in2 += idist;
		    }
		    for (; buf < nbuffers; ++buf) {
			 *cur_in2 = *cur_buffers;
			 cur_buffers += nb;
			 cur_in2 += idist;
		    }
	       }
	  }

	  /* 
	   * we skip howmany % nbuffers ffts at the end of the loop,
	   * so we have to go back and do them: 
	   */
	  nbuffers = howmany - i;
     } while (i < howmany);
}

/********************* Computing the N-Dimensional FFT *********************/

void fftwnd_aux(fftwnd_plan p, int cur_dim,
		fftw_complex *in, int istride,
		fftw_complex *out, int ostride,
		fftw_complex *work)
{
     int n_after = p->n_after[cur_dim], n = p->n[cur_dim];

     if (cur_dim == p->rank - 2) {
	  /* just do the last dimension directly: */
	  if (p->is_in_place)
	       fftw(p->plans[p->rank - 1], n,
		    in, istride, n_after * istride,
		    work, 1, 0);
	  else
	       fftw(p->plans[p->rank - 1], n,
		    in, istride, n_after * istride,
		    out, ostride, n_after * ostride);
     } else {			/* we have at least two dimensions to go */
	  int i;

	  /* 
	   * process the subsequent dimensions recursively, in hyperslabs,
	   * to get maximum locality: 
	   */
	  for (i = 0; i < n; ++i)
	       fftwnd_aux(p, cur_dim + 1,
			  in + i * n_after * istride, istride,
			  out + i * n_after * ostride, ostride, work);
     }

     /* do the current dimension (in-place): */
     if (p->nbuffers == 0) {
	  fftw(p->plans[cur_dim], n_after,
	       out, n_after * ostride, ostride,
	       work, 1, 0);
     } else			/* using contiguous copy buffers: */
	  fftw_buffered(p->plans[cur_dim], n_after,
			out, n_after * ostride, ostride,
			work, p->nbuffers, work + n);
}

/*
 * alternate version of fftwnd_aux -- this version pushes the howmany
 * loop down to the leaves of the computation, for greater locality in
 * cases where dist < stride
 */
void fftwnd_aux_howmany(fftwnd_plan p, int cur_dim,
			int howmany,
			fftw_complex *in, int istride, int idist,
			fftw_complex *out, int ostride, int odist,
			fftw_complex *work)
{
     int n_after = p->n_after[cur_dim], n = p->n[cur_dim];
     int k;

     if (cur_dim == p->rank - 2) {
	  /* just do the last dimension directly: */
	  if (p->is_in_place)
	       for (k = 0; k < n; ++k)
		    fftw(p->plans[p->rank - 1], howmany,
			 in + k * n_after * istride, istride, idist,
			 work, 1, 0);
	  else
	       for (k = 0; k < n; ++k)
		    fftw(p->plans[p->rank - 1], howmany,
			 in + k * n_after * istride, istride, idist,
			 out + k * n_after * ostride, ostride, odist);
     } else {			/* we have at least two dimensions to go */
	  int i;

	  /* 
	   * process the subsequent dimensions recursively, in
	   * hyperslabs, to get maximum locality:
	   */
	  for (i = 0; i < n; ++i)
	       fftwnd_aux_howmany(p, cur_dim + 1, howmany,
			      in + i * n_after * istride, istride, idist,
				  out + i * n_after * ostride, ostride, odist,
				  work);
     }

     /* do the current dimension (in-place): */
     if (p->nbuffers == 0)
	  for (k = 0; k < n_after; ++k)
	       fftw(p->plans[cur_dim], howmany,
		    out + k * ostride, n_after * ostride, odist,
		    work, 1, 0);
     else			/* using contiguous copy buffers: */
	  for (k = 0; k < n_after; ++k)
	       fftw_buffered(p->plans[cur_dim], howmany,
			     out + k * ostride, n_after * ostride, odist,
			     work, p->nbuffers, work + n);
}

void fftwnd(fftwnd_plan p, int howmany,
	    fftw_complex *in, int istride, int idist,
	    fftw_complex *out, int ostride, int odist)
{
     fftw_complex *work;

#ifdef FFTW_DEBUG
     if (p->rank > 0 && (p->plans[0]->flags & FFTW_THREADSAFE)
	 && p->nwork && p->work)
	  fftw_die("bug with FFTW_THREADSAFE flag");
#endif

     if (p->nwork && !p->work)
	  work = (fftw_complex *) fftw_malloc(p->nwork * sizeof(fftw_complex));
     else
	  work = p->work;

     switch (p->rank) {
	 case 0:
	      break;
	 case 1:
	      if (p->is_in_place)	/* fft is in-place */
		   fftw(p->plans[0], howmany, in, istride, idist,
			work, 1, 0);
	      else
		   fftw(p->plans[0], howmany, in, istride, idist,
			out, ostride, odist);
	      break;
	 default:		/* rank >= 2 */
	      {
		   if (p->is_in_place) {
			out = in;
			ostride = istride;
			odist = idist;
		   }
		   if (howmany > 1 && odist < ostride)
			fftwnd_aux_howmany(p, 0, howmany,
					   in, istride, idist,
					   out, ostride, odist,
					   work);
		   else {
			int i;

			for (i = 0; i < howmany; ++i)
			     fftwnd_aux(p, 0,
					in + i * idist, istride,
					out + i * odist, ostride,
					work);
		   }
	      }
     }

     if (p->nwork && !p->work)
	  fftw_free(work);

}

void fftwnd_one(fftwnd_plan p, fftw_complex *in, fftw_complex *out)
{
     fftwnd(p, 1, in, 1, 1, out, 1, 1);
}
