
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

/* $Id$ */

#include <fftw-int.h>
#include <rfftw.h>

/********************** prototypes for rexec2 routines **********************/

extern void rfftw_real2c_aux(fftw_plan plan, int howmany,
			     fftw_real *in, int istride, int idist,
			     fftw_complex *out, int ostride, int odist,
			     fftw_real *work);
extern void rfftw_c2real_aux(fftw_plan plan, int howmany,
			     fftw_complex *in, int istride, int idist,
			     fftw_real *out, int ostride, int odist,
			     fftw_real *work);
extern void rfftw_real2c_overlap_aux(fftw_plan plan, int howmany,
				   fftw_real *in, int istride, int idist,
			       fftw_complex *out, int ostride, int odist,
				     fftw_real *work);
extern void rfftw_c2real_overlap_aux(fftw_plan plan, int howmany,
				fftw_complex *in, int istride, int idist,
				  fftw_real *out, int ostride, int odist,
				     fftw_real *work);

/********************** Initializing the RFFTWND Plan ***********************/

/*
 * Create an fftwnd_plan specialized for specific arrays.  (These
 * arrays are ignored, however, if they are NULL or if the flags
 * do not include FFTW_MEASURE.)  The main advantage of being
 * provided arrays like this is that we can do runtime timing
 * measurements of our options, without worrying about allocating
 * excessive scratch space. 
 */
fftwnd_plan rfftwnd_create_plan_specific(int rank, const int *n,
					 fftw_direction dir, int flags,
					 fftw_real *in, int istride,
					 fftw_real *out, int ostride)
{
     fftwnd_plan p;
     int i;
     int rflags = flags & ~FFTW_IN_PLACE;
     /* note that we always do rfftw transforms out-of-place in rexec2.c */

     if (flags & FFTW_IN_PLACE) {
	  out = NULL;
	  ostride = istride;
     }
     istride = ostride = 1;	/* 
				 * strides don't work yet, since it is not 
				 * clear whether they apply to real 
				 * or complex data 
				 */

     if (!(p = fftwnd_create_plan_aux(rank, n, dir, flags)))
	  return 0;

     for (i = 0; i < rank - 1; ++i)
	  p->n_after[i] = (n[rank - 1]/2 + 1) * (p->n_after[i] / n[rank - 1]);
     if (rank > 0)
	  p->n[rank - 1] = n[rank - 1] / 2 + 1;

     p->plans = fftwnd_new_plan_array(rank);
     if (rank > 0 && !p->plans) {
	  rfftwnd_destroy_plan(p);
	  return 0;
     }
     if (rank > 0) {
	  p->plans[rank - 1] = rfftw_create_plan(n[rank - 1], dir, rflags);
	  if (!p->plans[rank - 1]) {
	       rfftwnd_destroy_plan(p);
	       return 0;
	  }
     }
     if (rank > 1) {
	  if (!(flags & FFTW_MEASURE) || in == 0
	      || (!p->is_in_place && out == 0)) {
	       if (!fftwnd_create_plans_generic(p->plans, rank - 1, n,
					   dir, flags | FFTW_IN_PLACE)) {
		    rfftwnd_destroy_plan(p);
		    return 0;
	       }
	  } else if (dir == FFTW_COMPLEX_TO_REAL || (flags & FFTW_IN_PLACE)) {
	       if (!fftwnd_create_plans_specific(p->plans, rank - 1, n,
						 p->n_after,
					      dir, flags | FFTW_IN_PLACE,
						 (fftw_complex *) in,
						 istride,
						 0, 0)) {
		    rfftwnd_destroy_plan(p);
		    return 0;
	       }
	  } else {
	       if (!fftwnd_create_plans_specific(p->plans, rank - 1, n,
						 p->n_after,
					      dir, flags | FFTW_IN_PLACE,
						 (fftw_complex *) out,
						 ostride,
						 0, 0)) {
		    rfftwnd_destroy_plan(p);
		    return 0;
	       }
	  }
     }
     p->nbuffers = 0;
     p->nwork = fftwnd_work_size(rank, p->n, flags | FFTW_IN_PLACE,
				 p->nbuffers + 1);
     if (p->nwork && !(flags & FFTW_THREADSAFE)) {
	  p->work = (fftw_complex *) fftw_malloc(p->nwork
						 * sizeof(fftw_complex));
	  if (!p->work) {
	       rfftwnd_destroy_plan(p);
	       return 0;
	  }
     }
     return p;
}

fftwnd_plan rfftw2d_create_plan_specific(int nx, int ny,
					 fftw_direction dir, int flags,
					 fftw_real *in, int istride,
					 fftw_real *out, int ostride)
{
     int n[2];

     n[0] = nx;
     n[1] = ny;

     return rfftwnd_create_plan_specific(2, n, dir, flags,
					 in, istride, out, ostride);
}

fftwnd_plan rfftw3d_create_plan_specific(int nx, int ny, int nz,
					 fftw_direction dir, int flags,
					 fftw_real *in, int istride,
					 fftw_real *out, int ostride)
{
     int n[3];

     n[0] = nx;
     n[1] = ny;
     n[2] = nz;

     return rfftwnd_create_plan_specific(3, n, dir, flags,
					 in, istride, out, ostride);
}

/* Create a generic fftwnd plan: */

fftwnd_plan rfftwnd_create_plan(int rank, const int *n,
				fftw_direction dir, int flags)
{
     return rfftwnd_create_plan_specific(rank, n, dir, flags, 0, 1, 0, 1);
}

fftwnd_plan rfftw2d_create_plan(int nx, int ny,
				fftw_direction dir, int flags)
{
     return rfftw2d_create_plan_specific(nx, ny, dir, flags, 0, 1, 0, 1);
}

fftwnd_plan rfftw3d_create_plan(int nx, int ny, int nz,
				fftw_direction dir, int flags)
{
     return rfftw3d_create_plan_specific(nx, ny, nz, dir, flags, 0, 1, 0, 1);
}

/************************ Freeing the RFFTWND Plan ************************/

void rfftwnd_destroy_plan(fftwnd_plan plan)
{
     fftwnd_destroy_plan(plan);
}

/************************ Printing the RFFTWND Plan ************************/

void rfftwnd_fprint_plan(FILE *f, fftwnd_plan plan)
{
     fftwnd_fprint_plan(f, plan);
}

void rfftwnd_print_plan(fftwnd_plan plan)
{
     rfftwnd_fprint_plan(stdout, plan);
}

/*********** Computing the N-Dimensional FFT: Auxiliary Routines ************/

void rfftwnd_real2c_aux(fftwnd_plan p, int cur_dim,
			fftw_real *in, int istride,
			fftw_complex *out, int ostride,
			fftw_real *work)
{
     int n_after = p->n_after[cur_dim], n = p->n[cur_dim];

     if (cur_dim == p->rank - 2) {
	  /* just do the last dimension directly: */
	  if (p->is_in_place)
	       rfftw_real2c_aux(p->plans[p->rank - 1], n,
				in, istride, (n_after * istride) * 2,
				out, istride, n_after * istride,
				work);
	  else
	       rfftw_real2c_aux(p->plans[p->rank - 1], n,
			 in, istride, p->plans[p->rank - 1]->n * istride,
				out, ostride, n_after * ostride,
				work);
     } else {			/* we have at least two dimensions to go */
	  int nr = p->plans[p->rank - 1]->n;
	  int n_after_r = p->is_in_place ? n_after * 2 
	       : nr * (n_after / (nr/2 + 1));
	  int i;

	  /* 
	   * process the subsequent dimensions recursively, in hyperslabs,
	   * to get maximum locality: 
	   */
	  for (i = 0; i < n; ++i)
	       rfftwnd_real2c_aux(p, cur_dim + 1,
				  in + i * n_after_r * istride, istride,
			     out + i * n_after * ostride, ostride, work);
     }

     /* do the current dimension (in-place): */
     fftw(p->plans[cur_dim], n_after,
	  out, n_after * ostride, ostride,
	  (fftw_complex *) work, 1, 0);
     /* I hate this cast */
}

void rfftwnd_c2real_aux(fftwnd_plan p, int cur_dim,
			fftw_complex *in, int istride,
			fftw_real *out, int ostride,
			fftw_real *work)
{
     int n_after = p->n_after[cur_dim], n = p->n[cur_dim];

     /* do the current dimension (in-place): */
     fftw(p->plans[cur_dim], n_after,
	  in, n_after * istride, istride,
	  (fftw_complex *) work, 1, 0);

     if (cur_dim == p->rank - 2) {
	  /* just do the last dimension directly: */
	  if (p->is_in_place)
	       rfftw_c2real_aux(p->plans[p->rank - 1], n,
				in, istride, n_after * istride,
				out, istride, (n_after * istride) * 2,
				work);
	  else
	       rfftw_c2real_aux(p->plans[p->rank - 1], n,
				in, istride, n_after * istride,
			out, ostride, p->plans[p->rank - 1]->n * ostride,
				work);
     } else {			/* we have at least two dimensions to go */
	  int nr = p->plans[p->rank - 1]->n;
	  int n_after_r = p->is_in_place ? n_after * 2 : 
	       nr * (n_after / (nr/2 + 1));
	  int i;

	  /* 
	   * process the subsequent dimensions recursively, in hyperslabs,
	   * to get maximum locality: 
	   */
	  for (i = 0; i < n; ++i)
	       rfftwnd_c2real_aux(p, cur_dim + 1,
				  in + i * n_after * istride, istride,
			   out + i * n_after_r * ostride, ostride, work);
     }
}

/*
 * alternate version of rfftwnd_aux -- this version pushes the howmany
 * loop down to the leaves of the computation, for greater locality
 * in cases where dist < stride.  It is also required for correctness
 * if in==out, and we must call a special version of the executor.
 * Note that work must point to 'howmany' copies of its data
 * if in == out. 
 */

void rfftwnd_real2c_aux_howmany(fftwnd_plan p, int cur_dim,
				int howmany,
				fftw_real *in, int istride, int idist,
				fftw_complex *out, int ostride, int odist,
				fftw_complex *work)
{
     int n_after = p->n_after[cur_dim], n = p->n[cur_dim];
     int k;

     if (cur_dim == p->rank - 2) {
	  /* just do the last dimension directly: */
	  if (p->is_in_place)
	       for (k = 0; k < n; ++k)
		    rfftw_real2c_overlap_aux(p->plans[p->rank - 1], howmany,
					in + (k * n_after * istride) * 2,
					     istride, idist,
					   out + (k * n_after * ostride),
					     ostride, odist,
					     (fftw_real *) work);
	  else {
	       int nlast = p->plans[p->rank - 1]->n;
	       for (k = 0; k < n; ++k)
		    rfftw_real2c_aux(p->plans[p->rank - 1], howmany,
				     in + k * nlast * istride,
				     istride, idist,
				     out + k * n_after * ostride,
				     ostride, odist,
				     (fftw_real *) work);
	  }
     } else {			/* we have at least two dimensions to go */
	  int nr = p->plans[p->rank - 1]->n;
	  int n_after_r = p->is_in_place ? n_after * 2 : 
	       nr * (n_after / (nr/2 + 1));
	  int i;

	  /* 
	   * process the subsequent dimensions recursively, in hyperslabs,
	   * to get maximum locality: 
	   */
	  for (i = 0; i < n; ++i)
	       rfftwnd_real2c_aux_howmany(p, cur_dim + 1, howmany,
			    in + i * n_after_r * istride, istride, idist,
			     out + i * n_after * ostride, ostride, odist,
					  work);
     }

     /* do the current dimension (in-place): */
     for (k = 0; k < n_after; ++k)
	  fftw(p->plans[cur_dim], howmany,
	       out + k * ostride, n_after * ostride, odist,
	       work, 1, 0);
}

void rfftwnd_c2real_aux_howmany(fftwnd_plan p, int cur_dim,
				int howmany,
				fftw_complex *in, int istride, int idist,
				fftw_real *out, int ostride, int odist,
				fftw_complex *work)
{
     int n_after = p->n_after[cur_dim], n = p->n[cur_dim];
     int k;

     /* do the current dimension (in-place): */
     for (k = 0; k < n_after; ++k)
	  fftw(p->plans[cur_dim], howmany,
	       in + k * istride, n_after * istride, idist,
	       work, 1, 0);

     if (cur_dim == p->rank - 2) {
	  /* just do the last dimension directly: */
	  if (p->is_in_place)
	       for (k = 0; k < n; ++k)
		    rfftw_c2real_overlap_aux(p->plans[p->rank - 1], howmany,
					     in + (k * n_after * istride),
					     istride, idist,
				       out + (k * n_after * ostride) * 2,
					     ostride, odist,
					     (fftw_real *) work);
	  else {
	       int nlast = p->plans[p->rank - 1]->n;
	       for (k = 0; k < n; ++k)
		    rfftw_c2real_aux(p->plans[p->rank - 1], howmany,
				     in + k * n_after * istride,
				     istride, idist,
				     out + k * nlast * ostride,
				     ostride, odist,
				     (fftw_real *) work);
	  }
     } else {			/* we have at least two dimensions to go */
	  int nr = p->plans[p->rank - 1]->n;
	  int n_after_r = p->is_in_place ? n_after * 2
	       : nr * (n_after / (nr/2 + 1));
	  int i;

	  /* 
	   * process the subsequent dimensions recursively, in hyperslabs,
	   * to get maximum locality: 
	   */
	  for (i = 0; i < n; ++i)
	       rfftwnd_c2real_aux_howmany(p, cur_dim + 1, howmany,
			      in + i * n_after * istride, istride, idist,
			   out + i * n_after_r * ostride, ostride, odist,
					  work);
     }
}

/********** Computing the N-Dimensional FFT: User-Visible Routines **********/

void rfftwnd_real_to_complex(fftwnd_plan p, int howmany,
			     fftw_real *in, int istride, int idist,
			     fftw_complex *out, int ostride, int odist)
{
     fftw_complex *work = p->work;
     int rank = p->rank;
     int free_work = 0;

     if (p->dir != FFTW_REAL_TO_COMPLEX)
	  fftw_die("rfftwnd_real_to_complex with complex-to-real plan");

#ifdef FFTW_DEBUG
     if (p->rank > 0 && (p->plans[0]->flags & FFTW_THREADSAFE)
	 && p->nwork && p->work)
	  fftw_die("bug with FFTW_THREADSAFE flag");
#endif

     if (p->is_in_place) {
	  ostride = istride;
	  odist = (idist == 1) ? 1 : (idist / 2);	/* ugh */
	  out = (fftw_complex *) in;
	  if (howmany > 1 && istride > idist && rank > 0) {
	       int new_nwork;

	       new_nwork = p->n[rank - 1] * howmany;
	       if (new_nwork > p->nwork) {
		    work = (fftw_complex *)
			fftw_malloc(sizeof(fftw_complex) * new_nwork);
		    if (!work)
			 fftw_die("error allocating work array");
		    free_work = 1;
	       }
	  }
     }
     if (p->nwork && !work) {
	  work = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * p->nwork);
	  free_work = 1;
     }
     switch (rank) {
	 case 0:
	      break;
	 case 1:
	      if (p->is_in_place && howmany > 1 && istride > idist)
		   rfftw_real2c_overlap_aux(p->plans[0], howmany,
					    in, istride, idist,
					    out, ostride, odist,
					    (fftw_real *) work);
	      else
		   rfftw_real2c_aux(p->plans[0], howmany,
				    in, istride, idist,
				    out, ostride, odist,
				    (fftw_real *) work);
	      break;
	 default:		/* rank >= 2 */
	      {
		   if (howmany > 1 && ostride > odist)
			rfftwnd_real2c_aux_howmany(p, 0, howmany,
						   in, istride, idist,
						   out, ostride, odist,
						   work);
		   else {
			int i;

			for (i = 0; i < howmany; ++i)
			     rfftwnd_real2c_aux(p, 0,
						in + i * idist, istride,
						out + i * odist, ostride,
						(fftw_real *) work);
		   }
	      }
     }

     if (free_work)
	  fftw_free(work);
}

void rfftwnd_complex_to_real(fftwnd_plan p, int howmany,
			     fftw_complex *in, int istride, int idist,
			     fftw_real *out, int ostride, int odist)
{
     fftw_complex *work = p->work;
     int rank = p->rank;
     int free_work = 0;

     if (p->dir != FFTW_COMPLEX_TO_REAL)
	  fftw_die("rfftwnd_complex_to_real with real-to-complex plan");

#ifdef FFTW_DEBUG
     if (p->rank > 0 && (p->plans[0]->flags & FFTW_THREADSAFE)
	 && p->nwork && p->work)
	  fftw_die("bug with FFTW_THREADSAFE flag");
#endif

     if (p->is_in_place) {
	  ostride = istride;
	  odist = idist;
	  odist = (idist == 1) ? 1 : (idist * 2);	/* ugh */
	  out = (fftw_real *) in;
	  if (howmany > 1 && istride > idist && rank > 0) {
	       int new_nwork = p->n[rank - 1] * howmany;
	       if (new_nwork > p->nwork) {
		    work = (fftw_complex *)
			fftw_malloc(sizeof(fftw_complex) * new_nwork);
		    if (!work)
			 fftw_die("error allocating work array");
		    free_work = 1;
	       }
	  }
     }
     if (p->nwork && !work) {
	  work = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * p->nwork);
	  free_work = 1;
     }
     switch (rank) {
	 case 0:
	      break;
	 case 1:
	      if (p->is_in_place && howmany > 1 && istride > idist)
		   rfftw_c2real_overlap_aux(p->plans[0], howmany,
					    in, istride, idist,
					    out, ostride, odist,
					    (fftw_real *) work);
	      else
		   rfftw_c2real_aux(p->plans[0], howmany,
				    in, istride, idist,
				    out, ostride, odist,
				    (fftw_real *) work);
	      break;
	 default:		/* rank >= 2 */
	      {
		   if (howmany > 1 && ostride > odist)
			rfftwnd_c2real_aux_howmany(p, 0, howmany,
						   in, istride, idist,
						   out, ostride, odist,
						   work);
		   else {
			int i;

			for (i = 0; i < howmany; ++i)
			     rfftwnd_c2real_aux(p, 0,
						in + i * idist, istride,
						out + i * odist, ostride,
						(fftw_real *) work);
		   }
	      }
     }

     if (free_work)
	  fftw_free(work);
}

void rfftwnd_one_real_to_complex(fftwnd_plan p,
				 fftw_real *in, fftw_complex *out)
{
     rfftwnd_real_to_complex(p, 1, in, 1, 1, out, 1, 1);
}

void rfftwnd_one_complex_to_real(fftwnd_plan p,
				 fftw_complex *in, fftw_real *out)
{
     rfftwnd_complex_to_real(p, 1, in, 1, 1, out, 1, 1);
}
