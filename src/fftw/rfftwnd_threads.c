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

#include <fftw_threads-int.h>
#include <rfftw_threads.h>

/****************** prototypes for rexec2_threads routines *****************/

extern void rfftw_real2c_threads_aux(fftw_plan plan, int howmany,
				     fftw_real *in, int istride, int idist,
				     fftw_complex *out, int ostride, int odist,
				     fftw_real *work,
				     int nthreads);
extern void rfftw_c2real_threads_aux(fftw_plan plan, int howmany,
				     fftw_complex *in, int istride, int idist,
				     fftw_real *out, int ostride, int odist,
				     fftw_real *work,
				     int nthreads);
extern void rfftw_real2c_overlap_threads_aux(fftw_plan plan, int howmany,
					     fftw_real *in, int istride,
					     int idist,
					     fftw_complex *out,
					     int ostride, int odist,
					     fftw_real *work,
					     int nthreads);
extern void rfftw_c2real_overlap_threads_aux(fftw_plan plan, int howmany,
					     fftw_complex *in,
					     int istride, int idist,
					     fftw_real *out, 
					     int ostride, int odist,
					     fftw_real *work,
					     int nthreads);

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

/****************** prototypes for rfftwnd routines *****************/

extern void rfftwnd_real2c_aux(fftwnd_plan p, int cur_dim,
			       fftw_real *in, int istride,
			       fftw_complex *out, int ostride,
			       fftw_real *work);
extern void rfftwnd_c2real_aux(fftwnd_plan p, int cur_dim,
			       fftw_complex *in, int istride,
			       fftw_real *out, int ostride,
			       fftw_real *work);
extern void rfftwnd_real2c_aux_howmany(fftwnd_plan p, int cur_dim,
                                int howmany,
                                fftw_real *in, int istride, int idist,
                                fftw_complex *out, int ostride, int odist,
                                fftw_complex *work);
extern void rfftwnd_c2real_aux_howmany(fftwnd_plan p, int cur_dim,
                                int howmany,
                                fftw_complex *in, int istride, int idist,
                                fftw_real *out, int ostride, int odist,
                                fftw_complex *work);


/*********** Computing the N-Dimensional FFT: Auxiliary Routines ************/

typedef struct {
     fftwnd_plan p;
     int cur_dim;
     void *in;
     int istride, idist;
     void *out;
     int ostride, odist;
     fftw_complex *work;
} aux_data;

static void *real2c_aux_thread(fftw_loop_data *ldata)
{
     int min = ldata->min, max = ldata->max;
     aux_data *d = (aux_data *) ldata->data;
     fftwnd_plan p = d->p;
     int cur_dim = d->cur_dim;
     fftw_real *in = (fftw_real *) d->in;
     int istride = d->istride, idist = d->idist;
     fftw_complex *out = (fftw_complex *) d->out;
     int ostride = d->ostride, odist = d->odist;
     fftw_real *work = (fftw_real*) (d->work + p->nwork * ldata->thread_num);

     for (; min < max; ++min)
	  rfftwnd_real2c_aux(p, cur_dim, in + idist * min, istride,
			     out + odist * min, ostride, work);
     return 0;
}

void rfftwnd_real2c_threads_aux(fftwnd_plan p, int cur_dim,
				fftw_real *in, int istride,
				fftw_complex *out, int ostride,
				fftw_complex *work,
				int nthreads)
{
     int n_after = p->n_after[cur_dim], n = p->n[cur_dim];

     if (cur_dim == p->rank - 2) {
	  /* just do the last dimension directly: */
	  if (p->is_in_place)
	       rfftw_real2c_threads_aux(p->plans[p->rank - 1], n,
					in, istride, (n_after * istride) * 2,
					out, istride, n_after * istride,
					(fftw_real *) work, nthreads);
	  else
	       rfftw_real2c_threads_aux(p->plans[p->rank - 1], n,
					in, istride,
					p->plans[p->rank - 1]->n * istride,
					out, ostride, n_after * ostride,
					(fftw_real *) work, nthreads);
     }
     else {    /* we have at least two dimensions to go */
	  int nr = p->plans[p->rank - 1]->n;
	  aux_data d;

	  d.p = p;
	  d.cur_dim = cur_dim + 1;
	  d.in = in;
	  d.istride = istride;
	  d.idist = istride * (p->is_in_place ? n_after * 2
			       : nr * (n_after / (nr/2 + 1)));
	  d.out = out;
	  d.ostride = ostride;
	  d.odist = ostride * n_after;
	  d.work = work;

	  fftw_thread_spawn_loop(n, nthreads, real2c_aux_thread, &d);
     }

     /* do the current dimension (in-place): */
     /* (Use internal function instead of fftw_threads so that we can
	pass our workspace array.) */
     fftw_executor_many_inplace_threads(p->plans[cur_dim]->n,
					out, work, p->plans[cur_dim]->root,
					n_after * ostride, n_after, ostride,
					nthreads);
}

static void *c2real_aux_thread(fftw_loop_data *ldata)
{
     int min = ldata->min, max = ldata->max;
     aux_data *d = (aux_data *) ldata->data;
     fftwnd_plan p = d->p;
     int cur_dim = d->cur_dim;
     fftw_complex *in = (fftw_complex *) d->in;
     int istride = d->istride, idist = d->idist;
     fftw_real *out = (fftw_real *) d->out;
     int ostride = d->ostride, odist = d->odist;
     fftw_real *work = (fftw_real*) (d->work + p->nwork * ldata->thread_num);

     for (; min < max; ++min)
	  rfftwnd_c2real_aux(p, cur_dim, in + idist * min, istride,
			     out + odist * min, ostride, work);
     return 0;
}

void rfftwnd_c2real_threads_aux(fftwnd_plan p, int cur_dim,
				fftw_complex *in, int istride,
				fftw_real *out, int ostride,
				fftw_complex *work, int nthreads)
{
     int n_after = p->n_after[cur_dim], n = p->n[cur_dim];

     /* do the current dimension (in-place): */
     /* (Use internal function instead of fftw_threads so that we can
	pass our workspace array.) */
     fftw_executor_many_inplace_threads(p->plans[cur_dim]->n,
					in, work, p->plans[cur_dim]->root,
					n_after * istride, n_after, istride,
					nthreads);

     if (cur_dim == p->rank - 2) {
	  /* just do the last dimension directly: */
	  if (p->is_in_place)
	       rfftw_c2real_threads_aux(p->plans[p->rank - 1], n,
					in, istride, n_after * istride,
					out, istride, (n_after * istride) * 2,
					(fftw_real *) work, nthreads);
	  else
	       rfftw_c2real_threads_aux(p->plans[p->rank - 1], n,
					in, istride, n_after * istride,
					out, ostride,
					p->plans[p->rank - 1]->n * ostride,
					(fftw_real *) work, nthreads);
     }
     else {	/* we have at least two dimensions to go */
	  int nr = p->plans[p->rank - 1]->n;
	  aux_data d;

	  d.p = p;
	  d.cur_dim = cur_dim + 1;
	  d.in = in;
	  d.istride = istride;
	  d.odist = ostride * (p->is_in_place ? n_after * 2
			       : nr * (n_after / (nr/2 + 1)));
	  d.out = out;
	  d.ostride = ostride;
	  d.idist = istride * n_after;
	  d.work = work;

	  fftw_thread_spawn_loop(n, nthreads, c2real_aux_thread, &d);
     }
}

typedef struct {
     fftw_plan p;
     int howmany;
     void *in;
     int istride, idist, idist0;
     void *out;
     int ostride, odist, odist0;
     fftw_real *work;
     int wdist;
} howmany_aux_data;

static void *r2c_overlap_howmany_thread(fftw_loop_data *ldata)
{
     int min = ldata->min, max = ldata->max;
     howmany_aux_data *d = (howmany_aux_data *) ldata->data;
     fftw_plan p = d->p;
     int howmany = d->howmany;
     fftw_real *in = (fftw_real *) d->in;
     int istride = d->istride, idist = d->idist, idist0 = d->idist0;
     fftw_complex *out = (fftw_complex *) d->out;
     int ostride = d->ostride, odist = d->odist, odist0 = d->odist0;
     fftw_real *work = d->work + d->wdist * ldata->thread_num;

     for (; min < max; ++min)
	  rfftw_real2c_overlap_aux(p, howmany,
				   in + min * idist0, istride, idist,
				   out + min * odist0, ostride, odist,
				   work);

     return 0;
}

static void *c2r_overlap_howmany_thread(fftw_loop_data *ldata)
{
     int min = ldata->min, max = ldata->max;
     howmany_aux_data *d = (howmany_aux_data *) ldata->data;
     fftw_plan p = d->p;
     int howmany = d->howmany;
     fftw_real *out = (fftw_real *) d->out;
     int istride = d->istride, idist = d->idist, idist0 = d->idist0;
     fftw_complex *in = (fftw_complex *) d->in;
     int ostride = d->ostride, odist = d->odist, odist0 = d->odist0;
     fftw_real *work = d->work + d->wdist * ldata->thread_num;

     for (; min < max; ++min)
	  rfftw_c2real_overlap_aux(p, howmany,
				   in + min * idist0, istride, idist,
				   out + min * odist0, ostride, odist,
				   work);

     return 0;
}

static void *r2c_howmany_thread(fftw_loop_data *ldata)
{
     int min = ldata->min, max = ldata->max;
     howmany_aux_data *d = (howmany_aux_data *) ldata->data;
     fftw_plan p = d->p;
     int howmany = d->howmany;
     fftw_real *in = (fftw_real *) d->in;
     int istride = d->istride, idist = d->idist, idist0 = d->idist0;
     fftw_complex *out = (fftw_complex *) d->out;
     int ostride = d->ostride, odist = d->odist, odist0 = d->odist0;
     fftw_real *work = d->work + d->wdist * ldata->thread_num;

     for (; min < max; ++min)
	  rfftw_real2c_aux(p, howmany,
			   in + min * idist0, istride, idist,
			   out + min * odist0, ostride, odist,
			   work);

     return 0;
}

static void *c2r_howmany_thread(fftw_loop_data *ldata)
{
     int min = ldata->min, max = ldata->max;
     howmany_aux_data *d = (howmany_aux_data *) ldata->data;
     fftw_plan p = d->p;
     int howmany = d->howmany;
     fftw_real *out = (fftw_real *) d->out;
     int istride = d->istride, idist = d->idist, idist0 = d->idist0;
     fftw_complex *in = (fftw_complex *) d->in;
     int ostride = d->ostride, odist = d->odist, odist0 = d->odist0;
     fftw_real *work = d->work + d->wdist * ldata->thread_num;

     for (; min < max; ++min)
	  rfftw_c2real_aux(p, howmany,
			   in + min * idist0, istride, idist,
			   out + min * odist0, ostride, odist,
			   work);

     return 0;
}

typedef struct {
     fftwnd_plan p;
     int cur_dim;
     int howmany;
     void *in;
     int istride, idist, idist0;
     void *out;
     int ostride, odist, odist0;
     fftw_complex *work;
     int wdist;
} howmany_hyperslab_aux_data;

static void *r2c_hyperslab_howmany_thread(fftw_loop_data *ldata)
{
     int min = ldata->min, max = ldata->max;
     howmany_hyperslab_aux_data *d = (howmany_hyperslab_aux_data*) ldata->data;
     fftwnd_plan p = d->p;
     int cur_dim = d->cur_dim;
     int howmany = d->howmany;
     fftw_real *in = (fftw_real *) d->in;
     int istride = d->istride, idist = d->idist, idist0 = d->idist0;
     fftw_complex *out = (fftw_complex *) d->out;
     int ostride = d->ostride, odist = d->odist, odist0 = d->odist0;
     fftw_complex *work = d->work + d->wdist * ldata->thread_num;

     for (; min < max; ++min)
	  rfftwnd_real2c_aux_howmany(p, cur_dim, howmany,
				     in + min * idist0, istride, idist,
				     out + min * odist0, ostride, odist,
				     work);

     return 0;
}

static void *c2r_hyperslab_howmany_thread(fftw_loop_data *ldata)
{
     int min = ldata->min, max = ldata->max;
     howmany_hyperslab_aux_data *d = (howmany_hyperslab_aux_data*) ldata->data;
     fftwnd_plan p = d->p;
     int cur_dim = d->cur_dim;
     int howmany = d->howmany;
     fftw_real *out = (fftw_real *) d->out;
     int istride = d->istride, idist = d->idist, idist0 = d->idist0;
     fftw_complex *in = (fftw_complex *) d->in;
     int ostride = d->ostride, odist = d->odist, odist0 = d->odist0;
     fftw_complex *work = d->work + d->wdist * ldata->thread_num;

     for (; min < max; ++min)
	  rfftwnd_c2real_aux_howmany(p, cur_dim, howmany,
				     in + min * idist0, istride, idist,
				     out + min * odist0, ostride, odist,
				     work);

     return 0;
}

typedef struct {
     fftw_plan p;
     int howmany;
     fftw_complex *io_data;
     int iostride, iodist, iodist0;
     fftw_complex *work;
     int wdist;
} fftw_howmany_data;

static void *fftw_howmany_thread(fftw_loop_data *ldata)
{
     int min = ldata->min, max = ldata->max;
     fftw_howmany_data *d = (fftw_howmany_data*) ldata->data;
     fftw_plan p = d->p;
     int howmany = d->howmany;
     fftw_complex *io_data = d->io_data;
     int iostride = d->iostride, iodist = d->iodist, iodist0 = d->iodist0;
     fftw_complex *work = d->work + d->wdist * ldata->thread_num;

     for (; min < max; ++min)
	  fftw(p, howmany, io_data + min*iodist0, iostride, iodist, work,1,0);

     return 0;
}

/*
 * alternate version of rfftwnd_aux -- this version pushes the howmany
 * loop down to the leaves of the computation, for greater locality
 * in cases where dist < stride.  It is also required for correctness
 * if in==out, and we must call a special version of the executor.
 * Note that work must point to 'howmany' copies of its data
 * if in == out. 
 */

void rfftwnd_real2c_aux_howmany_threads(fftwnd_plan p, int cur_dim,
					int howmany,
					fftw_real *in, int istride, int idist,
					fftw_complex *out,
					int ostride, int odist,
					fftw_complex *work, int nwork,
					int nthreads)
{
     int n_after = p->n_after[cur_dim], n = p->n[cur_dim];

     if (cur_dim == p->rank - 2) {
	  howmany_aux_data d;

	  d.p = p->plans[p->rank - 1];
	  d.howmany = howmany;
	  d.in = in;
	  d.istride = istride; d.idist = idist;
	  d.out = out;
	  d.ostride = ostride; d.odist = odist;
	  d.work = (fftw_real *) work;
	  d.wdist = nwork * 2;

	  /* just do the last dimension directly: */
	  if (p->is_in_place) {
	       d.idist0 =  n_after * istride * 2;
	       d.odist0 = n_after * ostride;
	       fftw_thread_spawn_loop(n, nthreads,
				      r2c_overlap_howmany_thread, &d);
	  }
	  else {
	       d.idist0 = p->plans[p->rank - 1]->n * istride;
	       d.odist0 = n_after * ostride;
	       fftw_thread_spawn_loop(n, nthreads,
				      r2c_howmany_thread, &d);
	  }
     } 
     else {      /* we have at least two dimensions to go */
	  /* 
	   * process the subsequent dimensions recursively, in hyperslabs,
	   * to get maximum locality: 
	   */

	  int nr = p->plans[p->rank - 1]->n;
	  int n_after_r = p->is_in_place ? n_after * 2 : 
	       nr * (n_after / (nr/2 + 1));
	  howmany_hyperslab_aux_data d;

	  d.p = p;
	  d.cur_dim = cur_dim + 1;
	  d.howmany = howmany;
	  d.in = in;
	  d.istride = istride;
	  d.idist = idist;
	  d.idist0 = n_after_r * istride;
	  d.out = out;
	  d.ostride = ostride;
	  d.odist = odist;
	  d.odist0 = n_after * ostride;
	  d.work = work;
	  d.wdist = nwork;

	  fftw_thread_spawn_loop(n, nthreads,
				 r2c_hyperslab_howmany_thread, &d);
     }

     /* do the current dimension (in-place): */
     {
	  fftw_howmany_data d;

	  d.p = p->plans[cur_dim];
	  d.howmany = howmany;
	  d.io_data = out;
	  d.iostride = n_after * ostride;
	  d.iodist = odist;
	  d.iodist0 = ostride;
	  d.work =  work;
	  d.wdist = nwork;

	  fftw_thread_spawn_loop(n_after, nthreads, fftw_howmany_thread, &d);
     }
}

void rfftwnd_c2real_aux_howmany_threads(fftwnd_plan p, int cur_dim,
					int howmany,
					fftw_complex *in,
					int istride, int idist,
					fftw_real *out,
					int ostride, int odist,
					fftw_complex *work, int nwork,
					int nthreads)
{
     int n_after = p->n_after[cur_dim], n = p->n[cur_dim];

     /* do the current dimension (in-place): */
     {
          fftw_howmany_data d;

          d.p = p->plans[cur_dim];
          d.howmany = howmany;
          d.io_data = in;
          d.iostride = n_after * istride;
          d.iodist = idist;
          d.iodist0 = istride;
          d.work = work;
	  d.wdist = nwork;

          fftw_thread_spawn_loop(n_after, nthreads, fftw_howmany_thread, &d);
     }

     if (cur_dim == p->rank - 2) {
          howmany_aux_data d;

          d.p = p->plans[p->rank - 1];
          d.howmany = howmany;
          d.in = in;
          d.istride = istride; d.idist = idist;
          d.out = out;
          d.ostride = ostride; d.odist = odist;
          d.work = (fftw_real *) work;
	  d.wdist = nwork * 2;

	  /* just do the last dimension directly: */
          if (p->is_in_place) {
               d.idist0 = n_after * istride;
               d.odist0 = n_after * ostride * 2;
               fftw_thread_spawn_loop(n, nthreads,
                                      c2r_overlap_howmany_thread, &d);
          }
          else {
               d.odist0 = p->plans[p->rank - 1]->n * ostride;
               d.idist0 = n_after * istride;
               fftw_thread_spawn_loop(n, nthreads,
                                      c2r_howmany_thread, &d);
          }
     } 
     else {			/* we have at least two dimensions to go */
          /*
           * process the subsequent dimensions recursively, in hyperslabs,
           * to get maximum locality:
           */

          int nr = p->plans[p->rank - 1]->n;
          int n_after_r = p->is_in_place ? n_after * 2 :
               nr * (n_after / (nr/2 + 1));
          howmany_hyperslab_aux_data d;

          d.p = p;
          d.cur_dim = cur_dim + 1;
          d.howmany = howmany;
          d.in = in;
          d.istride = istride;
          d.idist = idist;
          d.idist0 = n_after * istride;
          d.out = out;
          d.ostride = ostride;
          d.odist = odist;
          d.odist0 = n_after_r * ostride;
          d.work = work;
	  d.wdist = nwork;

          fftw_thread_spawn_loop(n, nthreads,
                                 c2r_hyperslab_howmany_thread, &d);
     }
}

/********** Computing the N-Dimensional FFT: User-Visible Routines **********/

void rfftwnd_threads_real_to_complex(int nthreads, fftwnd_plan p, int howmany,
				     fftw_real *in, int istride, int idist,
				     fftw_complex *out, int ostride, int odist)
{
     fftw_complex *work = 0;
     int rank = p->rank;
     int nwork = p->nwork, size_work = nwork * nthreads;

     if (p->dir != FFTW_REAL_TO_COMPLEX)
	  fftw_die("rfftwnd_real_to_complex with complex-to-real plan");

     if (p->is_in_place) {
	  ostride = istride;
	  odist = (idist == 1) ? 1 : (idist / 2);	/* ugh */
	  out = (fftw_complex *) in;
	  if (howmany > 1 && istride > idist && rank > 0) {
	       int new_nwork = p->n[rank - 1] * howmany;
	       if (new_nwork > nwork)
		    nwork = new_nwork;
	       if (rank != 1) {
		    if (nwork * nthreads > size_work)
			 size_work = nwork * nthreads;
	       }
	       else
		    size_work = nwork;
	  }
     }

     work = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * size_work);

     switch (rank) {
	 case 0:
	      break;
	 case 1:
	      if (p->is_in_place && howmany > 1 && istride > idist)
		   rfftw_real2c_overlap_threads_aux(p->plans[0], howmany,
						    in, istride, idist,
						    out, ostride, odist,
						    (fftw_real *) work,
						    nthreads);
	      else
		   rfftw_real2c_threads_aux(p->plans[0], howmany,
					    in, istride, idist,
					    out, ostride, odist,
					    (fftw_real *) work, nthreads);
	      break;
	 default:		/* rank >= 2 */
	      {
		   if (howmany > 1 && ostride > odist)
			rfftwnd_real2c_aux_howmany_threads(p, 0, howmany,
							   in, istride, idist,
							   out, ostride, odist,
							   work, nwork,
							   nthreads);
		   else {
			int i;

			for (i = 0; i < howmany; ++i)
			     rfftwnd_real2c_threads_aux(p, 0,
							in + i * idist,
							istride,
							out + i * odist,
							ostride,
							work,
							nthreads);
		   }
	      }
     }

     fftw_free(work);
}

void rfftwnd_threads_complex_to_real(int nthreads, fftwnd_plan p, int howmany,
				     fftw_complex *in, int istride, int idist,
				     fftw_real *out, int ostride, int odist)
{
     fftw_complex *work = 0;
     int rank = p->rank;
     int nwork = p->nwork, size_work = nwork * nthreads;

     if (p->dir != FFTW_COMPLEX_TO_REAL)
	  fftw_die("rfftwnd_complex_to_real with real-to-complex plan");

     if (p->is_in_place) {
	  ostride = istride;
	  odist = idist;
	  odist = (idist == 1) ? 1 : (idist * 2);	/* ugh */
	  out = (fftw_real *) in;
	  if (howmany > 1 && istride > idist && rank > 0) {
	       int new_nwork = p->n[rank - 1] * howmany;
	       if (new_nwork > nwork)
		    nwork = new_nwork;
	       if (rank != 1) {
		    if (nwork * nthreads > size_work)
			 size_work = nwork * nthreads;
	       }
	       else
		    size_work = nwork;
	  }
     }

     work = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * size_work);

     switch (rank) {
	 case 0:
	      break;
	 case 1:
              if (p->is_in_place && howmany > 1 && istride > idist)
                   rfftw_c2real_overlap_threads_aux(p->plans[0], howmany,
                                                    in, istride, idist,
                                                    out, ostride, odist,
                                                    (fftw_real *) work,
                                                    nthreads);
              else
                   rfftw_c2real_threads_aux(p->plans[0], howmany,
                                            in, istride, idist,
                                            out, ostride, odist,
                                            (fftw_real *) work, nthreads);
	      break;
	 default:		/* rank >= 2 */
	      {
		   if (howmany > 1 && ostride > odist)
                       rfftwnd_c2real_aux_howmany_threads(p, 0, howmany,
                                                           in, istride, idist,
                                                           out, ostride, odist,
                                                           work, nwork,
                                                           nthreads);
		   else {
			int i;

			for (i = 0; i < howmany; ++i)
                             rfftwnd_c2real_threads_aux(p, 0,
                                                        in + i * idist,
                                                        istride,
                                                        out + i * odist,
                                                        ostride,
                                                        work,
                                                        nthreads);
		   }
	      }
     }

     fftw_free(work);
}

void rfftwnd_threads_one_real_to_complex(int nthreads, fftwnd_plan p,
					 fftw_real *in, fftw_complex *out)
{
     rfftwnd_threads_real_to_complex(nthreads, p, 1, in, 1, 1, out, 1, 1);
}

void rfftwnd_threads_one_complex_to_real(int nthreads, fftwnd_plan p,
					 fftw_complex *in, fftw_real *out)
{
     rfftwnd_threads_complex_to_real(nthreads, p, 1, in, 1, 1, out, 1, 1);
}
