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
 * executor_threads.c -- execute the fft in parallel using threads
 */

#include <stdio.h>
#include <stdlib.h>

#include <fftw_threads-int.h>

static void executor_simple_threads(int n, const fftw_complex *in,
				    fftw_complex *out,
				    fftw_plan_node *p,
				    int istride,
				    int ostride,
				    int nthreads);

typedef struct {
     int m,r;
     const fftw_complex *in;
     fftw_complex *out;
     fftw_plan_node *p;
     int istride, ostride;
     int nthreads;
} executor_simple_data;

static void *executor_simple_thread(fftw_loop_data *loop_data)
{
     int min = loop_data->min, max = loop_data->max;
     executor_simple_data *d = (executor_simple_data *) loop_data->data;
     int m = d->m, r = d->r;
     const fftw_complex *in = d->in;
     fftw_complex *out = d->out;
     fftw_plan_node *p = d->p;
     int istride = d->istride, ostride = d->ostride;
     int nthreads = d->nthreads;
     
     for (; min < max; ++min)
	  executor_simple_threads(m, in + min * istride,
				  out + min * (m * ostride),
				  p,
				  istride * r, ostride,
				  nthreads);

     return 0;
}

typedef struct {
     fftw_twiddle_codelet *codelet;
     int m, ntwiddle, ostride;
     fftw_complex *out, *W;
} twiddle_thread_data;

static void *twiddle_thread(fftw_loop_data *loop_data)
{
     twiddle_thread_data *d = (twiddle_thread_data *) loop_data->data;     
     HACK_ALIGN_STACK_EVEN();
     (d->codelet)(d->out + d->ostride * loop_data->min,
		  d->W + d->ntwiddle * loop_data->min,
		  d->m * d->ostride,
		  loop_data->max - loop_data->min,
		  d->ostride);
     return 0;
}

static void executor_simple_threads(int n, const fftw_complex *in,
				    fftw_complex *out,
				    fftw_plan_node *p,
				    int istride,
				    int ostride,
				    int nthreads)
{
     switch (p->type) {
	 case FFTW_NOTW:
              HACK_ALIGN_STACK_ODD();
	      (p->nodeu.notw.codelet) (in, out, istride, ostride);
	      break;

	 case FFTW_TWIDDLE:
	      {
		   int r = p->nodeu.twiddle.size;
		   int m = n / r;
		   int i;

		   if (nthreads <= 1) {
			fftw_twiddle_codelet *codelet;
			fftw_complex *W;

			for (i = 0; i < r; ++i) {
			     fftw_executor_simple(m, in + i * istride,
						  out + i * (m * ostride),
						  p->nodeu.twiddle.recurse,
						  istride * r, ostride);
			}
			codelet = p->nodeu.twiddle.codelet;
			W = p->nodeu.twiddle.tw->twarray;
			HACK_ALIGN_STACK_EVEN();
			codelet(out, W, m * ostride, m, ostride);
		   }
		   else {
			{
			     executor_simple_data d;
			     
			     d.m = m; d.r = r;
			     d.in = in; d.out = out;
			     d.p = p->nodeu.twiddle.recurse;
			     d.istride = istride;
			     d.ostride = ostride;
			     d.nthreads = nthreads / r;
			     
			     fftw_thread_spawn_loop(r, nthreads,
						    executor_simple_thread,&d);
			}
			{
			     twiddle_thread_data d;

			     d.codelet = p->nodeu.twiddle.codelet;
			     d.m = m;
			     d.ntwiddle =
				  p->nodeu.twiddle.codelet_desc->ntwiddle;
			     d.ostride = ostride;
			     d.out = out;
			     d.W = p->nodeu.twiddle.tw->twarray;

			     fftw_thread_spawn_loop(m, nthreads,
						    twiddle_thread, &d);
			}
		   }

		   break;
	      }

	 case FFTW_RADER:
	      {
		   int r = p->nodeu.twiddle.size;
		   int m = n / r;
		   int i;

		   if (nthreads <= 1) {
			fftw_rader_codelet *codelet;
			fftw_complex *W;

			for (i = 0; i < r; ++i) {
			     fftw_executor_simple(m, in + i * istride,
						  out + i * (m * ostride),
						  p->nodeu.rader.recurse,
						  istride * r, ostride);
			}
			codelet = p->nodeu.rader.codelet;
			W = p->nodeu.rader.tw->twarray;
			codelet(out, W, m, r, ostride,
				p->nodeu.rader.rader_data);
		   }
		   else {
			{
			     executor_simple_data d;

			     d.m = m; d.r = r;
			     d.in = in; d.out = out;
			     d.p = p->nodeu.rader.recurse;
			     d.istride = istride;
			     d.ostride = ostride;
			     d.nthreads = nthreads / r;
			     
			     fftw_thread_spawn_loop(r, nthreads,
						    executor_simple_thread,&d);
			}
			{
			     fftw_rader_codelet *codelet;
			     fftw_complex *W;
			     
			     codelet = p->nodeu.rader.codelet;
			     W = p->nodeu.rader.tw->twarray;
			     codelet(out, W, m, r, ostride,
				     p->nodeu.rader.rader_data);
			}
		   }

		   break;
	      }

      	 case FFTW_GENERIC:
	      {
		   int r = p->nodeu.generic.size;
		   int m = n / r;
		   int i;
		   fftw_generic_codelet *codelet;
		   fftw_complex *W;

		   if (nthreads <= 1)
			for (i = 0; i < r; ++i) {
			     fftw_executor_simple(m, in + i * istride,
						  out + i * (m * ostride),
						  p->nodeu.generic.recurse,
						  istride * r, ostride);
			}
		   else {
			executor_simple_data d;

			d.m = m; d.r = r;
			d.in = in; d.out = out;
			d.p = p->nodeu.generic.recurse;
			d.istride = istride;
			d.ostride = ostride;
			d.nthreads = nthreads / r;

			fftw_thread_spawn_loop(r, nthreads,
					       executor_simple_thread, &d);
		   }

		   codelet = p->nodeu.generic.codelet;
		   W = p->nodeu.generic.tw->twarray;
		   codelet(out, W, m, r, n, ostride);

		   break;
	      }

	 default:
	      fftw_die("BUG in executor: invalid plan\n");
	      break;
     }
}

static void executor_simple_inplace_threads(int n, fftw_complex *in,
					    fftw_complex *out,
					    fftw_plan_node *p,
					    int istride,
					    int nthreads)
{
     switch (p->type) {
	 case FFTW_NOTW:
              HACK_ALIGN_STACK_ODD();
	      (p->nodeu.notw.codelet) (in, in, istride, istride);
	      break;

	 default:
	      {
		   fftw_complex *tmp;
		   
		   if (out)
			tmp = out;
		   else
			tmp = (fftw_complex *)
			     fftw_malloc(n * sizeof(fftw_complex));
		   
		   executor_simple_threads(n, in, tmp, p, istride, 1, 
					   nthreads);
		   fftw_strided_copy(n, tmp, istride, in);

		   if (!out)
			fftw_free(tmp);
	      }
     }
}

typedef struct {
     union {
	  fftw_notw_codelet *codelet;
	  struct {
	       int n;
	       fftw_plan_node *p;
	  } plan;
     } u;
     const fftw_complex *in;
     fftw_complex *out;
     int idist, odist, istride, ostride;
} executor_many_data;

static void *executor_many_codelet_thread(fftw_loop_data *loop_data)
{
     int min = loop_data->min, max = loop_data->max;
     executor_many_data *d = (executor_many_data *) loop_data->data;
     fftw_notw_codelet *codelet = d->u.codelet;
     const fftw_complex *in = d->in;
     fftw_complex *out = d->out;
     int idist = d->idist, odist = d->odist;
     int istride = d->istride, ostride = d->ostride;
     
     HACK_ALIGN_STACK_ODD();
     for (; min < max; ++min)
	  codelet(in + min * idist,
		  out + min * odist,
		  istride, ostride);

     return 0;
}

static void *executor_many_simple_thread(fftw_loop_data *loop_data)
{
     int min = loop_data->min, max = loop_data->max;
     executor_many_data *d = (executor_many_data *) loop_data->data;
     int n = d->u.plan.n;
     fftw_plan_node *p = d->u.plan.p;
     const fftw_complex *in = d->in;
     fftw_complex *out = d->out;
     int idist = d->idist, odist = d->odist;
     int istride = d->istride, ostride = d->ostride;

     for (; min < max; ++min)
	  fftw_executor_simple(n, in + min * idist,
			       out + min * odist,
			       p, istride, ostride);

     return 0;
}

static void executor_many_threads(int n, const fftw_complex *in,
				  fftw_complex *out,
				  fftw_plan_node *p,
				  int istride,
				  int ostride,
				  int howmany, int idist, int odist,
				  int nthreads)
{
     switch (p->type) {
	 case FFTW_NOTW:
	      {
		   int s;
		   
		   if (nthreads <= 1) {
			fftw_notw_codelet *codelet = p->nodeu.notw.codelet;
			HACK_ALIGN_STACK_ODD();
			for (s = 0; s < howmany; ++s)
			     codelet(in + s * idist,
				     out + s * odist,
				     istride, ostride);
		   }
		   else {
			executor_many_data d;

			d.in = in;
			d.out = out;
			d.u.codelet = p->nodeu.notw.codelet;
			d.istride = istride;
			d.ostride = ostride;
			d.idist = idist;
			d.odist = odist;
			fftw_thread_spawn_loop(howmany, nthreads,
				   executor_many_codelet_thread, &d);
		   }

		   break;
	      }

	 default:
	      {
		   int s;

		   if (nthreads <= 1)
			for (s = 0; s < howmany; ++s) {
			     fftw_executor_simple(n, in + s * idist,
						  out + s * odist,
						  p, istride, ostride);
			}
		   else {
			executor_many_data d;

			d.in = in; d.out = out;
			d.u.plan.n = n;
			d.u.plan.p = p;
			d.istride = istride;
			d.ostride = ostride;
			d.idist = idist;
			d.odist = odist;
			fftw_thread_spawn_loop(howmany, nthreads,
				   executor_many_simple_thread, &d);
		   }
	      }
     }
}

typedef struct {
     union {
	  fftw_notw_codelet *codelet;
	  struct {
	       int n;
	       fftw_plan_node *p;
	       fftw_complex *tmp;
	  } plan;
     } u;
     fftw_complex *in;
     int idist, istride;
} executor_many_inplace_data;

static void *executor_many_inplace_codelet_thread(fftw_loop_data *loop_data)
{
     int min = loop_data->min, max = loop_data->max;
     executor_many_inplace_data
	  *d = (executor_many_inplace_data *) loop_data->data;
     fftw_notw_codelet *codelet = d->u.codelet;
     fftw_complex *in = d->in;
     int idist = d->idist, istride = d->istride;
     
     HACK_ALIGN_STACK_ODD();
     for (; min < max; ++min)
	  codelet(in + min * idist,
		  in + min * idist,
		  istride, istride);

     return 0;
}

static void *executor_many_inplace_simple_thread(fftw_loop_data *loop_data)
{
     int min = loop_data->min, max = loop_data->max;
     executor_many_inplace_data
          *d = (executor_many_inplace_data *) loop_data->data;
     int n = d->u.plan.n;
     fftw_plan_node *p = d->u.plan.p;
     fftw_complex *tmp = d->u.plan.tmp + n * loop_data->thread_num;
     fftw_complex *in = d->in;
     int idist = d->idist, istride = d->istride;

     for (; min < max; ++min) {
	  fftw_executor_simple(n, in + min * idist,
			       tmp,
			       p, istride, 1);
	  fftw_strided_copy(n, tmp, istride, in + min * idist);
     }

     return 0;
}

void fftw_executor_many_inplace_threads(int n, fftw_complex *in,
					fftw_complex *work,
					fftw_plan_node *p,
					int istride,
					int howmany, int idist,
					int nthreads)
{
     switch (p->type) {
	 case FFTW_NOTW:
	      {
		   int s;

		   if (nthreads <= 1) {
			fftw_notw_codelet *codelet = p->nodeu.notw.codelet;
			HACK_ALIGN_STACK_ODD();
			for (s = 0; s < howmany; ++s)
			     codelet(in + s * idist,
				     in + s * idist,
				     istride, istride);
		   }
		   else {
			executor_many_inplace_data d;
			
			d.in = in;
			d.u.codelet = p->nodeu.notw.codelet;
			d.istride = istride;
			d.idist = idist;

			fftw_thread_spawn_loop(howmany, nthreads,
				   executor_many_inplace_codelet_thread, &d);
		   }
		   break;
	      }

	 default:
	      {
		   int s;
		   fftw_complex *tmp;

		   if (nthreads <= 1) {
			if (work)
			     tmp = work;
			else
			     tmp = (fftw_complex *)
				  fftw_malloc(n * sizeof(fftw_complex));

			for (s = 0; s < howmany; ++s) {
			     fftw_executor_simple(n,
						  in + s * idist,
						  tmp,
						  p, istride, 1);
			     fftw_strided_copy(n, tmp, istride, 
					       in + s * idist);
			}
		   }
		   else {
			executor_many_inplace_data d;
			
			if (work)
			     tmp = work;
			else
			     tmp = (fftw_complex *)
				  fftw_malloc((nthreads > howmany ? 
					       howmany : nthreads) 
					      * n * sizeof(fftw_complex));
			
			d.in = in;
			d.u.plan.n = n;
			d.u.plan.p = p;
			d.u.plan.tmp = tmp;
			d.istride = istride;
			d.idist = idist;
			fftw_thread_spawn_loop(howmany, nthreads,
				   executor_many_inplace_simple_thread, &d);
		   }

		   if (!work)
			fftw_free(tmp);
	      }
     }
}

/* user interface */
void fftw_threads(int nthreads,
		  fftw_plan plan, int howmany, fftw_complex *in, int istride,
		  int idist, fftw_complex *out, int ostride, int odist)
{
     int n = plan->n;

     if (plan->flags & FFTW_IN_PLACE) {
	  if (howmany == 1) {
	       executor_simple_inplace_threads(n, in, out, plan->root, 
					       istride, nthreads);
	  } else {
	       fftw_executor_many_inplace_threads(n, in, NULL,
						  plan->root, istride, 
						  howmany, idist, nthreads);
	  }
     } else {
	  if (howmany == 1) {
	       executor_simple_threads(n, in, out, plan->root, 
				       istride, ostride, nthreads);
	  } else {
	       executor_many_threads(n, in, out, plan->root, istride, ostride,
				     howmany, idist, odist, nthreads);
	  }
     }
}

void fftw_threads_one(int nthreads,
		      fftw_plan plan, fftw_complex *in, fftw_complex *out)
{
     if (plan->flags & FFTW_IN_PLACE)
	  executor_simple_inplace_threads(plan->n, in, out, plan->root,
					  1, nthreads);
     else
	  executor_simple_threads(plan->n, in, out, plan->root,
				  1, 1, nthreads);
}
