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
 * rexec_threads.c -- execute the fft in parallel
 */

#include <stdio.h>
#include <stdlib.h>

#include <fftw_threads-int.h>
#include <rfftw_threads.h>

extern void rfftw_strided_copy(int n, fftw_real *in, int ostride,
			       fftw_real *out);
extern void rfftw_executor_simple(int n, fftw_real *in,
				  fftw_real *out,
				  fftw_plan_node *p,
				  int istride,
				  int ostride);

static void rexec_simple_threads(int n, fftw_real *in,
				 fftw_real *out,
				 fftw_plan_node *p,
				 int istride,
				 int ostride,
				 int nthreads);

typedef struct {
     int m,r;
     fftw_real *in;
     fftw_real *out;
     fftw_plan_node *p;
     int istride, ostride;
     int nthreads;
} rexec_simple_data;

static void *rexec_simple_thread_r2c(fftw_loop_data *ldata)
{
     int min = ldata->min, max = ldata->max;
     rexec_simple_data *d = (rexec_simple_data *) ldata->data;
     int m = d->m, r = d->r;
     fftw_real *in = d->in;
     fftw_real *out = d->out;
     fftw_plan_node *p = d->p;
     int istride = d->istride, ostride = d->ostride;
     int nthreads = d->nthreads;
     
     for (; min < max; ++min)
	  rexec_simple_threads(m, in + min * istride,
			       out + min * (m * ostride),
			       p,
			       istride * r, ostride,
			       nthreads);

     return 0;
}

static void *rexec_simple_thread_c2r(fftw_loop_data *ldata)
{
     int min = ldata->min, max = ldata->max;
     rexec_simple_data *d = (rexec_simple_data *) ldata->data;
     int m = d->m, r = d->r;
     fftw_real *in = d->in;
     fftw_real *out = d->out;
     fftw_plan_node *p = d->p;
     int istride = d->istride, ostride = d->ostride;
     int nthreads = d->nthreads;
     
     for (; min < max; ++min)
	  rexec_simple_threads(m, in + min * (m * istride),
			       out + min * ostride,
			       p,
			       istride, ostride * r,
			       nthreads);

     return 0;
}

static void spawn_h2hc_recurse_threads(int m, int r,
				       fftw_real *in,
				       fftw_real *out,
				       fftw_plan_node *p,
				       int istride,
				       int ostride,
				       int nthreads)
{
     rexec_simple_data d;
     
     d.m = m; d.r = r;
     d.in = in; d.out = out;
     d.p = p->nodeu.hc2hc.recurse;
     d.istride = istride;
     d.ostride = ostride;
     d.nthreads = nthreads / r;
     
     switch (p->nodeu.hc2hc.dir) {
	 case FFTW_REAL_TO_COMPLEX:
	      fftw_thread_spawn_loop(r, nthreads,
				     rexec_simple_thread_r2c, &d);
	      break;
	 case FFTW_COMPLEX_TO_REAL:
	      fftw_thread_spawn_loop(r, nthreads,
				     rexec_simple_thread_c2r, &d);
	      break;
     }
}

static void rexec_simple_threads(int n, fftw_real *in, fftw_real *out,
				 fftw_plan_node *p,
				 int istride,
				 int ostride,
				 int nthreads)
{
     switch (p->type) {
	 case FFTW_REAL2HC:
	      HACK_ALIGN_STACK_ODD();
	      (p->nodeu.real2hc.codelet) (in, out, out + n * ostride,
					  istride, ostride, -ostride);
	      break;

	 case FFTW_HC2REAL:
	      HACK_ALIGN_STACK_ODD();
	      (p->nodeu.hc2real.codelet) (in, in + n * istride, out,
					  istride, -istride, ostride);
	      break;

	 case FFTW_HC2HC:
	      {
		   int r = p->nodeu.hc2hc.size;
		   int m = n / r;
		   int i;
		   fftw_hc2hc_codelet *codelet;
		   fftw_complex *W;

		   if (nthreads <= 1) { 
			switch (p->nodeu.hc2hc.dir) {
			    case FFTW_REAL_TO_COMPLEX:
				 for (i = 0; i < r; ++i)
				      rfftw_executor_simple(m,
						    in + i * istride,
						    out + i * (m * ostride),
						    p->nodeu.hc2hc.recurse,
						    istride * r, ostride);

				 W = p->nodeu.hc2hc.tw->twarray;
				 codelet = p->nodeu.hc2hc.codelet;
				 HACK_ALIGN_STACK_EVEN();
				 codelet(out, W, m * ostride, m, ostride);
				 break;
			    case FFTW_COMPLEX_TO_REAL:
				 W = p->nodeu.hc2hc.tw->twarray;
				 codelet = p->nodeu.hc2hc.codelet;
				 HACK_ALIGN_STACK_EVEN();
				 codelet(in, W, m * istride, m, istride);
				 
				 for (i = 0; i < r; ++i)
				      rfftw_executor_simple(m,
						    in + i * (m * istride),
						    out + i * ostride,
						    p->nodeu.hc2hc.recurse,
						    istride, ostride * r);
				 break;
			    default:
				 goto bug;
			}
		   }
		   else
			switch (p->nodeu.hc2hc.dir) {
			    case FFTW_REAL_TO_COMPLEX:
				 spawn_h2hc_recurse_threads(m, r, in, out, p,
							    istride, ostride,
							    nthreads);

				 W = p->nodeu.hc2hc.tw->twarray;
				 codelet = p->nodeu.hc2hc.codelet;
				 HACK_ALIGN_STACK_EVEN();
				 codelet(out, W, m * ostride, m, ostride);

				 break;
			    case FFTW_COMPLEX_TO_REAL:
				 W = p->nodeu.hc2hc.tw->twarray;
				 codelet = p->nodeu.hc2hc.codelet;
				 HACK_ALIGN_STACK_EVEN();
				 codelet(in, W, m * istride, m, istride);
				 
				 spawn_h2hc_recurse_threads(m, r, in, out, p,
							    istride, ostride,
							    nthreads);
				 break;
			}
		   
		   break;
	      }

	 case FFTW_RGENERIC:
	      {
		   int r = p->nodeu.rgeneric.size;
		   int m = n / r;
		   int i;
		   fftw_rgeneric_codelet *codelet = p->nodeu.rgeneric.codelet;
		   fftw_complex *W = p->nodeu.rgeneric.tw->twarray;

		   if (nthreads <= 1)
			switch (p->nodeu.rgeneric.dir) {
			    case FFTW_REAL_TO_COMPLEX:
				 for (i = 0; i < r; ++i)
				      rfftw_executor_simple(m,
						    in + i * istride,
						    out + i * (m * ostride),
						    p->nodeu.rgeneric.recurse,
						    istride * r, ostride);
				 
				 codelet(out, W, m, r, n, ostride);
				 break;
			    case FFTW_COMPLEX_TO_REAL:
				 codelet(in, W, m, r, n, istride);
				 
				 for (i = 0; i < r; ++i)
				      rfftw_executor_simple(m,
						    in + i * m * istride,
						    out + i * ostride,
						    p->nodeu.rgeneric.recurse,
						    istride, ostride * r);
				 break;
			    default:
				 goto bug;
			}
		   else
			switch (p->nodeu.hc2hc.dir) {
			    case FFTW_REAL_TO_COMPLEX:
				 spawn_h2hc_recurse_threads(m, r, in, out, p,
							    istride, ostride,
							    nthreads);
				 codelet(out, W, m, r, n, ostride);
				 break;
			    case FFTW_COMPLEX_TO_REAL:
				 codelet(in, W, m, r, n, istride);
				 spawn_h2hc_recurse_threads(m, r, in, out, p,
							    istride, ostride,
							    nthreads);
				 break;
			}

		   break;
	      }

	 default:
	    bug:
	      fftw_die("BUG in rexecutor: invalid plan\n");
	      break;
     }
}

static void rexecutor_simple_inplace_threads(int n, fftw_real *in,
					     fftw_real *out,
					     fftw_plan_node *p,
					     int istride,
					     int nthreads)
{
     switch (p->type) {
	 case FFTW_REAL2HC:
	      HACK_ALIGN_STACK_ODD();
	      (p->nodeu.real2hc.codelet) (in, in, in + n * istride,
					  istride, istride, -istride);
	      break;

	 case FFTW_HC2REAL:
	      HACK_ALIGN_STACK_ODD();
	      (p->nodeu.hc2real.codelet) (in, in + n * istride, in,
					  istride, -istride, istride);
	      break;

	 default:
	      {
		   fftw_real *tmp;

		   if (out)
			tmp = out;
		   else
			tmp = (fftw_real *) fftw_malloc(n * sizeof(fftw_real));

		   rexec_simple_threads(n, in, tmp, p, istride, 1, nthreads);
		   rfftw_strided_copy(n, tmp, istride, in);

		   if (!out)
			fftw_free(tmp);
	      }
     }
}

typedef struct {
     union {
          fftw_real2hc_codelet *r2c_codelet;
          fftw_hc2real_codelet *c2r_codelet;
	  fftw_plan_node *p;
     } u;
     int n;
     fftw_real *in;
     fftw_real *out;
     int idist, odist, istride, ostride;
} rexec_many_data;

static void *rexec_many_r2c_codelet_thread(fftw_loop_data *ldata)
{
     int min = ldata->min, max = ldata->max;
     rexec_many_data *d = (rexec_many_data *) ldata->data;
     fftw_real2hc_codelet *r2c_codelet = d->u.r2c_codelet;
     int n = d->n;
     fftw_real *in = d->in;
     fftw_real *out = d->out;
     int idist = d->idist, odist = d->odist;
     int istride = d->istride, ostride = d->ostride;

     HACK_ALIGN_STACK_ODD();
     for (; min < max; ++min)
          r2c_codelet(in + min * idist,
		      out + min * odist,
		      out + n * ostride + min * odist,
		      istride, ostride, -ostride);
     return 0;
}

static void *rexec_many_c2r_codelet_thread(fftw_loop_data *ldata)
{
     int min = ldata->min, max = ldata->max;
     rexec_many_data *d = (rexec_many_data *) ldata->data;
     fftw_hc2real_codelet *c2r_codelet = d->u.c2r_codelet;
     int n = d->n;
     fftw_real *in = d->in;
     fftw_real *out = d->out;
     int idist = d->idist, odist = d->odist;
     int istride = d->istride, ostride = d->ostride;

     HACK_ALIGN_STACK_ODD();
     for (; min < max; ++min)
          c2r_codelet(in + min * idist,
		      in + n * istride + min * idist,
		      out + min * odist,
		      istride, -istride, ostride);
     return 0;
}

static void *rexec_many_simple_thread(fftw_loop_data *ldata)
{
     int min = ldata->min, max = ldata->max;
     rexec_many_data *d = (rexec_many_data *) ldata->data;
     fftw_plan_node *p = d->u.p;
     int n = d->n;
     fftw_real *in = d->in;
     fftw_real *out = d->out;
     int idist = d->idist, odist = d->odist;
     int istride = d->istride, ostride = d->ostride;

     for (; min < max; ++min)
          rfftw_executor_simple(n, in + min * idist,
				out + min * odist,
				p, istride, ostride);

     return 0;
}

static void rexecutor_many_threads(int n, fftw_real *in,
				   fftw_real *out,
				   fftw_plan_node *p,
				   int istride,
				   int ostride,
				   int howmany, int idist, int odist,
				   int nthreads)
{
     if (nthreads > howmany)
	  nthreads = howmany;

     switch (p->type) {
	 case FFTW_REAL2HC:
	 {
	      int s;
	      fftw_real2hc_codelet *codelet = p->nodeu.real2hc.codelet;
	      
	      if (nthreads <= 1) {
		   HACK_ALIGN_STACK_ODD();
		   for (s = 0; s < howmany; ++s)
			codelet(in + s * idist, out + s * odist,
				out + n * ostride + s * odist,
				istride, ostride, -ostride);
	      }
	      else {
		   rexec_many_data d;
		   
		   d.n = n;
		   d.in = in;
		   d.out = out;
		   d.u.r2c_codelet = codelet;
		   d.istride = istride;
		   d.ostride = ostride;
		   d.idist = idist;
		   d.odist = odist;
		   fftw_thread_spawn_loop(howmany, nthreads,
					  rexec_many_r2c_codelet_thread, &d);
	      }
	      
	      break;
	 }
	 
	 case FFTW_HC2REAL:
	 {
	      int s;
	      fftw_hc2real_codelet *codelet = p->nodeu.hc2real.codelet;
	      
	      if (nthreads <= 1) {
		   HACK_ALIGN_STACK_ODD();
		   for (s = 0; s < howmany; ++s)
			codelet(in + s * idist,
				in + n * istride + s * idist,
				out + s * odist,
				istride, -istride, ostride);
	      }
	      else {
		   rexec_many_data d;

		   d.n = n;
		   d.in = in;
		   d.out = out;
		   d.u.c2r_codelet = codelet;
		   d.istride = istride;
		   d.ostride = ostride;
		   d.idist = idist;
		   d.odist = odist;
		   fftw_thread_spawn_loop(howmany, nthreads,
					  rexec_many_c2r_codelet_thread, &d);
	      }
	      
	      break;
	 }
	 
	 default:
	 {
	      int s;		   
	      
	      if (nthreads <= 1)
		   for (s = 0; s < howmany; ++s) {
			rfftw_executor_simple(n, in + s * idist,
					      out + s * odist,
					      p, istride, ostride);
		   }
	      else {
		   rexec_many_data d;

		   d.in = in; d.out = out;
		   d.n = n;
		   d.u.p = p;
		   d.istride = istride;
		   d.ostride = ostride;
		   d.idist = idist;
		   d.odist = odist;
		   fftw_thread_spawn_loop(howmany, nthreads,
					  rexec_many_simple_thread, &d);
	      }
	 }
     }
}

static void *rexec_many_simple_inplace_thread(fftw_loop_data *ldata)
{
     int min = ldata->min, max = ldata->max;
     rexec_many_data *d = (rexec_many_data *) ldata->data;
     fftw_plan_node *p = d->u.p;
     int n = d->n;
     fftw_real *in = d->in;
     fftw_real *out = d->out + n * ldata->thread_num;
     int idist = d->idist;
     int istride = d->istride;

     for (; min < max; ++min) {
          rfftw_executor_simple(n, in + min * idist, out, p, istride, 1);
	  rfftw_strided_copy(n, out, istride, in + min * idist);
     }

     return 0;
}

static void rexecutor_many_inplace_threads(int n, fftw_real *in,
					   fftw_real *out,
					   fftw_plan_node *p,
					   int istride,
					   int howmany, int idist,
					   int nthreads)
{
     switch (p->type) {
	 case FFTW_REAL2HC:
	 {
	      int s;
	      fftw_real2hc_codelet *codelet = p->nodeu.real2hc.codelet;
	      
	      if (nthreads <= 1) {
		   HACK_ALIGN_STACK_ODD();
		   for (s = 0; s < howmany; ++s)
			codelet(in + s * idist, in + s * idist,
				in + n * istride + s * idist,
				istride, istride, -istride);
	      }
	      else {
		   rexec_many_data d;
		   
		   d.n = n;
		   d.in = in;
		   d.out = in;
		   d.u.r2c_codelet = codelet;
		   d.istride = istride;
		   d.ostride = istride;
		   d.idist = idist;
		   d.odist = idist;
		   fftw_thread_spawn_loop(howmany, nthreads,
					  rexec_many_r2c_codelet_thread, &d);
	      }
	      
	      break;
	 }
	 
	 case FFTW_HC2REAL:
	 {
	      int s;
	      fftw_hc2real_codelet *codelet = p->nodeu.hc2real.codelet;
	      
	      if (nthreads <= 1) {
		   HACK_ALIGN_STACK_ODD();
		   for (s = 0; s < howmany; ++s)
			codelet(in + s * idist,
				in + n * istride + s * idist,
				in + s * idist,
				istride, -istride, istride);
	      }
	      else {
		   rexec_many_data d;
		   
		   d.n = n;
		   d.in = in;
		   d.out = in;
		   d.u.c2r_codelet = codelet;
		   d.istride = istride;
		   d.ostride = istride;
		   d.idist = idist;
		   d.odist = idist;
		   fftw_thread_spawn_loop(howmany, nthreads,
					  rexec_many_c2r_codelet_thread, &d);
	      }
	      
	      break;
	 }
	 
	 default:
	 {
	      int s;
	      fftw_real *tmp;
	      
	      if (nthreads > howmany)
		   nthreads = howmany;
	      
	      if (nthreads <= 1) {
		   if (out)
			tmp = out;
		   else
			tmp =(fftw_real *) fftw_malloc(n * 
						       sizeof(fftw_real));
		   
		   for (s = 0; s < howmany; ++s) {
			rfftw_executor_simple(n,
					      in + s * idist,
					      tmp,
					      p, istride, 1);
			rfftw_strided_copy(n, tmp, istride,
					   in + s * idist);
		   }
		   
		   if (!out)
			fftw_free(tmp);
	      }
	      else {
		   rexec_many_data d;
		   
		   tmp = (fftw_real *)
			fftw_malloc(nthreads * n * sizeof(fftw_real));

		   d.in = in; 
		   d.out = tmp;
		   d.n = n;
		   d.u.p = p;
		   d.istride = istride;
		   d.ostride = 1;
		   d.idist = idist;
		   d.odist = 0;
		   fftw_thread_spawn_loop(howmany, nthreads,
					  rexec_many_simple_inplace_thread,&d);
		   
		   fftw_free(tmp);
	      }
	 }
     }
}

/* user interface */
void rfftw_threads(int nthreads,
		   fftw_plan plan, int howmany, fftw_real *in, int istride,
		   int idist, fftw_real *out, int ostride, int odist)
{
     int n = plan->n;

     if (plan->flags & FFTW_IN_PLACE) {
	  if (howmany == 1) {
	       rexecutor_simple_inplace_threads(n, in, out, plan->root, 
						istride, nthreads);
	  } else {
	       rexecutor_many_inplace_threads(n, in, out, plan->root,
					      istride, howmany, idist,
					      nthreads);
	  }
     } else {
	  if (howmany == 1) {
	       rexec_simple_threads(n, in, out, plan->root, istride, ostride,
				    nthreads);
	  } else {
	       rexecutor_many_threads(n, in, out, plan->root, istride, ostride,
				      howmany, idist, odist, nthreads);
	  }
     }
}

void rfftw_threads_one(int nthreads,
		       fftw_plan plan, fftw_real *in, fftw_real *out)
{
     int n = plan->n;

     if (plan->flags & FFTW_IN_PLACE)
	  rexecutor_simple_inplace_threads(n, in, out, plan->root, 1,
					   nthreads);
     else
	  rexec_simple_threads(n, in, out, plan->root, 1, 1, nthreads);
}
