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

/********************** prototypes for rexec2 routines **********************/

extern void rfftw_hc2c(int n, fftw_real *in, fftw_complex *out, int ostride);
extern void rfftw_c2hc(int n, fftw_complex *in, int istride, fftw_real *out);
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


/***************************************************************************/

typedef struct {
     fftw_plan plan;
     void *in;
     int istride, idist;
     void *out;
     int ostride, odist;
     fftw_real *work;
} rexec2_thread_data;

static void *real2c_aux_thread(fftw_loop_data *ldata)
{
     rexec2_thread_data *d = (rexec2_thread_data *) ldata->data;
     rfftw_real2c_aux(d->plan, ldata->max - ldata->min,
		      ldata->min * d->idist + (fftw_real *) d->in,
		      d->istride, d->idist,
		      ldata->min * d->odist + (fftw_complex *) d->out,
		      d->ostride, d->odist,
		      d->work + d->plan->n * ldata->thread_num);
     return 0;
}

static void *c2real_aux_thread(fftw_loop_data *ldata)
{
     rexec2_thread_data *d = (rexec2_thread_data *) ldata->data;
     rfftw_c2real_aux(d->plan, ldata->max - ldata->min,
		      ldata->min * d->idist + (fftw_complex *) d->in,
		      d->istride, d->idist,
		      ldata->min * d->odist + (fftw_real *) d->out,
		      d->ostride, d->odist,
		      d->work + d->plan->n * ldata->thread_num);
     return 0;
}

static void *real2c_overlap_aux_thread1(fftw_loop_data *ldata)
{
     rexec2_thread_data *d = (rexec2_thread_data *) ldata->data;
     int n = d->plan->n;
     rfftw(d->plan, ldata->max - ldata->min,
	   ldata->min * d->idist + (fftw_real *) d->in,
	   d->istride, d->idist,
	   d->work + n * ldata->min, 1, n);
     return 0;
}

static void *real2c_overlap_aux_thread2(fftw_loop_data *ldata)
{
     rexec2_thread_data *d = (rexec2_thread_data *) ldata->data;
     int min = ldata->min, max = ldata->max;
     int n = d->plan->n;
     fftw_complex *out = (fftw_complex *) d->out;
     int ostride = d->ostride, odist = d->odist;
     fftw_real *work = d->work;

     for (; min < max; ++min)
	  rfftw_hc2c(n, work + min*n, out + min*odist, ostride);
     return 0;
}

static void *c2real_overlap_aux_thread2(fftw_loop_data *ldata)
{
     rexec2_thread_data *d = (rexec2_thread_data *) ldata->data;
     int n = d->plan->n;
     rfftw(d->plan, ldata->max - ldata->min,
	   d->work + n * ldata->min, 1, n,
	   ldata->min * d->odist + (fftw_real *) d->out,
	   d->ostride, d->odist);
     return 0;
}

static void *c2real_overlap_aux_thread1(fftw_loop_data *ldata)
{
     rexec2_thread_data *d = (rexec2_thread_data *) ldata->data;
     int min = ldata->min, max = ldata->max;
     int n = d->plan->n;
     fftw_complex *in = (fftw_complex *) d->in;
     int istride = d->istride, idist = d->idist;
     fftw_real *work = d->work;

     for (; min < max; ++min)
	  rfftw_c2hc(n, in + min*idist, istride, work + min*n);
     return 0;
}

void rfftw_real2c_threads_aux(fftw_plan plan, int howmany,
			      fftw_real *in, int istride, int idist,
			      fftw_complex *out, int ostride, int odist,
			      fftw_real *work,
			      int nthreads)
{
     rexec2_thread_data d;
     
     d.plan = plan;
     d.in = in;
     d.istride = istride;
     d.idist = idist;
     d.out = out;
     d.ostride = ostride;
     d.odist = odist;
     d.work = work;
     
     fftw_thread_spawn_loop(howmany, nthreads, real2c_aux_thread, &d);
}

void rfftw_c2real_threads_aux(fftw_plan plan, int howmany,
			      fftw_complex *in, int istride, int idist,
			      fftw_real *out, int ostride, int odist,
			      fftw_real *work,
			      int nthreads)
{
     rexec2_thread_data d;
     
     d.plan = plan;
     d.in = in;
     d.istride = istride;
     d.idist = idist;
     d.out = out;
     d.ostride = ostride;
     d.odist = odist;
     d.work = work;
     
     fftw_thread_spawn_loop(howmany, nthreads, c2real_aux_thread, &d);
}

void rfftw_real2c_overlap_threads_aux(fftw_plan plan, int howmany,
			      fftw_real *in, int istride, int idist,
			      fftw_complex *out, int ostride, int odist,
			      fftw_real *work,
			      int nthreads)
{
     rexec2_thread_data d;
     
     d.plan = plan;
     d.in = in;
     d.istride = istride;
     d.idist = idist;
     d.out = out;
     d.ostride = ostride;
     d.odist = odist;
     d.work = work;

     fftw_thread_spawn_loop(howmany, nthreads, real2c_overlap_aux_thread1, &d);
     fftw_thread_spawn_loop(howmany, nthreads, real2c_overlap_aux_thread2, &d);
}

void rfftw_c2real_overlap_threads_aux(fftw_plan plan, int howmany,
			      fftw_complex *in, int istride, int idist,
			      fftw_real *out, int ostride, int odist,
			      fftw_real *work,
			      int nthreads)
{
     rexec2_thread_data d;
     
     d.plan = plan;
     d.in = in;
     d.istride = istride;
     d.idist = idist;
     d.out = out;
     d.ostride = ostride;
     d.odist = odist;
     d.work = work;

     fftw_thread_spawn_loop(howmany, nthreads, c2real_overlap_aux_thread1, &d);
     fftw_thread_spawn_loop(howmany, nthreads, c2real_overlap_aux_thread2, &d);
}
