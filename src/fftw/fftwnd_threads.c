
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

#include <stdlib.h>

#include <fftw_threads-int.h>
#include <fftw_threads.h>

typedef struct {
     fftwnd_plan plan;
     int cur_dim;
     int distance;
     fftw_complex *in, *out;
     int istride, ostride;
     fftw_complex *work;
} fftwnd_aux_many_data;

static void *fftwnd_aux_many_thread(fftw_loop_data *loop_data)
{
     int min = loop_data->min, max = loop_data->max;
     fftwnd_aux_many_data *d = (fftwnd_aux_many_data *) loop_data->data;
     int distance = d->distance, cur_dim = d->cur_dim;
     fftwnd_plan plan = d->plan;
     fftw_complex *in = d->in, *out = d->out;
     int istride = d->istride, ostride = d->ostride;
     fftw_complex *work = d->work + loop_data->thread_num * plan->nwork;

     for (; min < max; ++min)
	  fftwnd_aux(plan,cur_dim,
		     in + min*istride*distance,istride,
		     out + min*ostride*distance,ostride,
		     work);

     return 0;
}

static void fftwnd_aux_many_threads(int nthreads, int n, int n_after,
				    fftwnd_plan plan, int cur_dim,
				    fftw_complex *in, int istride,
				    fftw_complex *out, int ostride)
{
     fftw_complex *tmp;
     fftwnd_aux_many_data d;

     if (nthreads > n)
	  nthreads = n;

     tmp = (fftw_complex *) fftw_malloc(nthreads * plan->nwork
					* sizeof(fftw_complex));
     
     d.plan = plan;
     d.cur_dim = cur_dim;
     d.distance = n_after;
     d.in = in;
     d.out = out;
     d.istride = istride;
     d.ostride = ostride;
     d.work = tmp;
     
     fftw_thread_spawn_loop(n, nthreads, fftwnd_aux_many_thread, &d);

     fftw_free(tmp);
}

static void fftwnd_threads_aux(int nthreads, fftwnd_plan p, int cur_dim, 
			       fftw_complex *in, int istride,
			       fftw_complex *out, int ostride)
{
     int n_after = p->n_after[cur_dim], n = p->n[cur_dim];

     if (cur_dim == p->rank - 2) {
	  /* just do the last dimension directly: */
	  if (p->is_in_place)
	       fftw_threads(nthreads, p->plans[p->rank - 1], n,
			    in, istride, n_after * istride,
			    (fftw_complex*)NULL, 0, 0);
	  else
	       fftw_threads(nthreads, p->plans[p->rank - 1], n,
			    in,  istride, n_after * istride,
			    out, ostride, n_after * ostride);
     }
     else { /* we have at least two dimensions to go */
	  /* process the subsequent dimensions recursively, in hyperslabs,
	     to get maximum locality: */
	  fftwnd_aux_many_threads(nthreads, n, n_after,
				  p, cur_dim + 1,
				  in, istride, out, ostride);
     }

     /* do the current dimension (in-place): */
     fftw_threads(nthreads, p->plans[cur_dim], n_after,
		  out, n_after * ostride, ostride,
		  (fftw_complex*)NULL, 0, 0);
}

void fftwnd_threads(int nthreads, fftwnd_plan p, int howmany,
		    fftw_complex *in, int istride, int idist,
		    fftw_complex *out, int ostride, int odist)
{
     switch (p->rank) {
	 case 0:
	      break;
	 case 1:
	      if (p->is_in_place)	/* fft is in-place */
		   fftw_threads(nthreads, p->plans[0], howmany,
				in, istride, idist,
				(fftw_complex*) NULL, 0, 0);
	      else
		   fftw_threads(nthreads, p->plans[0], howmany,
				in, istride, idist,
				out, ostride, odist);
	      break;
	 default: /* rank >= 2 */
	 {
	      int i;
	      
	      if (p->is_in_place) {
		   out = in;
		   ostride = istride;
		   odist = idist;
	      }

	      if (nthreads <= 1)
		   fftwnd(p, howmany, in, istride, idist, out, ostride, odist);
	      else
		   for (i = 0; i < howmany; ++i)
			fftwnd_threads_aux(nthreads, p, 0,
					   in + i*idist, istride,
					   out + i*odist, ostride);
	 }
     }
}

void fftwnd_threads_one(int nthreads, fftwnd_plan p,
			fftw_complex *in, fftw_complex *out)
{
     fftwnd_threads(nthreads, p, 1, in, 1, 0, out, 1, 0);
}
