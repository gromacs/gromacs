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
 * executor.c -- execute the fft
 */

/* $Id$ */
#include <fftw-int.h>
#include <stdio.h>
#include <stdlib.h>

const char *fftw_version = "FFTW V" FFTW_VERSION " ($Id$)";

/*
 * This function is called in other files, so we cannot declare
 * it static. 
 */
void fftw_strided_copy(int n, fftw_complex *in, int ostride,
		       fftw_complex *out)
{
     int i;
     fftw_real r0, r1, i0, i1;
     fftw_real r2, r3, i2, i3;

     i = 0;

     for (; i < (n & 3); ++i) {
	  out[i * ostride] = in[i];
     }

     for (; i < n; i += 4) {
	  r0 = c_re(in[i]);
	  i0 = c_im(in[i]);
	  r1 = c_re(in[i + 1]);
	  i1 = c_im(in[i + 1]);
	  r2 = c_re(in[i + 2]);
	  i2 = c_im(in[i + 2]);
	  r3 = c_re(in[i + 3]);
	  i3 = c_im(in[i + 3]);
	  c_re(out[i * ostride]) = r0;
	  c_im(out[i * ostride]) = i0;
	  c_re(out[(i + 1) * ostride]) = r1;
	  c_im(out[(i + 1) * ostride]) = i1;
	  c_re(out[(i + 2) * ostride]) = r2;
	  c_im(out[(i + 2) * ostride]) = i2;
	  c_re(out[(i + 3) * ostride]) = r3;
	  c_im(out[(i + 3) * ostride]) = i3;
     }
}


/*
 * Do *not* declare simple executor static--we need to call it
 * from executor_cilk.cilk...also, preface its name with "fftw_"
 * to avoid any possible name collisions. 
 */
void fftw_executor_simple(int n, const fftw_complex *in,
			  fftw_complex *out,
			  fftw_plan_node *p,
			  int istride,
			  int ostride)
{
     switch (p->type) {
	 case FFTW_NOTW:
	      HACK_ALIGN_STACK_ODD();
	      (p->nodeu.notw.codelet)(in, out, istride, ostride);
	      break;

	 case FFTW_TWIDDLE:
	      {
		   int r = p->nodeu.twiddle.size;
		   int m = n / r;
		   int i;
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

		   break;
	      }

	 case FFTW_GENERIC:
	      {
		   int r = p->nodeu.generic.size;
		   int m = n / r;
		   int i;
		   fftw_generic_codelet *codelet;
		   fftw_complex *W;

		   for (i = 0; i < r; ++i) {
			fftw_executor_simple(m, in + i * istride,
					     out + i * (m * ostride),
					     p->nodeu.generic.recurse,
					     istride * r, ostride);
		   }

		   codelet = p->nodeu.generic.codelet;
		   W = p->nodeu.generic.tw->twarray;
		   codelet(out, W, m, r, n, ostride);

		   break;
	      }

	 case FFTW_RADER:
	      {
		   int r = p->nodeu.rader.size;
		   int m = n / r;
		   int i;
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

		   break;
	      }

	 default:
	      fftw_die("BUG in executor: invalid plan\n");
	      break;
     }
}

static void executor_simple_inplace(int n, fftw_complex *in,
				    fftw_complex *out,
				    fftw_plan_node *p,
				    int istride)
{
     switch (p->type) {
	 case FFTW_NOTW:
	      HACK_ALIGN_STACK_ODD();
	      (p->nodeu.notw.codelet)(in, in, istride, istride);
	      break;

	 default:
	      {
		   fftw_complex *tmp;

		   if (out)
			tmp = out;
		   else
			tmp = (fftw_complex *)
			    fftw_malloc(n * sizeof(fftw_complex));

		   fftw_executor_simple(n, in, tmp, p, istride, 1);
		   fftw_strided_copy(n, tmp, istride, in);

		   if (!out)
			fftw_free(tmp);
	      }
     }
}

static void executor_many(int n, const fftw_complex *in,
			  fftw_complex *out,
			  fftw_plan_node *p,
			  int istride,
			  int ostride,
			  int howmany, int idist, int odist)
{
     switch (p->type) {
	 case FFTW_NOTW:
	      {
		   fftw_notw_codelet *codelet = p->nodeu.notw.codelet;
		   int s;

		   HACK_ALIGN_STACK_ODD();
		   for (s = 0; s < howmany; ++s)
			codelet(in + s * idist,
				out + s * odist,
				istride, ostride);
		   break;
	      }

	 default:
	      {
		   int s;
		   for (s = 0; s < howmany; ++s) {
			fftw_executor_simple(n, in + s * idist,
					     out + s * odist,
					     p, istride, ostride);
		   }
	      }
     }
}

static void executor_many_inplace(int n, fftw_complex *in,
				  fftw_complex *out,
				  fftw_plan_node *p,
				  int istride,
				  int howmany, int idist)
{
     switch (p->type) {
	 case FFTW_NOTW:
	      {
		   fftw_notw_codelet *codelet = p->nodeu.notw.codelet;
		   int s;

		   HACK_ALIGN_STACK_ODD();
		   for (s = 0; s < howmany; ++s)
			codelet(in + s * idist,
				in + s * idist,
				istride, istride);
		   break;
	      }

	 default:
	      {
		   int s;
		   fftw_complex *tmp;
		   if (out)
			tmp = out;
		   else
			tmp = (fftw_complex *)
			    fftw_malloc(n * sizeof(fftw_complex));

		   for (s = 0; s < howmany; ++s) {
			fftw_executor_simple(n,
					     in + s * idist,
					     tmp,
					     p, istride, 1);
			fftw_strided_copy(n, tmp, istride, in + s * idist);
		   }

		   if (!out)
			fftw_free(tmp);
	      }
     }
}

/* user interface */
void fftw(fftw_plan plan, int howmany, fftw_complex *in, int istride,
	  int idist, fftw_complex *out, int ostride, int odist)
{
     int n = plan->n;

     if (plan->flags & FFTW_IN_PLACE) {
	  if (howmany == 1) {
	       executor_simple_inplace(n, in, out, plan->root, istride);
	  } else {
	       executor_many_inplace(n, in, out, plan->root, istride, howmany,
				     idist);
	  }
     } else {
	  if (howmany == 1) {
	       fftw_executor_simple(n, in, out, plan->root, istride, ostride);
	  } else {
	       executor_many(n, in, out, plan->root, istride, ostride,
			     howmany, idist, odist);
	  }
     }
}

void fftw_one(fftw_plan plan, fftw_complex *in, fftw_complex *out)
{
     int n = plan->n;

     if (plan->flags & FFTW_IN_PLACE)
	  executor_simple_inplace(n, in, out, plan->root, 1);
     else
	  fftw_executor_simple(n, in, out, plan->root, 1, 1);
}
