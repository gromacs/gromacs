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

/*
 * executor.c -- execute the fft
 */

/* $Id$ */
#include <fftw.h>
#include <stdio.h>
#include <stdlib.h>

char *fftw_version = "FFTW V1.1 ($Id$)";

/*
 * This function is called in other files, so we cannot declare
 * it as static. 
 */

void fftw_strided_copy(int n, FFTW_COMPLEX *in, int ostride,
		       FFTW_COMPLEX *out)
{
     int i;
     FFTW_REAL r0, r1, i0, i1;
     FFTW_REAL r2, r3, i2, i3;

     i = 0;
     if (n & 3)
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
 * Do *not* declare simple executor as static--we need to call it
 * from executor_cilk.cilk...also, preface its name with "fftw_"
 * to avoid any possible name collisions. 
 */
void fftw_executor_simple(int n, const FFTW_COMPLEX *in,
			  FFTW_COMPLEX *out,
			  fftw_plan_node *p,
			  int istride,
			  int ostride)
{
     switch (p->type) {
	 case FFTW_NOTW:
	      (p->nodeu.notw.codelet) (in, out, istride, ostride);
	      break;

	 case FFTW_TWIDDLE:
	      {
		   int r = p->nodeu.twiddle.size;
		   int m = n / r;
		   int i;
		   twiddle_codelet *codelet;
		   FFTW_COMPLEX *W;

		   for (i = 0; i < r; ++i) {
			fftw_executor_simple(m, in + i * istride,
					     out + i * (m * ostride),
					     p->nodeu.twiddle.recurse,
					     istride * r, ostride);
		   }

		   codelet = p->nodeu.twiddle.codelet;
		   W = p->nodeu.twiddle.tw->twarray;
		   codelet(out, W, m * ostride, m, ostride);

		   break;
	      }

	 case FFTW_GENERIC:
	      {
		   int r = p->nodeu.generic.size;
		   int m = n / r;
		   int i;
		   generic_codelet *codelet;
		   FFTW_COMPLEX *W;

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

	 default:
	      fftw_die("BUG in executor: illegal plan\n");
	      break;
     }
}

static void executor_simple_inplace(int n, FFTW_COMPLEX *in,
				    FFTW_COMPLEX *out,
				    fftw_plan_node *p,
				    int istride)
{
     switch (p->type) {
	 case FFTW_NOTW:
	      (p->nodeu.notw.codelet) (in, in, istride, istride);
	      break;

	 default:
	      {
		   FFTW_COMPLEX *tmp;

		   if (out)
			tmp = out;
		   else
			tmp = (FFTW_COMPLEX *)
			    fftw_malloc(n * sizeof(FFTW_COMPLEX));

		   fftw_executor_simple(n, in, tmp, p, istride, 1);
		   fftw_strided_copy(n, tmp, istride, in);

		   if (!out)
			fftw_free(tmp);
	      }
     }
}

static void executor_many(int n, const FFTW_COMPLEX *in,
			  FFTW_COMPLEX *out,
			  fftw_plan_node *p,
			  int istride,
			  int ostride,
			  int howmany, int idist, int odist)
{
     switch (p->type) {
	 case FFTW_NOTW:
	      {
		   int s;
		   notw_codelet *codelet = p->nodeu.notw.codelet;
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

static void executor_many_inplace(int n, FFTW_COMPLEX *in,
				  FFTW_COMPLEX *out,
				  fftw_plan_node *p,
				  int istride,
				  int howmany, int idist)
{
     switch (p->type) {
	 case FFTW_NOTW:
	      {
		   int s;
		   notw_codelet *codelet = p->nodeu.notw.codelet;
		   for (s = 0; s < howmany; ++s)
			codelet(in + s * idist,
				in + s * idist,
				istride, istride);
		   break;
	      }

	 default:
	      {
		   int s;
		   FFTW_COMPLEX *tmp;
		   if (out)
			tmp = out;
		   else
			tmp = (FFTW_COMPLEX *)
			    fftw_malloc(n * sizeof(FFTW_COMPLEX));

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
void fftw(fftw_plan plan, int howmany, FFTW_COMPLEX *in, int istride,
	  int idist, FFTW_COMPLEX *out, int ostride, int odist)
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
