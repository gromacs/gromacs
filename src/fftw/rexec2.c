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
/*
 * rexec2.c -- alternate rfftw executor, specifically designed for the
 *             multidimensional transforms.  Given an extra work array,
 *             expects complex data in FFTW_COMPLEX format, and does
 *             not destroy the input in hc2real transforms.
 */

#include <fftw-int.h>
#include <rfftw.h>

/* copies halfcomplex array in (contiguous) to fftw_complex array out. */
void rfftw_hc2c(int n, fftw_real *in, fftw_complex *out, int ostride)
{
     int n2 = (n + 1) / 2;
     int i = 1;

     c_re(out[0]) = in[0];
     c_im(out[0]) = 0.0;
     for (; i < ((n2 - 1) & 3) + 1; ++i) {
	  c_re(out[i * ostride]) = in[i];
	  c_im(out[i * ostride]) = in[n - i];
     }
     for (; i < n2; i += 4) {
	  fftw_real r0, r1, r2, r3;
	  fftw_real i0, i1, i2, i3;
	  r0 = in[i];
	  r1 = in[i + 1];
	  r2 = in[i + 2];
	  r3 = in[i + 3];
	  i3 = in[n - (i + 3)];
	  i2 = in[n - (i + 2)];
	  i1 = in[n - (i + 1)];
	  i0 = in[n - i];
	  c_re(out[i * ostride]) = r0;
	  c_im(out[i * ostride]) = i0;
	  c_re(out[(i + 1) * ostride]) = r1;
	  c_im(out[(i + 1) * ostride]) = i1;
	  c_re(out[(i + 2) * ostride]) = r2;
	  c_im(out[(i + 2) * ostride]) = i2;
	  c_re(out[(i + 3) * ostride]) = r3;
	  c_im(out[(i + 3) * ostride]) = i3;
     }
     if ((n & 1) == 0) {	/* store the Nyquist frequency */
	  c_re(out[n2 * ostride]) = in[n2];
	  c_im(out[n2 * ostride]) = 0.0;
     }
}

/* reverse of rfftw_hc2c */
void rfftw_c2hc(int n, fftw_complex *in, int istride, fftw_real *out)
{
     int n2 = (n + 1) / 2;
     int i = 1;

     out[0] = c_re(in[0]);
     for (; i < ((n2 - 1) & 3) + 1; ++i) {
	  out[i] = c_re(in[i * istride]);
	  out[n - i] = c_im(in[i * istride]);
     }
     for (; i < n2; i += 4) {
	  fftw_real r0, r1, r2, r3;
	  fftw_real i0, i1, i2, i3;
	  r0 = c_re(in[i * istride]);
	  i0 = c_im(in[i * istride]);
	  r1 = c_re(in[(i + 1) * istride]);
	  i1 = c_im(in[(i + 1) * istride]);
	  r2 = c_re(in[(i + 2) * istride]);
	  i2 = c_im(in[(i + 2) * istride]);
	  r3 = c_re(in[(i + 3) * istride]);
	  i3 = c_im(in[(i + 3) * istride]);
	  out[i] = r0;
	  out[i + 1] = r1;
	  out[i + 2] = r2;
	  out[i + 3] = r3;
	  out[n - (i + 3)] = i3;
	  out[n - (i + 2)] = i2;
	  out[n - (i + 1)] = i1;
	  out[n - i] = i0;
     }
     if ((n & 1) == 0)		/* store the Nyquist frequency */
	  out[n2] = c_re(in[n2 * istride]);
}

/* 
 * in: array of n real numbers (* howmany).
 * out: array of n/2 + 1 complex numbers (* howmany).
 * work: array of n real numbers (stride 1) 
 * 
 * We must have out != in if dist < stride. 
 */
void rfftw_real2c_aux(fftw_plan plan, int howmany,
		      fftw_real *in, int istride, int idist,
		      fftw_complex *out, int ostride, int odist,
		      fftw_real *work)
{
     fftw_plan_node *p = plan->root;
     int j;

     switch (p->type) {
	 case FFTW_REAL2HC:
	      {
		   fftw_real2hc_codelet *codelet = p->nodeu.real2hc.codelet;
		   int n = plan->n;
		   int n2 = (n & 1) ? 0 : (n + 1) / 2;

		   HACK_ALIGN_STACK_ODD();
		   for (j = 0; j < howmany; ++j, out += odist) {
			codelet(in + j * idist,
				&c_re(*out),
				&c_im(*out),
				istride, ostride * 2, ostride * 2);
			c_im(out[0]) = 0.0;
			c_im(out[n2 * ostride]) = 0.0;
		   }
		   break;
	      }

	 default:
	      {
		   int n = plan->n;

		   for (j = 0; j < howmany; ++j, in += idist, out += odist) {
			rfftw_executor_simple(n, in, work, p, istride, 1);
			rfftw_hc2c(n, work, out, ostride);
		   }
		   break;
	      }
     }
}

/*
 * in: array of n/2 + 1 complex numbers (* howmany).
 * out: array of n real numbers (* howmany).
 * work: array of n real numbers (stride 1)
 * 
 * We must have out != in if dist < stride.  
 */
void rfftw_c2real_aux(fftw_plan plan, int howmany,
		      fftw_complex *in, int istride, int idist,
		      fftw_real *out, int ostride, int odist,
		      fftw_real *work)
{
     fftw_plan_node *p = plan->root;

     switch (p->type) {
	 case FFTW_HC2REAL:
	      {
		   fftw_hc2real_codelet *codelet = p->nodeu.hc2real.codelet;
		   int j;

		   HACK_ALIGN_STACK_ODD();
		   for (j = 0; j < howmany; ++j)
			codelet(&c_re(*(in + j * idist)),
				&c_im(*(in + j * idist)),
				out + j * odist,
				istride * 2, istride * 2, ostride);
		   break;
	      }

	 default:
	      {
		   int j, n = plan->n;

		   for (j = 0; j < howmany; ++j, in += idist, out += odist) {
			rfftw_c2hc(n, in, istride, work);
			rfftw_executor_simple(n, work, out, p, 1, ostride);
		   }
		   break;
	      }
     }
}

/*
 * The following two functions are similar to the ones above, BUT:
 * 
 * work must contain n * howmany elements (stride 1)
 * 
 * Can handle out == in for any stride/dist. 
 */
void rfftw_real2c_overlap_aux(fftw_plan plan, int howmany,
			      fftw_real *in, int istride, int idist,
			      fftw_complex *out, int ostride, int odist,
			      fftw_real *work)
{
     int n = plan->n;
     int j;

     rfftw(plan, howmany, in, istride, idist, work, 1, n);

     /* copy from work to out: */
     for (j = 0; j < howmany; ++j, work += n, out += odist)
	  rfftw_hc2c(n, work, out, ostride);
}

void rfftw_c2real_overlap_aux(fftw_plan plan, int howmany,
			      fftw_complex *in, int istride, int idist,
			      fftw_real *out, int ostride, int odist,
			      fftw_real *work)
{
     int n = plan->n;
     int j;

     /* copy from in to work: */
     for (j = 0; j < howmany; ++j, in += idist)
	  rfftw_c2hc(n, in, istride, work + j * n);

     rfftw(plan, howmany, work, 1, n, out, ostride, odist);
}
