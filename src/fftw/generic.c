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
 *
 * generic.c -- "generic" codelets.  They work for all n (and they are
 * slow)
 */
#include <fftw-int.h>
#include <stdlib.h>

void fftw_twiddle_generic(fftw_complex *A, const fftw_complex *W,
			  int m, int r, int n, int stride)
{
     int i, j, k;
     const fftw_complex *jp;
     fftw_complex *kp;
     fftw_complex *tmp = (fftw_complex *)
     fftw_malloc(r * sizeof(fftw_complex));

     for (i = 0; i < m; ++i) {
	  for (k = 0, kp = tmp; k < r; ++k, kp++) {
	       fftw_real r0, i0, rt, it, rw, iw;
	       int l1 = i + m * k;
	       int l0;

	       r0 = i0 = 0.0;
	       for (j = 0, jp = A + i * stride, l0 = 0; j < r; ++j,
		    jp += m * stride) {
		    rw = c_re(W[l0]);
		    iw = c_im(W[l0]);
		    rt = c_re(*jp);
		    it = c_im(*jp);
		    r0 += rt * rw - it * iw;
		    i0 += rt * iw + it * rw;
		    l0 += l1;
		    if (l0 >= n)
			 l0 -= n;
	       }
	       c_re(*kp) = r0;
	       c_im(*kp) = i0;
	  }
	  for (k = 0, kp = A + i * stride; k < r; ++k, kp += m * stride)
	       *kp = tmp[k];
     }

     fftw_free(tmp);
}

void fftwi_twiddle_generic(fftw_complex *A, const fftw_complex *W,
			   int m, int r, int n, int stride)
{
     int i, j, k;
     const fftw_complex *jp;
     fftw_complex *kp;
     fftw_complex *tmp = (fftw_complex *)
     fftw_malloc(r * sizeof(fftw_complex));

     for (i = 0; i < m; ++i) {
	  for (k = 0, kp = tmp; k < r; ++k, kp++) {
	       fftw_real r0, i0, rt, it, rw, iw;
	       int l1 = i + m * k;
	       int l0;

	       r0 = i0 = 0.0;
	       for (j = 0, jp = A + i * stride, l0 = 0; j < r; ++j,
		    jp += m * stride) {
		    rw = c_re(W[l0]);
		    iw = c_im(W[l0]);
		    rt = c_re(*jp);
		    it = c_im(*jp);
		    r0 += rt * rw + it * iw;
		    i0 += it * rw - rt * iw;
		    l0 += l1;
		    if (l0 >= n)
			 l0 -= n;
	       }
	       c_re(*kp) = r0;
	       c_im(*kp) = i0;
	  }
	  for (k = 0, kp = A + i * stride; k < r; ++k, kp += m * stride)
	       *kp = tmp[k];
     }

     fftw_free(tmp);
}
