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
 * generic.c -- "generic" solvers.  They work for all
 * n (and are slow)
 */
#include <fftw.h>
#include <math.h>
#include <stdlib.h>

void fftw_twiddle_generic(FFTW_COMPLEX *A, const FFTW_COMPLEX *W,
			  int m, int r, int n, int stride)
{
     int i, j, k;
     const FFTW_COMPLEX *jp;
     FFTW_COMPLEX *kp;
     FFTW_COMPLEX *tmp = (FFTW_COMPLEX *)
     fftw_malloc(r * sizeof(FFTW_COMPLEX));

     for (i = 0; i < m; ++i) {
	  for (k = 0, kp = tmp; k < r; ++k, kp++) {
	       FFTW_REAL r0, i0, rt, it, rw, iw;
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
		    if (l0 > n)
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

void fftwi_twiddle_generic(FFTW_COMPLEX *A, const FFTW_COMPLEX *W,
			   int m, int r, int n, int stride)
{
     int i, j, k;
     const FFTW_COMPLEX *jp;
     FFTW_COMPLEX *kp;
     FFTW_COMPLEX *tmp = (FFTW_COMPLEX *)
     fftw_malloc(r * sizeof(FFTW_COMPLEX));

     for (i = 0; i < m; ++i) {
	  for (k = 0, kp = tmp; k < r; ++k, kp++) {
	       FFTW_REAL r0, i0, rt, it, rw, iw;
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
		    if (l0 > n)
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
