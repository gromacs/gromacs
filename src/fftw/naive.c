
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

/* $Id$ */
#include <fftw.h>
#include <math.h>

/*
 * Naive O(n^2) algorithm, used for testing purposes
 */
void fftw_naive(int n, FFTW_COMPLEX *in, FFTW_COMPLEX *out)
{
     int i, j;
     FFTW_COMPLEX sum;
     FFTW_COMPLEX w;
     FFTW_REAL pi = 3.1415926535897932384626434;

     for (j = 0; j < n; ++j) {
	  c_re(sum) = c_im(sum) = 0.0;
	  for (i = 0; i < n; ++i) {
	       c_re(w) = cos((2.0 * pi * (i * j % n)) / n);
	       c_im(w) = -sin((2.0 * pi * (i * j % n)) / n);
	       c_re(sum) += c_re(in[i]) * c_re(w) - c_im(in[i]) * c_im(w);
	       c_im(sum) += c_im(in[i]) * c_re(w) + c_re(in[i]) * c_im(w);
	  }
	  out[j] = sum;
     }
     return;
}

/*
 * Naive O(n^2) algorithm, for the inverse.
 */
void fftwi_naive(int n, FFTW_COMPLEX *in, FFTW_COMPLEX *out)
{
     int i, j;
     FFTW_COMPLEX sum;
     FFTW_COMPLEX w;
     FFTW_REAL pi = 3.1415926535897932384626434;

     for (j = 0; j < n; ++j) {
	  c_re(sum) = c_im(sum) = 0.0;
	  for (i = 0; i < n; ++i) {
	       c_re(w) = cos((2.0 * pi * (i * j % n)) / n);
	       c_im(w) = sin((2.0 * pi * (i * j % n)) / n);
	       c_re(sum) += c_re(in[i]) * c_re(w) - c_im(in[i]) * c_im(w);
	       c_im(sum) += c_im(in[i]) * c_re(w) + c_re(in[i]) * c_im(w);
	  }
	  out[j] = sum;
     }
     return;
}
