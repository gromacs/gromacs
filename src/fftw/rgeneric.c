
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
 * rgeneric.c -- "generic" rfftw codelets.  They work for all n (and
 * they are slow)
 */
#include <fftw-int.h>
#include <rfftw.h>

/* this code assumes that r and m are both odd */
void fftw_hc2hc_forward_generic(fftw_real *A, const fftw_complex *W,
				int m, int r, int n, int dist)
{
     int i, j, k;
     fftw_complex *tmp = (fftw_complex *)
     fftw_malloc(r * sizeof(fftw_complex));
     fftw_real rsum, isum;
     fftw_real *X, *YO, *YI;
     int wp, wincr;
     int iostride = m * dist;
     X = A;
     YO = A + r * iostride;
     YI = A + iostride;

     /* compute the transform of the r 0th elements (which are real) */
     for (i = 0; i + i < r; ++i) {
	  rsum = 0.0;
	  isum = 0.0;
	  wincr = m * i;
	  for (j = 0, wp = 0; j < r; ++j) {
	       fftw_real tw_r = c_re(W[wp]);
	       fftw_real tw_i = c_im(W[wp]);
	       fftw_real re = X[j * iostride];
	       rsum += re * tw_r;
	       isum += re * tw_i;
	       wp += wincr;
	       if (wp >= n)
		    wp -= n;
	  }
	  c_re(tmp[i]) = rsum;
	  c_im(tmp[i]) = isum;
     }

     /* store the transform back onto the A array */
     X[0] = c_re(tmp[0]);
     for (i = 1; i + i < r; ++i) {
	  X[i * iostride] = c_re(tmp[i]);
	  YO[-i * iostride] = c_im(tmp[i]);
     }

     X += dist;
     YI -= dist;
     YO -= dist;

     /* compute the transform of the middle elements (which are complex) */
     for (k = 1; k + k < m; ++k, X += dist, YI -= dist, YO -= dist) {
	  for (i = 0; i < r; ++i) {
	       rsum = 0.0;
	       isum = 0.0;
	       wincr = k + m * i;
	       for (j = 0, wp = 0; j < r; ++j) {
		    fftw_real tw_r = c_re(W[wp]);
		    fftw_real tw_i = c_im(W[wp]);
		    fftw_real re = X[j * iostride];
		    fftw_real im = YI[j * iostride];
		    rsum += re * tw_r - im * tw_i;
		    isum += re * tw_i + im * tw_r;
		    wp += wincr;
		    if (wp >= n)
			 wp -= n;
	       }
	       c_re(tmp[i]) = rsum;
	       c_im(tmp[i]) = isum;
	  }

	  /* store the transform back onto the A array */
	  for (i = 0; i + i < r; ++i) {
	       X[i * iostride] = c_re(tmp[i]);
	       YO[-i * iostride] = c_im(tmp[i]);
	  }
	  for (; i < r; ++i) {
	       X[i * iostride] = -c_im(tmp[i]);
	       YO[-i * iostride] = c_re(tmp[i]);
	  }
     }

     /* no final element, since m is odd */
     fftw_free(tmp);
}

void fftw_hc2hc_backward_generic(fftw_real *A, const fftw_complex *W,
				 int m, int r, int n, int dist)
{
     int i, j, k;
     int wp, wincr;
     fftw_complex *tmp = (fftw_complex *)
     fftw_malloc(r * sizeof(fftw_complex));
     fftw_real rsum, isum;
     fftw_real *X, *YO, *YI;
     int iostride = m * dist;
     X = A;
     YO = A + iostride;
     YI = A + r * iostride;

     /* 
      * compute the transform of the r 0th elements (which are halfcomplex)
      * yielding real numbers
      */
     /* copy the input into the temporary array */
     c_re(tmp[0]) = X[0];
     for (i = 1; i + i < r; ++i) {
	  c_re(tmp[i]) = X[i * iostride];
	  c_im(tmp[i]) = YI[-i * iostride];
     }

     for (i = 0; i < r; ++i) {
	  rsum = 0.0;
	  wincr = m * i;
	  for (j = 1, wp = wincr; j + j < r; ++j) {
	       fftw_real tw_r = c_re(W[wp]);
	       fftw_real tw_i = c_im(W[wp]);
	       fftw_real re = c_re(tmp[j]);
	       fftw_real im = c_im(tmp[j]);
	       rsum += re * tw_r + im * tw_i;
	       wp += wincr;
	       if (wp >= n)
		    wp -= n;
	  }
	  X[i * iostride] = 2.0 * rsum + c_re(tmp[0]);
     }

     X += dist;
     YI -= dist;
     YO -= dist;

     /* compute the transform of the middle elements (which are complex) */
     for (k = 1; k + k < m; ++k, X += dist, YI -= dist, YO -= dist) {
	  /* copy the input into the temporary array */
	  for (i = 0; i + i < r; ++i) {
	       c_re(tmp[i]) = X[i * iostride];
	       c_im(tmp[i]) = YI[-i * iostride];
	  }
	  for (; i < r; ++i) {
	       c_im(tmp[i]) = -X[i * iostride];
	       c_re(tmp[i]) = YI[-i * iostride];
	  }

	  for (i = 0; i < r; ++i) {
	       rsum = 0.0;
	       isum = 0.0;
	       wincr = m * i;
	       for (j = 0, wp = k * i; j < r; ++j) {
		    fftw_real tw_r = c_re(W[wp]);
		    fftw_real tw_i = c_im(W[wp]);
		    fftw_real re = c_re(tmp[j]);
		    fftw_real im = c_im(tmp[j]);
		    rsum += re * tw_r + im * tw_i;
		    isum += im * tw_r - re * tw_i;
		    wp += wincr;
		    if (wp >= n)
			 wp -= n;
	       }
	       X[i * iostride] = rsum;
	       YO[i * iostride] = isum;
	  }
     }

     /* no final element, since m is odd */
     fftw_free(tmp);
}
