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

#ifndef RFFTW_THREADS_H
#define RFFTW_THREADS_H

#include <rfftw.h>
#include <fftw_threads.h>

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/******************** User Interface *********************/

extern void rfftw_threads(int nthreads,
                   fftw_plan plan, int howmany, fftw_real *in, int istride,
                   int idist, fftw_real *out, int ostride, int odist);
extern void rfftw_threads_one(int nthread, fftw_plan plan,
			      fftw_real *in, fftw_real *out);

extern void rfftwnd_threads_real_to_complex(int nthreads, fftwnd_plan p,
					    int howmany,
					    fftw_real *in,
					    int istride, int idist,
					    fftw_complex *out,
					    int ostride, int odist);
extern void rfftwnd_threads_complex_to_real(int nthreads, fftwnd_plan p,
					    int howmany,
					    fftw_complex *in,
					    int istride, int idist,
					    fftw_real *out,
					    int ostride, int odist);
extern void rfftwnd_threads_one_real_to_complex(int nthreads, fftwnd_plan p,
						fftw_real *in,
						fftw_complex *out);
extern void rfftwnd_threads_one_complex_to_real(int nthreads, fftwnd_plan p,
						fftw_complex *in,
						fftw_real *out);

#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* RFFTW_THREADS_H */
