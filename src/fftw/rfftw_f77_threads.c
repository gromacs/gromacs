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

#include <rfftw_threads.h>
#include <f77_func.h>

#ifdef F77_FUNC_ /* only compile wrappers if fortran mangling is known */

#ifdef __cplusplus
extern "C" {
#endif                          /* __cplusplus */

/************************************************************************/

void F77_FUNC_(rfftw_f77_threads,RFFTW_F77_THREADS)
(int *nthreads, fftw_plan *p,
 int *howmany, fftw_real *in, int *istride, int *idist,
 fftw_real *out, int *ostride, int *odist)
{
     rfftw_threads(*nthreads,*p,
		   *howmany,in,*istride,*idist,out,*ostride,*odist);
}

void F77_FUNC_(rfftw_f77_threads_one,RFFTW_F77_THREADS_ONE)
(int *nthreads, fftw_plan *p, fftw_real *in, fftw_real *out)
{
     rfftw_threads_one(*nthreads,*p,in,out);
}

void F77_FUNC_(rfftwnd_f77_threads_real_to_complex,RFFTWND_F77_THREADS_REAL_TO_COMPLEX)
(int *nthreads, fftwnd_plan *p,
 int *howmany, fftw_real *in, int *istride, int *idist,
 fftw_complex *out, int *ostride, int *odist)
{
     rfftwnd_threads_real_to_complex(*nthreads,*p,*howmany,in,*istride,*idist,
				     out,*ostride,*odist);
}

void F77_FUNC_(rfftwnd_f77_threads_one_real_to_complex,RFFTWND_F77_THREADS_ONE_REAL_TO_COMPLEX)
(int *nthreads, fftwnd_plan *p, fftw_real *in, fftw_complex *out)
{
     rfftwnd_threads_one_real_to_complex(*nthreads,*p,in,out);
}

void F77_FUNC_(rfftwnd_f77_threads_complex_to_real,RFFTWND_F77_THREADS_COMPLEX_TO_REAL)
(int *nthreads, fftwnd_plan *p,
 int *howmany, fftw_complex *in, int *istride, int *idist,
 fftw_real *out, int *ostride, int *odist)
{
     rfftwnd_threads_complex_to_real(*nthreads,*p,*howmany,in,*istride,*idist,
				     out,*ostride,*odist);
}

void F77_FUNC_(rfftwnd_f77_threads_one_complex_to_real,RFFTWND_F77_THREADS_ONE_COMPLEX_TO_REAL)
(int *nthreads, fftwnd_plan *p, fftw_complex *in, fftw_real *out)
{
     rfftwnd_threads_one_complex_to_real(*nthreads,*p,in,out);
}

/****************************************************************************/

#ifdef __cplusplus
}                               /* extern "C" */
#endif                          /* __cplusplus */

#endif /* defined(F77_FUNC_) */
