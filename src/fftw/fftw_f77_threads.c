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

#include <fftw_threads.h>
#include <f77_func.h>

#ifdef F77_FUNC_ /* only compile wrappers if fortran mangling is known */

#ifdef __cplusplus
extern "C" {
#endif                          /* __cplusplus */

/************************************************************************/

void F77_FUNC_(fftw_f77_threads_init,FFTW_F77_THREADS_INIT) (int *ierr)
{
     *ierr = fftw_threads_init();
}

void F77_FUNC_(fftw_f77_threads,FFTW_F77_THREADS)
(int *nthreads, fftw_plan *p,
 int *howmany, fftw_complex *in, int *istride, int *idist,
 fftw_complex *out, int *ostride, int *odist)
{
     fftw_threads(*nthreads,*p,
		  *howmany,in,*istride,*idist,out,*ostride,*odist);
}

void F77_FUNC_(fftw_f77_threads_one,FFTW_F77_THREADS_ONE)
(int *nthreads, fftw_plan *p, fftw_complex *in, fftw_complex *out)
{
     fftw_threads_one(*nthreads,*p,in,out);
}

void F77_FUNC_(fftwnd_f77_threads,FFTWND_F77_THREADS)
(int *nthreads, fftwnd_plan *p,
 int *howmany, fftw_complex *in, int *istride, int *idist,
 fftw_complex *out, int *ostride, int *odist)
{
     fftwnd_threads(*nthreads,*p,
		    *howmany,in,*istride,*idist,out,*ostride,*odist);
}

void F77_FUNC_(fftwnd_f77_threads_one,FFTWND_F77_THREADS_ONE)
(int *nthreads, fftwnd_plan *p, fftw_complex *in, fftw_complex *out)
{
     fftwnd_threads_one(*nthreads,*p,in,out);
}

/****************************************************************************/

#ifdef __cplusplus
}                               /* extern "C" */
#endif                          /* __cplusplus */

#endif /* defined(F77_FUNC_) */
