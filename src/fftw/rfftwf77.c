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

#include <rfftw.h>
#include <f77_func.h>

#ifdef F77_FUNC_ /* only compile wrappers if fortran mangling is known */

/* rfftwf77.c: FORTRAN-callable "wrappers" for some of the RFFTW routines.

   See also fftwf77.c. */

#ifdef __cplusplus
extern "C" {
#endif                          /* __cplusplus */

/************************************************************************/

void F77_FUNC_(rfftw_f77_create_plan,RFFTW_F77_CREATE_PLAN)
(fftw_plan *p, int *n, int *idir, int *flags)
{
     fftw_direction dir = *idir < 0 ? FFTW_FORWARD : FFTW_BACKWARD;

     *p = rfftw_create_plan(*n,dir,*flags);
}

void F77_FUNC_(rfftw_f77_destroy_plan,RFFTW_F77_DESTROY_PLAN)
(fftw_plan *p)
{
     rfftw_destroy_plan(*p);
}

void F77_FUNC_(rfftw_f77,RFFTW_F77)
(fftw_plan *p, int *howmany, fftw_real *in, int *istride, int *idist,
 fftw_real *out, int *ostride, int *odist)
{
     rfftw(*p,*howmany,in,*istride,*idist,out,*ostride,*odist);
}

void F77_FUNC_(rfftw_f77_one,RFFTW_F77_ONE)
(fftw_plan *p, fftw_real *in, fftw_real *out)
{
     rfftw_one(*p,in,out);
}

extern void fftw_reverse_int_array(int *a, int n);

void F77_FUNC_(rfftwnd_f77_create_plan,RFFTWND_F77_CREATE_PLAN)
(fftwnd_plan *p, int *rank, int *n, int *idir, int *flags)
{
     fftw_direction dir = *idir < 0 ? FFTW_FORWARD : FFTW_BACKWARD;

     fftw_reverse_int_array(n,*rank);  /* column-major -> row-major */
     *p = rfftwnd_create_plan(*rank,n,dir,*flags);
     fftw_reverse_int_array(n,*rank);  /* reverse back */
}

void F77_FUNC_(rfftw2d_f77_create_plan,RFFTW2D_F77_CREATE_PLAN)
(fftwnd_plan *p, int *nx, int *ny, int *idir, int *flags)
{
     fftw_direction dir = *idir < 0 ? FFTW_FORWARD : FFTW_BACKWARD;

     *p = rfftw2d_create_plan(*ny,*nx,dir,*flags);
}

void F77_FUNC_(rfftw3d_f77_create_plan,RFFTW3D_F77_CREATE_PLAN)
(fftwnd_plan *p, int *nx, int *ny, int *nz, int *idir, int *flags)
{
     fftw_direction dir = *idir < 0 ? FFTW_FORWARD : FFTW_BACKWARD;

     *p = rfftw3d_create_plan(*nz,*ny,*nx,dir,*flags);
}

void F77_FUNC_(rfftwnd_f77_destroy_plan,RFFTWND_F77_DESTROY_PLAN)
(fftwnd_plan *p)
{
     rfftwnd_destroy_plan(*p);
}

void F77_FUNC_(rfftwnd_f77_real_to_complex,RFFTWND_F77_REAL_TO_COMPLEX)
(fftwnd_plan *p, int *howmany, fftw_real *in, int *istride, int *idist,
 fftw_complex *out, int *ostride, int *odist)
{
     rfftwnd_real_to_complex(*p,*howmany,in,*istride,*idist,
			     out,*ostride,*odist);
}

void F77_FUNC_(rfftwnd_f77_one_real_to_complex,RFFTWND_F77_ONE_REAL_TO_COMPLEX)
(fftwnd_plan *p, fftw_real *in, fftw_complex *out)
{
     rfftwnd_one_real_to_complex(*p,in,out);
}

void F77_FUNC_(rfftwnd_f77_complex_to_real,RFFTWND_F77_COMPLEX_TO_REAL)
(fftwnd_plan *p, int *howmany, fftw_complex *in, int *istride, int *idist,
 fftw_real *out, int *ostride, int *odist)
{
     rfftwnd_complex_to_real(*p,*howmany,in,*istride,*idist,
			     out,*ostride,*odist);
}

void F77_FUNC_(rfftwnd_f77_one_complex_to_real,RFFTWND_F77_ONE_COMPLEX_TO_REAL)
(fftwnd_plan *p, fftw_complex *in, fftw_real *out)
{
     rfftwnd_one_complex_to_real(*p,in,out);
}

/****************************************************************************/

#ifdef __cplusplus
}                               /* extern "C" */
#endif                          /* __cplusplus */

#endif /* defined(F77_FUNC_) */
