/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Gromacs 4.0                         Copyright (c) 1991-2003
 * David van der Spoel, Erik Lindahl, University of Groningen.
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef GMX_FFT_FFTPACK

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>


#include "gmx_fft.h"
#include "gmx_fatal.h"
#include "fftpack.h"

/** Contents of the FFTPACK fft datatype.
 *
 *  FFTPACK only does 1d transforms, so we use a pointers to another fft for
 *  the transform in the next dimension.
 * Thus, a 3d-structure contains a pointer to a 2d one, which in turns contains
 * a pointer to a 1d. The 1d structure has next==NULL.
 */
struct gmx_fft
{
    int             ndim;     /**< Dimensions, including our subdimensions.  */
    int             n;        /**< Number of points in this dimension.       */
    int             ifac[15]; /**< 15 bytes needed for cfft and rfft         */
    struct gmx_fft *next;     /**< Pointer to next dimension, or NULL.       */
    real  *         work;     /**< 1st 4n reserved for cfft, 1st 2n for rfft */
};

#include <math.h>
#include <stdio.h>


int
gmx_fft_init_1d(gmx_fft_t *        pfft,
                int                nx,
                int                flags)
{
    gmx_fft_t    fft;

    if (pfft == NULL)
    {
        gmx_fatal(FARGS, "Invalid FFT opaque type pointer.");
        return EINVAL;
    }
    *pfft = NULL;

    if ( (fft = (struct gmx_fft *)malloc(sizeof(struct gmx_fft))) == NULL)
    {
        return ENOMEM;
    }

    fft->next = NULL;
    fft->n    = nx;

    /* Need 4*n storage for 1D complex FFT */
    if ( (fft->work = (real *)malloc(sizeof(real)*(4*nx))) == NULL)
    {
        free(fft);
        return ENOMEM;
    }

    if (fft->n > 1)
    {
        fftpack_cffti1(nx, fft->work, fft->ifac);
    }

    *pfft = fft;
    return 0;
};



int
gmx_fft_init_1d_real(gmx_fft_t *        pfft,
                     int                nx,
                     int                flags)
{
    gmx_fft_t    fft;

    if (pfft == NULL)
    {
        gmx_fatal(FARGS, "Invalid FFT opaque type pointer.");
        return EINVAL;
    }
    *pfft = NULL;

    if ( (fft = (struct gmx_fft *)malloc(sizeof(struct gmx_fft))) == NULL)
    {
        return ENOMEM;
    }

    fft->next = NULL;
    fft->n    = nx;

    /* Need 2*n storage for 1D real FFT */
    if ((fft->work = (real *)malloc(sizeof(real)*(2*nx))) == NULL)
    {
        free(fft);
        return ENOMEM;
    }

    if (fft->n > 1)
    {
        fftpack_rffti1(nx, fft->work, fft->ifac);
    }

    *pfft = fft;
    return 0;
}



int
gmx_fft_init_2d(gmx_fft_t *        pfft,
                int                nx,
                int                ny,
                int                flags)
{
    gmx_fft_t     fft;
    int           rc;

    if (pfft == NULL)
    {
        gmx_fatal(FARGS, "Invalid FFT opaque type pointer.");
        return EINVAL;
    }
    *pfft = NULL;

    /* Create the X transform */
    if ( (rc = gmx_fft_init_1d(&fft, nx, flags)) != 0)
    {
        return rc;
    }

    /* Create Y transform as a link from X */
    if ( (rc = gmx_fft_init_1d(&(fft->next), ny, flags)) != 0)
    {
        free(fft);
        return rc;
    }

    *pfft = fft;
    return 0;
};


int
gmx_fft_init_2d_real(gmx_fft_t *        pfft,
                     int                nx,
                     int                ny,
                     int                flags)
{
    gmx_fft_t     fft;
    int           nyc = (ny/2 + 1);
    int           rc;

    if (pfft == NULL)
    {
        gmx_fatal(FARGS, "Invalid FFT opaque type pointer.");
        return EINVAL;
    }
    *pfft = NULL;

    /* Create the X transform */
    if ( (fft = (struct gmx_fft *)malloc(sizeof(struct gmx_fft))) == NULL)
    {
        return ENOMEM;
    }

    fft->n    = nx;

    /* Need 4*nx storage for 1D complex FFT, and another
     * 2*nx*nyc elements for complex-to-real storage in our high-level routine.
     */
    if ( (fft->work = (real *)malloc(sizeof(real)*(4*nx+2*nx*nyc))) == NULL)
    {
        free(fft);
        return ENOMEM;
    }
    fftpack_cffti1(nx, fft->work, fft->ifac);

    /* Create real Y transform as a link from X */
    if ( (rc = gmx_fft_init_1d_real(&(fft->next), ny, flags)) != 0)
    {
        free(fft);
        return rc;
    }

    *pfft = fft;
    return 0;
}


int
gmx_fft_init_3d(gmx_fft_t *        pfft,
                int                nx,
                int                ny,
                int                nz,
                int                flags)
{
    gmx_fft_t     fft;
    int           rc;

    if (pfft == NULL)
    {
        gmx_fatal(FARGS, "Invalid FFT opaque type pointer.");
        return EINVAL;
    }
    *pfft = NULL;

    /* Create the X transform */

    if ( (fft = (struct gmx_fft *)malloc(sizeof(struct gmx_fft))) == NULL)
    {
        return ENOMEM;
    }

    fft->n    = nx;

    /* Need 4*nx storage for 1D complex FFT, and another
     * 2*nz elements for gmx_fft_transpose_2d_nelem() storage.
     */
    if ( (fft->work = (real *)malloc(sizeof(real)*(4*nx+2*nz))) == NULL)
    {
        free(fft);
        return ENOMEM;
    }

    fftpack_cffti1(nx, fft->work, fft->ifac);


    /* Create 2D Y/Z transforms as a link from X */
    if ( (rc = gmx_fft_init_2d(&(fft->next), ny, nz, flags)) != 0)
    {
        free(fft);
        return rc;
    }

    *pfft = fft;
    return 0;
};


int
gmx_fft_init_3d_real(gmx_fft_t *        pfft,
                     int                nx,
                     int                ny,
                     int                nz,
                     int                flags)
{
    gmx_fft_t     fft;
    int           nzc = (nz/2 + 1);
    int           rc;

    if (pfft == NULL)
    {
        gmx_fatal(FARGS, "Invalid FFT opaque type pointer.");
        return EINVAL;
    }
    *pfft = NULL;

    /* Create the X transform */
    if ( (fft = (struct gmx_fft *)malloc(sizeof(struct gmx_fft))) == NULL)
    {
        return ENOMEM;
    }

    fft->n    = nx;

    /* Need 4*nx storage for 1D complex FFT, another
     * 2*nx*ny*nzc elements to copy the entire 3D matrix when
     * doing out-of-place complex-to-real FFTs, and finally
     * 2*nzc elements for transpose work space.
     */
    if ( (fft->work = (real *)malloc(sizeof(real)*(4*nx+2*nx*ny*nzc+2*nzc))) == NULL)
    {
        free(fft);
        return ENOMEM;
    }
    fftpack_cffti1(nx, fft->work, fft->ifac);

    /* Create 2D real Y/Z transform as a link from X */
    if ( (rc = gmx_fft_init_2d_real(&(fft->next), ny, nz, flags)) != 0)
    {
        free(fft);
        return rc;
    }

    *pfft = fft;
    return 0;
}


int
gmx_fft_1d               (gmx_fft_t                  fft,
                          enum gmx_fft_direction     dir,
                          void *                     in_data,
                          void *                     out_data)
{
    int             i, n;
    real       *    p1;
    real       *    p2;

    n = fft->n;

    if (n == 1)
    {
        p1    = (real *)in_data;
        p2    = (real *)out_data;
        p2[0] = p1[0];
        p2[1] = p1[1];
    }

    /* FFTPACK only does in-place transforms, so emulate out-of-place
     * by copying data to the output array first.
     */
    if (in_data != out_data)
    {
        p1 = (real *)in_data;
        p2 = (real *)out_data;

        /* n complex = 2*n real elements */
        for (i = 0; i < 2*n; i++)
        {
            p2[i] = p1[i];
        }
    }

    /* Elements 0   .. 2*n-1 in work are used for ffac values,
     * Elements 2*n .. 4*n-1 are internal FFTPACK work space.
     */

    if (dir == GMX_FFT_FORWARD)
    {
        fftpack_cfftf1(n, (real *)out_data, fft->work+2*n, fft->work, fft->ifac, -1);
    }
    else if (dir == GMX_FFT_BACKWARD)
    {
        fftpack_cfftf1(n, (real *)out_data, fft->work+2*n, fft->work, fft->ifac, 1);
    }
    else
    {
        gmx_fatal(FARGS, "FFT plan mismatch - bad plan or direction.");
        return EINVAL;
    }

    return 0;
}



int
gmx_fft_1d_real          (gmx_fft_t                  fft,
                          enum gmx_fft_direction     dir,
                          void *                     in_data,
                          void *                     out_data)
{
    int           i, n;
    real       *  p1;
    real       *  p2;

    n = fft->n;

    if (n == 1)
    {
        p1    = (real *)in_data;
        p2    = (real *)out_data;
        p2[0] = p1[0];
        if (dir == GMX_FFT_REAL_TO_COMPLEX)
        {
            p2[1] = 0.0;
        }
    }

    if (dir == GMX_FFT_REAL_TO_COMPLEX)
    {
        /* FFTPACK only does in-place transforms, so emulate out-of-place
         * by copying data to the output array first. This works fine, since
         * the complex array must be larger than the real.
         */
        if (in_data != out_data)
        {
            p1 = (real *)in_data;
            p2 = (real *)out_data;

            for (i = 0; i < 2*(n/2+1); i++)
            {
                p2[i] = p1[i];
            }
        }

        /* Elements 0 ..   n-1 in work are used for ffac values,
         * Elements n .. 2*n-1 are internal FFTPACK work space.
         */
        fftpack_rfftf1(n, (real *)out_data, fft->work+n, fft->work, fft->ifac);

        /*
         * FFTPACK has a slightly more compact storage than we, time to
         * convert it: ove most of the array one step up to make room for
         * zero imaginary parts.
         */
        p2 = (real *)out_data;
        for (i = n-1; i > 0; i--)
        {
            p2[i+1] = p2[i];
        }
        /* imaginary zero freq. */
        p2[1] = 0;

        /* Is n even? */
        if ( (n & 0x1) == 0)
        {
            p2[n+1] = 0;
        }

    }
    else if (dir == GMX_FFT_COMPLEX_TO_REAL)
    {
        /* FFTPACK only does in-place transforms, and we cannot just copy
         * input to output first here since our real array is smaller than
         * the complex one. However, since the FFTPACK complex storage format
         * is more compact than ours (2 reals) it will fit, so compact it
         * and copy on-the-fly to the output array.
         */
        p1 = (real *) in_data;
        p2 = (real *)out_data;

        p2[0] = p1[0];
        for (i = 1; i < n; i++)
        {
            p2[i] = p1[i+1];
        }
        fftpack_rfftb1(n, (real *)out_data, fft->work+n, fft->work, fft->ifac);
    }
    else
    {
        gmx_fatal(FARGS, "FFT plan mismatch - bad plan or direction.");
        return EINVAL;
    }

    return 0;
}


int
gmx_fft_2d               (gmx_fft_t                  fft,
                          enum gmx_fft_direction     dir,
                          void *                     in_data,
                          void *                     out_data)
{
    int                i, nx, ny;
    t_complex     *    data;

    nx = fft->n;
    ny = fft->next->n;

    /* FFTPACK only does in-place transforms, so emulate out-of-place
     * by copying data to the output array first.
     * For 2D there is likely enough data to benefit from memcpy().
     */
    if (in_data != out_data)
    {
        memcpy(out_data, in_data, sizeof(t_complex)*nx*ny);
    }

    /* Much easier to do pointer arithmetic when base has the correct type */
    data = (t_complex *)out_data;

    /* y transforms */
    for (i = 0; i < nx; i++)
    {
        gmx_fft_1d(fft->next, dir, data+i*ny, data+i*ny);
    }

    /* Transpose in-place to get data in place for x transform now */
    gmx_fft_transpose_2d(data, data, nx, ny);

    /* x transforms */
    for (i = 0; i < ny; i++)
    {
        gmx_fft_1d(fft, dir, data+i*nx, data+i*nx);
    }

    /* Transpose in-place to get data back in original order */
    gmx_fft_transpose_2d(data, data, ny, nx);

    return 0;
}



int
gmx_fft_2d_real          (gmx_fft_t                  fft,
                          enum gmx_fft_direction     dir,
                          void *                     in_data,
                          void *                     out_data)
{
    int                i, j, nx, ny, nyc;
    t_complex     *    data;
    real       *       work;
    real       *       p1;
    real       *       p2;

    nx = fft->n;
    ny = fft->next->n;
    /* Number of complex elements in y direction */
    nyc = (ny/2+1);

    work = fft->work+4*nx;

    if (dir == GMX_FFT_REAL_TO_COMPLEX)
    {
        /* If we are doing an in-place transform the 2D array is already
         * properly padded by the user, and we are all set.
         *
         * For out-of-place there is no array padding, but FFTPACK only
         * does in-place FFTs internally, so we need to start by copying
         * data from the input to the padded (larger) output array.
         */
        if (in_data != out_data)
        {
            p1 = (real *)in_data;
            p2 = (real *)out_data;

            for (i = 0; i < nx; i++)
            {
                for (j = 0; j < ny; j++)
                {
                    p2[i*nyc*2+j] = p1[i*ny+j];
                }
            }
        }
        data = (t_complex *)out_data;

        /* y real-to-complex FFTs */
        for (i = 0; i < nx; i++)
        {
            gmx_fft_1d_real(fft->next, GMX_FFT_REAL_TO_COMPLEX, data+i*nyc, data+i*nyc);
        }

        /* Transform to get X data in place */
        gmx_fft_transpose_2d(data, data, nx, nyc);

        /* Complex-to-complex X FFTs */
        for (i = 0; i < nyc; i++)
        {
            gmx_fft_1d(fft, GMX_FFT_FORWARD, data+i*nx, data+i*nx);
        }

        /* Transpose back */
        gmx_fft_transpose_2d(data, data, nyc, nx);

    }
    else if (dir == GMX_FFT_COMPLEX_TO_REAL)
    {
        /* An in-place complex-to-real transform is straightforward,
         * since the output array must be large enough for the padding to fit.
         *
         * For out-of-place complex-to-real transforms we cannot just copy
         * data to the output array, since it is smaller than the input.
         * In this case there's nothing to do but employing temporary work data,
         * starting at work+4*nx and using nx*nyc*2 elements.
         */
        if (in_data != out_data)
        {
            memcpy(work, in_data, sizeof(t_complex)*nx*nyc);
            data = (t_complex *)work;
        }
        else
        {
            /* in-place */
            data = (t_complex *)out_data;
        }

        /* Transpose to get X arrays */
        gmx_fft_transpose_2d(data, data, nx, nyc);

        /* Do X iFFTs */
        for (i = 0; i < nyc; i++)
        {
            gmx_fft_1d(fft, GMX_FFT_BACKWARD, data+i*nx, data+i*nx);
        }

        /* Transpose to get Y arrays */
        gmx_fft_transpose_2d(data, data, nyc, nx);

        /* Do Y iFFTs */
        for (i = 0; i < nx; i++)
        {
            gmx_fft_1d_real(fft->next, GMX_FFT_COMPLEX_TO_REAL, data+i*nyc, data+i*nyc);
        }

        if (in_data != out_data)
        {
            /* Output (pointed to by data) is now in padded format.
             * Pack it into out_data if we were doing an out-of-place transform.
             */
            p1 = (real *)data;
            p2 = (real *)out_data;

            for (i = 0; i < nx; i++)
            {
                for (j = 0; j < ny; j++)
                {
                    p2[i*ny+j] = p1[i*nyc*2+j];
                }
            }
        }
    }
    else
    {
        gmx_fatal(FARGS, "FFT plan mismatch - bad plan or direction.");
        return EINVAL;
    }

    return 0;
}



int
gmx_fft_3d          (gmx_fft_t                  fft,
                     enum gmx_fft_direction     dir,
                     void *                     in_data,
                     void *                     out_data)
{
    int              i, nx, ny, nz, rc;
    t_complex     *  data;
    t_complex     *  work;
    nx = fft->n;
    ny = fft->next->n;
    nz = fft->next->next->n;

    /* First 4*nx positions are FFTPACK workspace, then ours starts */
    work = (t_complex *)(fft->work+4*nx);

    /* FFTPACK only does in-place transforms, so emulate out-of-place
     * by copying data to the output array first.
     * For 3D there is likely enough data to benefit from memcpy().
     */
    if (in_data != out_data)
    {
        memcpy(out_data, in_data, sizeof(t_complex)*nx*ny*nz);
    }

    /* Much easier to do pointer arithmetic when base has the correct type */
    data = (t_complex *)out_data;

    /* Perform z transforms */
    for (i = 0; i < nx*ny; i++)
    {
        gmx_fft_1d(fft->next->next, dir, data+i*nz, data+i*nz);
    }

    /* For each X slice, transpose the y & z dimensions inside the slice */
    for (i = 0; i < nx; i++)
    {
        gmx_fft_transpose_2d(data+i*ny*nz, data+i*ny*nz, ny, nz);
    }

    /* Array is now (nx,nz,ny) - perform y transforms */
    for (i = 0; i < nx*nz; i++)
    {
        gmx_fft_1d(fft->next, dir, data+i*ny, data+i*ny);
    }

    /* Transpose back to (nx,ny,nz) */
    for (i = 0; i < nx; i++)
    {
        gmx_fft_transpose_2d(data+i*ny*nz, data+i*ny*nz, nz, ny);
    }

    /* Transpose entire x & y slices to go from
     * (nx,ny,nz) to (ny,nx,nz).
     * Use work data elements 4*n .. 4*n+2*nz-1.
     */
    rc = gmx_fft_transpose_2d_nelem(data, data, nx, ny, nz, work);
    if (rc != 0)
    {
        gmx_fatal(FARGS, "Cannot transpose X & Y/Z in gmx_fft_3d().");
        return rc;
    }

    /* Then go from (ny,nx,nz) to (ny,nz,nx) */
    for (i = 0; i < ny; i++)
    {
        gmx_fft_transpose_2d(data+i*nx*nz, data+i*nx*nz, nx, nz);
    }

    /* Perform x transforms */
    for (i = 0; i < ny*nz; i++)
    {
        gmx_fft_1d(fft, dir, data+i*nx, data+i*nx);
    }

    /* Transpose back from (ny,nz,nx) to (ny,nx,nz) */
    for (i = 0; i < ny; i++)
    {
        gmx_fft_transpose_2d(data+i*nz*nx, data+i*nz*nx, nz, nx);
    }

    /* Transpose from (ny,nx,nz) to (nx,ny,nz)
     * Use work data elements 4*n .. 4*n+2*nz-1.
     */
    rc = gmx_fft_transpose_2d_nelem(data, data, ny, nx, nz, work);
    if (rc != 0)
    {
        gmx_fatal(FARGS, "Cannot transpose Y/Z & X in gmx_fft_3d().");
        return rc;
    }

    return 0;
}


int
gmx_fft_3d_real          (gmx_fft_t                  fft,
                          enum gmx_fft_direction     dir,
                          void *                     in_data,
                          void *                     out_data)
{
    int              i, j, k;
    int              nx, ny, nz, nzc, rc;
    t_complex     *  data;
    t_complex     *  work_transp;
    t_complex     *  work_c2r;
    real       *     p1;
    real       *     p2;

    nx  = fft->n;
    ny  = fft->next->n;
    nz  = fft->next->next->n;
    nzc = (nz/2+1);


    /* First 4*nx positions are FFTPACK workspace, then ours starts.
     * We have 2*nx*ny*nzc elements for temp complex-to-real storage when
     * doing out-of-place transforms, and another 2*nzc for transpose data.
     */
    work_c2r    = (t_complex *)(fft->work+4*nx);
    work_transp = (t_complex *)(fft->work+4*nx+2*nx*ny*nzc);

    /* Much easier to do pointer arithmetic when base has the correct type */
    data = (t_complex *)out_data;

    if (dir == GMX_FFT_REAL_TO_COMPLEX)
    {
        /* FFTPACK only does in-place transforms, so emulate out-of-place
         * by copying data to the output array first. This is guaranteed to
         * work for real-to-complex since complex data is larger than the real.
         * For 3D there is likely enough data to benefit from memcpy().
         */
        if (in_data != out_data)
        {
            p1 = (real *)in_data;
            p2 = (real *)out_data;

            for (i = 0; i < nx; i++)
            {
                for (j = 0; j < ny; j++)
                {
                    for (k = 0; k < nz; k++)
                    {
                        p2[(i*ny+j)*2*nzc+k] = p1[(i*ny+j)*nz+k];
                    }
                }
            }
        }
        data = (t_complex *)out_data;

        /* Transform the Y/Z slices real-to-complex */
        for (i = 0; i < nx; i++)
        {
            gmx_fft_2d_real(fft->next, dir, data+i*ny*nzc, data+i*ny*nzc);
        }

        /* Transpose x & y slices to go from
         * (nx,ny,nzc) to (ny,nx,nzc).
         */
        rc = gmx_fft_transpose_2d_nelem(data, data, nx, ny, nzc, work_transp);
        if (rc != 0)
        {
            gmx_fatal(FARGS, "Cannot transpose X & Y/Z gmx_fft_3d_real().");
            return rc;
        }

        /* Then transpose from (ny,nx,nzc) to (ny,nzc,nx) */
        for (i = 0; i < ny; i++)
        {
            gmx_fft_transpose_2d(data+i*nx*nzc, data+i*nx*nzc, nx, nzc);
        }

        /* Perform x transforms */
        for (i = 0; i < ny*nzc; i++)
        {
            gmx_fft_1d(fft, GMX_FFT_FORWARD, data+i*nx, data+i*nx);
        }

        /* Transpose from (ny,nzc,nx) back to (ny,nx,nzc) */
        for (i = 0; i < ny; i++)
        {
            gmx_fft_transpose_2d(data+i*nzc*nx, data+i*nzc*nx, nzc, nx);
        }

        /* Transpose back from (ny,nx,nzc) to (nx,ny,nz) */
        rc = gmx_fft_transpose_2d_nelem(data, data, ny, nx, nzc, work_transp);
        if (rc != 0)
        {
            gmx_fatal(FARGS, "Cannot transpose Y/Z & X in gmx_fft_3d_real().");
            return rc;
        }

    }
    else if (dir == GMX_FFT_COMPLEX_TO_REAL)
    {
        /* An in-place complex-to-real transform is straightforward,
         * since the output array must be large enough for the padding to fit.
         *
         * For out-of-place complex-to-real transforms we cannot just copy
         * data to the output array, since it is smaller than the input.
         * In this case there's nothing to do but employing temporary work data.
         */
        if (in_data != out_data)
        {
            memcpy(work_c2r, in_data, sizeof(t_complex)*nx*ny*nzc);
            data = (t_complex *)work_c2r;
        }
        else
        {
            /* in-place */
            data = (t_complex *)out_data;
        }

        /* Transpose x & y slices to go from
         * (nx,ny,nz) to (ny,nx,nz).
         */
        gmx_fft_transpose_2d_nelem(data, data, nx, ny, nzc, work_transp);

        /* Then go from (ny,nx,nzc) to (ny,nzc,nx) */
        for (i = 0; i < ny; i++)
        {
            gmx_fft_transpose_2d(data+i*nx*nzc, data+i*nx*nzc, nx, nzc);
        }


        /* Perform x transforms */
        for (i = 0; i < ny*nzc; i++)
        {
            gmx_fft_1d(fft, GMX_FFT_BACKWARD, data+i*nx, data+i*nx);
        }

        /* Transpose back from (ny,nzc,nx) to (ny,nx,nzc) */
        for (i = 0; i < ny; i++)
        {
            gmx_fft_transpose_2d(data+i*nzc*nx, data+i*nzc*nx, nzc, nx);
        }

        /* Transpose back from (ny,nx,nzc) to (nx,ny,nz) */
        gmx_fft_transpose_2d_nelem(data, data, ny, nx, nzc, work_transp);


        /* Do 2D complex-to-real */
        for (i = 0; i < nx; i++)
        {
            gmx_fft_2d_real(fft->next, dir, data+i*ny*nzc, data+i*ny*nzc);
        }

        if (in_data != out_data)
        {
            /* Output (pointed to by data) is now in padded format.
             * Pack it into out_data if we were doing an out-of-place transform.
             */
            p1 = (real *)data;
            p2 = (real *)out_data;

            for (i = 0; i < nx; i++)
            {
                for (j = 0; j < ny; j++)
                {
                    for (k = 0; k < nz; k++)
                    {
                        p2[(i*ny+j)*nz+k] = p1[(i*ny+j)*nzc*2+k];
                    }
                }
            }
        }

    }
    else
    {
        gmx_fatal(FARGS, "FFT plan mismatch - bad plan or direction.");
        return EINVAL;
    }

    return 0;
}




void
gmx_fft_destroy(gmx_fft_t      fft)
{
    if (fft != NULL)
    {
        free(fft->work);
        if (fft->next != NULL)
        {
            gmx_fft_destroy(fft->next);
        }
        free(fft);
    }
}
#endif /* GMX_FFT_FFTPACK */
