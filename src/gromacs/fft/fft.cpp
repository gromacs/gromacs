/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2003 Erik Lindahl, David van der Spoel, University of Groningen.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
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
#include "gmxpre.h"

#include "fft.h"

#include "config.h"

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gromacs/math/gmxcomplex.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/real.h"

/* This file contains common fft utility functions, but not
 * the actual transform implementations. Check the
 * files like fft_fftw3.c or fft_mkl.c for that.
 */

#ifndef GMX_FFT_FFTW3

struct gmx_many_fft {
    int       howmany;
    int       dist;
    gmx_fft_t fft;
};

typedef struct gmx_many_fft* gmx_many_fft_t;

int
gmx_fft_init_many_1d(gmx_fft_t *        pfft,
                     int                nx,
                     int                howmany,
                     gmx_fft_flag       flags)
{
    gmx_many_fft_t fft;
    if (pfft == NULL)
    {
        gmx_fatal(FARGS, "Invalid opaque FFT datatype pointer.");
        return EINVAL;
    }
    *pfft = NULL;

    if ( (fft = (gmx_many_fft_t)malloc(sizeof(struct gmx_many_fft))) == NULL)
    {
        return ENOMEM;
    }

    gmx_fft_init_1d(&fft->fft, nx, flags);
    fft->howmany = howmany;
    fft->dist    = 2*nx;

    *pfft = (gmx_fft_t)fft;
    return 0;
}

int
gmx_fft_init_many_1d_real(gmx_fft_t *        pfft,
                          int                nx,
                          int                howmany,
                          gmx_fft_flag       flags)
{
    gmx_many_fft_t fft;
    if (pfft == NULL)
    {
        gmx_fatal(FARGS, "Invalid opaque FFT datatype pointer.");
        return EINVAL;
    }
    *pfft = NULL;

    if ( (fft = (gmx_many_fft_t)malloc(sizeof(struct gmx_many_fft))) == NULL)
    {
        return ENOMEM;
    }

    gmx_fft_init_1d_real(&fft->fft, nx, flags);
    fft->howmany = howmany;
    fft->dist    = 2*(nx/2+1);

    *pfft = (gmx_fft_t)fft;
    return 0;
}

int
gmx_fft_many_1d     (gmx_fft_t                  fft,
                     enum gmx_fft_direction     dir,
                     void *                     in_data,
                     void *                     out_data)
{
    gmx_many_fft_t mfft = (gmx_many_fft_t)fft;
    int            i, ret;
    for (i = 0; i < mfft->howmany; i++)
    {
        ret = gmx_fft_1d(mfft->fft, dir, in_data, out_data);
        if (ret != 0)
        {
            return ret;
        }
        in_data  = (real*)in_data+mfft->dist;
        out_data = (real*)out_data+mfft->dist;
    }
    return 0;
}

int
gmx_fft_many_1d_real     (gmx_fft_t                  fft,
                          enum gmx_fft_direction     dir,
                          void *                     in_data,
                          void *                     out_data)
{
    gmx_many_fft_t mfft = (gmx_many_fft_t)fft;
    int            i, ret;
    for (i = 0; i < mfft->howmany; i++)
    {
        ret = gmx_fft_1d_real(mfft->fft, dir, in_data, out_data);
        if (ret != 0)
        {
            return ret;
        }
        in_data  = (real*)in_data+mfft->dist;
        out_data = (real*)out_data+mfft->dist;
    }
    return 0;
}


void
gmx_many_fft_destroy(gmx_fft_t    fft)
{
    gmx_many_fft_t mfft = (gmx_many_fft_t)fft;
    if (mfft != NULL)
    {
        if (mfft->fft != NULL)
        {
            gmx_fft_destroy(mfft->fft);
        }
        free(mfft);
    }
}

#endif //not GMX_FFT_FFTW3

int gmx_fft_transpose_2d(t_complex *          in_data,
                         t_complex *          out_data,
                         int                  nx,
                         int                  ny)
{
    int        i, j, k, im, n, ncount, done1, done2;
    int        i1, i1c, i2, i2c, kmi, max;

    t_complex  tmp1, tmp2, tmp3;
    t_complex *data;

    /* Use 500 bytes on stack to indicate moves.
     * This is just for optimization, it does not limit any dimensions.
     */
    char          move[500];
    int           nmove = 500;

    if (nx < 2 || ny < 2)
    {
        if (in_data != out_data)
        {
            memcpy(out_data, in_data, sizeof(t_complex)*nx*ny);
        }
        return 0;
    }

    /* Out-of-place transposes are easy */
    if (in_data != out_data)
    {
        for (i = 0; i < nx; i++)
        {
            for (j = 0; j < ny; j++)
            {
                out_data[j*nx+i].re = in_data[i*ny+j].re;
                out_data[j*nx+i].im = in_data[i*ny+j].im;
            }
        }
        return 0;
    }

    /* In-place transform. in_data=out_data=data */
    data = in_data;

    if (nx == ny)
    {
        /* trivial case, just swap elements */
        for (i = 0; i < nx; i++)
        {
            for (j = i+1; j < nx; j++)
            {
                tmp1.re         = data[i*nx+j].re;
                tmp1.im         = data[i*nx+j].im;
                data[i*nx+j].re = data[j*nx+i].re;
                data[i*nx+j].im = data[j*nx+i].im;
                data[j*nx+i].re = tmp1.re;
                data[j*nx+i].im = tmp1.im;
            }
        }
        return 0;
    }

    for (i = 0; i < nmove; i++)
    {
        move[i] = 0;
    }

    ncount = 2;

    if (nx > 2 && ny > 2)
    {
        i = nx-1;
        j = ny-1;
        do
        {
            k = i % j;
            i = j;
            j = k;
        }
        while (k);
        ncount += i-1;
    }

    n  = nx*ny;
    k  = n - 1;
    i  = 1;
    im = ny;

    done1 = 0;
    do
    {
        i1      = i;
        kmi     = k-i;
        tmp1.re = data[i1].re;
        tmp1.im = data[i1].im;
        i1c     = kmi;
        tmp2.re = data[i1c].re;
        tmp2.im = data[i1c].im;

        done2 = 0;
        do
        {
            i2  = ny*i1-k*(i1/nx);
            i2c = k-i2;
            if (i1 < nmove)
            {
                move[i1] = 1;
            }
            if (i1c < nmove)
            {
                move[i1c] = 1;
            }
            ncount += 2;
            if (i2 == i)
            {
                done2 = 1;
            }
            else if (i2 == kmi)
            {
                tmp3.re = tmp1.re;
                tmp3.im = tmp1.im;
                tmp1.re = tmp2.re;
                tmp1.im = tmp2.im;
                tmp2.re = tmp3.re;
                tmp2.im = tmp3.im;
                done2   = 1;
            }
            else
            {
                data[i1].re  = data[i2].re;
                data[i1].im  = data[i2].im;
                data[i1c].re = data[i2c].re;
                data[i1c].im = data[i2c].im;
                i1           = i2;
                i1c          = i2c;
            }
        }
        while (!done2);

        data[i1].re  = tmp1.re;
        data[i1].im  = tmp1.im;
        data[i1c].re = tmp2.re;
        data[i1c].im = tmp2.im;

        if (ncount >= n)
        {
            done1 = 1;
        }
        else
        {
            done2 = 0;
            do
            {
                max = k-i;
                i++;
                im += ny;
                if (im > k)
                {
                    im -= k;
                }
                i2 = im;
                if (i != i2)
                {
                    if (i >= nmove)
                    {
                        while (i2 > i && i2 < max)
                        {
                            i1 = i2;
                            i2 = ny*i1-k*(i1/nx);
                        }
                        if (i2 == i)
                        {
                            done2 = 1;
                        }
                    }
                    else if (!move[i])
                    {
                        done2 = 1;
                    }
                }
            }
            while (!done2);
        }
    }
    while (!done1);

    return 0;
}
