/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2003 David van der Spoel, Erik Lindahl, University of Groningen.
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

#include <errno.h>
#include <stdlib.h>

#include <mkl_dfti.h>
#include <mkl_service.h>

#include "gromacs/fft/fft.h"
#include "gromacs/utility/fatalerror.h"


/* For MKL version (<10.0), we should define MKL_LONG. */
#ifndef MKL_LONG
#define MKL_LONG long int
#endif


#ifdef GMX_DOUBLE
#define GMX_DFTI_PREC  DFTI_DOUBLE
#else
#define GMX_DFTI_PREC  DFTI_SINGLE
#endif

/*! \internal
 * \brief
 * Contents of the Intel MKL FFT fft datatype.
 *
 * Note that this is one of several possible implementations of gmx_fft_t.
 *
 *  The MKL _API_ supports 1D,2D, and 3D transforms, including real-to-complex.
 *  Unfortunately the actual library implementation does not support 3D real
 *  transforms as of version 7.2, and versions before 7.0 don't support 2D real
 *  either. In addition, the multi-dimensional storage format for real data
 *  is not compatible with our padding.
 *
 *  To work around this we roll our own 2D and 3D real-to-complex transforms,
 *  using separate X/Y/Z handles defined to perform (ny*nz), (nx*nz), and
 *  (nx*ny) transforms at once when necessary. To perform strided multiple
 *  transforms out-of-place (i.e., without padding in the last dimension)
 *  on the fly we also need to separate the forward and backward
 *  handles for real-to-complex/complex-to-real data permutation.
 *
 *  This makes it necessary to define 3 handles for in-place FFTs, and 4 for
 *  the out-of-place transforms. Still, whenever possible we try to use
 *  a single 3D-transform handle instead.
 *
 *  So, the handles are enumerated as follows:
 *
 *  1D FFT (real too):    Index 0 is the handle for the entire FFT
 *  2D complex FFT:       Index 0 is the handle for the entire FFT
 *  3D complex FFT:       Index 0 is the handle for the entire FFT
 *  2D, inplace real FFT: 0=FFTx, 1=FFTy handle
 *  2D, ooplace real FFT: 0=FFTx, 1=real-to-complex FFTy, 2=complex-to-real FFTy
 *  3D, inplace real FFT: 0=FFTx, 1=FFTy, 2=FFTz handle
 *  3D, ooplace real FFT: 0=FFTx, 1=FFTy, 2=r2c FFTz, 3=c2r FFTz
 *
 *  Intel people reading this: Learn from FFTW what a good interface looks like :-)
 */
#ifdef DOXYGEN
struct gmx_fft_mkl
#else
struct gmx_fft
#endif
{
    int                ndim;          /**< Number of dimensions in FFT  */
    int                nx;            /**< Length of X transform        */
    int                ny;            /**< Length of Y transform        */
    int                nz;            /**< Length of Z transform        */
    int                real_fft;      /**< 1 if real FFT, otherwise 0   */
    DFTI_DESCRIPTOR *  inplace[3];    /**< in-place FFT                 */
    DFTI_DESCRIPTOR *  ooplace[4];    /**< out-of-place FFT             */
    t_complex     *    work;          /**< Enable out-of-place c2r FFT  */
};



int
gmx_fft_init_1d(gmx_fft_t *        pfft,
                int                nx,
                gmx_fft_flag       flags)
{
    gmx_fft_t      fft;
    int            d;
    int            status;

    if (pfft == NULL)
    {
        gmx_fatal(FARGS, "Invalid opaque FFT datatype pointer.");
        return EINVAL;
    }
    *pfft = NULL;

    if ( (fft = (gmx_fft_t)malloc(sizeof(struct gmx_fft))) == NULL)
    {
        return ENOMEM;
    }

    /* Mark all handles invalid */
    for (d = 0; d < 3; d++)
    {
        fft->inplace[d] = fft->ooplace[d] = NULL;
    }
    fft->ooplace[3] = NULL;


    status = DftiCreateDescriptor(&fft->inplace[0], GMX_DFTI_PREC, DFTI_COMPLEX, 1, (MKL_LONG)nx);

    if (status == 0)
    {
        status = DftiSetValue(fft->inplace[0], DFTI_PLACEMENT, DFTI_INPLACE);
    }

    if (status == 0)
    {
        status = DftiCommitDescriptor(fft->inplace[0]);
    }


    if (status == 0)
    {
        status = DftiCreateDescriptor(&fft->ooplace[0], GMX_DFTI_PREC, DFTI_COMPLEX, 1, (MKL_LONG)nx);
    }

    if (status == 0)
    {
        DftiSetValue(fft->ooplace[0], DFTI_PLACEMENT, DFTI_NOT_INPLACE);
    }

    if (status == 0)
    {
        DftiCommitDescriptor(fft->ooplace[0]);
    }


    if (status != 0)
    {
        gmx_fatal(FARGS, "Error initializing Intel MKL FFT; status=%d", status);
        gmx_fft_destroy(fft);
        return status;
    }

    fft->ndim     = 1;
    fft->nx       = nx;
    fft->real_fft = 0;
    fft->work     = NULL;

    *pfft = fft;
    return 0;
}



int
gmx_fft_init_1d_real(gmx_fft_t *        pfft,
                     int                nx,
                     gmx_fft_flag       flags)
{
    gmx_fft_t      fft;
    int            d;
    int            status;

    if (pfft == NULL)
    {
        gmx_fatal(FARGS, "Invalid opaque FFT datatype pointer.");
        return EINVAL;
    }
    *pfft = NULL;

    if ( (fft = (gmx_fft_t)malloc(sizeof(struct gmx_fft))) == NULL)
    {
        return ENOMEM;
    }

    /* Mark all handles invalid */
    for (d = 0; d < 3; d++)
    {
        fft->inplace[d] = fft->ooplace[d] = NULL;
    }
    fft->ooplace[3] = NULL;

    status = DftiCreateDescriptor(&fft->inplace[0], GMX_DFTI_PREC, DFTI_REAL, 1, (MKL_LONG)nx);

    if (status == 0)
    {
        status = DftiSetValue(fft->inplace[0], DFTI_PLACEMENT, DFTI_INPLACE);
    }

    if (status == 0)
    {
        status = DftiCommitDescriptor(fft->inplace[0]);
    }


    if (status == 0)
    {
        status = DftiCreateDescriptor(&fft->ooplace[0], GMX_DFTI_PREC, DFTI_REAL, 1, (MKL_LONG)nx);
    }

    if (status == 0)
    {
        status = DftiSetValue(fft->ooplace[0], DFTI_PLACEMENT, DFTI_NOT_INPLACE);
    }

    if (status == 0)
    {
        status = DftiCommitDescriptor(fft->ooplace[0]);
    }


    if (status == DFTI_UNIMPLEMENTED)
    {
        gmx_fatal(FARGS,
                  "The linked Intel MKL version (<6.0?) cannot do real FFTs.");
        gmx_fft_destroy(fft);
        return status;
    }


    if (status != 0)
    {
        gmx_fatal(FARGS, "Error initializing Intel MKL FFT; status=%d", status);
        gmx_fft_destroy(fft);
        return status;
    }

    fft->ndim     = 1;
    fft->nx       = nx;
    fft->real_fft = 1;
    fft->work     = NULL;

    *pfft = fft;
    return 0;
}



int
gmx_fft_init_2d_real(gmx_fft_t *        pfft,
                     int                nx,
                     int                ny,
                     gmx_fft_flag       flags)
{
    gmx_fft_t      fft;
    int            d;
    int            status;
    MKL_LONG       stride[2];
    MKL_LONG       nyc;

    if (pfft == NULL)
    {
        gmx_fatal(FARGS, "Invalid opaque FFT datatype pointer.");
        return EINVAL;
    }
    *pfft = NULL;

    if ( (fft = (gmx_fft_t)malloc(sizeof(struct gmx_fft))) == NULL)
    {
        return ENOMEM;
    }

    nyc = (ny/2 + 1);

    /* Mark all handles invalid */
    for (d = 0; d < 3; d++)
    {
        fft->inplace[d] = fft->ooplace[d] = NULL;
    }
    fft->ooplace[3] = NULL;

    /* Roll our own 2D real transform using multiple transforms in MKL,
     * since the current MKL versions does not support our storage format,
     * and all but the most recent don't even have 2D real FFTs.
     */

    /* In-place X FFT */
    status = DftiCreateDescriptor(&fft->inplace[0], GMX_DFTI_PREC, DFTI_COMPLEX, 1, (MKL_LONG)nx);

    if (status == 0)
    {
        stride[0]  = 0;
        stride[1]  = nyc;

        status =
            (DftiSetValue(fft->inplace[0], DFTI_PLACEMENT, DFTI_INPLACE)    ||
             DftiSetValue(fft->inplace[0], DFTI_NUMBER_OF_TRANSFORMS, nyc)  ||
             DftiSetValue(fft->inplace[0], DFTI_INPUT_DISTANCE, 1)          ||
             DftiSetValue(fft->inplace[0], DFTI_INPUT_STRIDES, stride)      ||
             DftiSetValue(fft->inplace[0], DFTI_OUTPUT_DISTANCE, 1)         ||
             DftiSetValue(fft->inplace[0], DFTI_OUTPUT_STRIDES, stride));
    }

    if (status == 0)
    {
        status = DftiCommitDescriptor(fft->inplace[0]);
    }

    /* Out-of-place X FFT */
    if (status == 0)
    {
        status = DftiCreateDescriptor(&(fft->ooplace[0]), GMX_DFTI_PREC, DFTI_COMPLEX, 1, (MKL_LONG)nx);
    }

    if (status == 0)
    {
        stride[0] = 0;
        stride[1] = nyc;

        status =
            (DftiSetValue(fft->ooplace[0], DFTI_PLACEMENT, DFTI_NOT_INPLACE) ||
             DftiSetValue(fft->ooplace[0], DFTI_NUMBER_OF_TRANSFORMS, nyc)   ||
             DftiSetValue(fft->ooplace[0], DFTI_INPUT_DISTANCE, 1)           ||
             DftiSetValue(fft->ooplace[0], DFTI_INPUT_STRIDES, stride)       ||
             DftiSetValue(fft->ooplace[0], DFTI_OUTPUT_DISTANCE, 1)          ||
             DftiSetValue(fft->ooplace[0], DFTI_OUTPUT_STRIDES, stride));
    }

    if (status == 0)
    {
        status = DftiCommitDescriptor(fft->ooplace[0]);
    }


    /* In-place Y FFT  */
    if (status == 0)
    {
        status = DftiCreateDescriptor(&fft->inplace[1], GMX_DFTI_PREC, DFTI_REAL, 1, (MKL_LONG)ny);
    }

    if (status == 0)
    {
        stride[0] = 0;
        stride[1] = 1;

        status =
            (DftiSetValue(fft->inplace[1], DFTI_PLACEMENT, DFTI_INPLACE)             ||
             DftiSetValue(fft->inplace[1], DFTI_NUMBER_OF_TRANSFORMS, (MKL_LONG)nx)  ||
             DftiSetValue(fft->inplace[1], DFTI_INPUT_DISTANCE, 2*nyc)               ||
             DftiSetValue(fft->inplace[1], DFTI_INPUT_STRIDES, stride)               ||
             DftiSetValue(fft->inplace[1], DFTI_OUTPUT_DISTANCE, 2*nyc)              ||
             DftiSetValue(fft->inplace[1], DFTI_OUTPUT_STRIDES, stride)              ||
             DftiCommitDescriptor(fft->inplace[1]));
    }


    /* Out-of-place real-to-complex (affects output distance) Y FFT */
    if (status == 0)
    {
        status = DftiCreateDescriptor(&fft->ooplace[1], GMX_DFTI_PREC, DFTI_REAL, 1, (MKL_LONG)ny);
    }

    if (status == 0)
    {
        stride[0] = 0;
        stride[1] = 1;

        status =
            (DftiSetValue(fft->ooplace[1], DFTI_PLACEMENT, DFTI_NOT_INPLACE)           ||
             DftiSetValue(fft->ooplace[1], DFTI_NUMBER_OF_TRANSFORMS, (MKL_LONG)nx)    ||
             DftiSetValue(fft->ooplace[1], DFTI_INPUT_DISTANCE, (MKL_LONG)ny)          ||
             DftiSetValue(fft->ooplace[1], DFTI_INPUT_STRIDES, stride)                 ||
             DftiSetValue(fft->ooplace[1], DFTI_OUTPUT_DISTANCE, 2*nyc)                ||
             DftiSetValue(fft->ooplace[1], DFTI_OUTPUT_STRIDES, stride)                ||
             DftiCommitDescriptor(fft->ooplace[1]));
    }


    /* Out-of-place complex-to-real (affects output distance) Y FFT */
    if (status == 0)
    {
        status = DftiCreateDescriptor(&fft->ooplace[2], GMX_DFTI_PREC, DFTI_REAL, 1, (MKL_LONG)ny);
    }

    if (status == 0)
    {
        stride[0] = 0;
        stride[1] = 1;

        status =
            (DftiSetValue(fft->ooplace[2], DFTI_PLACEMENT, DFTI_NOT_INPLACE)           ||
             DftiSetValue(fft->ooplace[2], DFTI_NUMBER_OF_TRANSFORMS, (MKL_LONG)nx)    ||
             DftiSetValue(fft->ooplace[2], DFTI_INPUT_DISTANCE, 2*nyc)                 ||
             DftiSetValue(fft->ooplace[2], DFTI_INPUT_STRIDES, stride)                 ||
             DftiSetValue(fft->ooplace[2], DFTI_OUTPUT_DISTANCE, (MKL_LONG)ny)         ||
             DftiSetValue(fft->ooplace[2], DFTI_OUTPUT_STRIDES, stride)                ||
             DftiCommitDescriptor(fft->ooplace[2]));
    }


    if (status == 0)
    {
        if ((fft->work = (t_complex *)malloc(sizeof(t_complex)*(nx*(ny/2+1)))) == NULL)
        {
            status = ENOMEM;
        }
    }

    if (status != 0)
    {
        gmx_fatal(FARGS, "Error initializing Intel MKL FFT; status=%d", status);
        gmx_fft_destroy(fft);
        return status;
    }

    fft->ndim     = 2;
    fft->nx       = nx;
    fft->ny       = ny;
    fft->real_fft = 1;

    *pfft = fft;
    return 0;
}

int
gmx_fft_1d(gmx_fft_t                  fft,
           enum gmx_fft_direction     dir,
           void *                     in_data,
           void *                     out_data)
{
    int inplace = (in_data == out_data);
    int status  = 0;

    if ( (fft->real_fft == 1) || (fft->ndim != 1) ||
         ((dir != GMX_FFT_FORWARD) && (dir != GMX_FFT_BACKWARD)) )
    {
        gmx_fatal(FARGS, "FFT plan mismatch - bad plan or direction.");
        return EINVAL;
    }

    if (dir == GMX_FFT_FORWARD)
    {
        if (inplace)
        {
            status = DftiComputeForward(fft->inplace[0], in_data);
        }
        else
        {
            status = DftiComputeForward(fft->ooplace[0], in_data, out_data);
        }
    }
    else
    {
        if (inplace)
        {
            status = DftiComputeBackward(fft->inplace[0], in_data);
        }
        else
        {
            status = DftiComputeBackward(fft->ooplace[0], in_data, out_data);
        }
    }

    if (status != 0)
    {
        gmx_fatal(FARGS, "Error executing Intel MKL FFT.");
        status = -1;
    }

    return status;
}



int
gmx_fft_1d_real(gmx_fft_t                  fft,
                enum gmx_fft_direction     dir,
                void *                     in_data,
                void *                     out_data)
{
    int inplace = (in_data == out_data);
    int status  = 0;

    if ( (fft->real_fft != 1) || (fft->ndim != 1) ||
         ((dir != GMX_FFT_REAL_TO_COMPLEX) && (dir != GMX_FFT_COMPLEX_TO_REAL)) )
    {
        gmx_fatal(FARGS, "FFT plan mismatch - bad plan or direction.");
        return EINVAL;
    }

    if (dir == GMX_FFT_REAL_TO_COMPLEX)
    {
        if (inplace)
        {
            status = DftiComputeForward(fft->inplace[0], in_data);
        }
        else
        {
            status = DftiComputeForward(fft->ooplace[0], in_data, out_data);
        }
    }
    else
    {
        if (inplace)
        {
            status = DftiComputeBackward(fft->inplace[0], in_data);
        }
        else
        {
            status = DftiComputeBackward(fft->ooplace[0], in_data, out_data);
        }
    }

    if (status != 0)
    {
        gmx_fatal(FARGS, "Error executing Intel MKL FFT.");
        status = -1;
    }

    return status;
}


int
gmx_fft_2d_real(gmx_fft_t                  fft,
                enum gmx_fft_direction     dir,
                void *                     in_data,
                void *                     out_data)
{
    int inplace = (in_data == out_data);
    int status  = 0;

    if ( (fft->real_fft != 1) || (fft->ndim != 2) ||
         ((dir != GMX_FFT_REAL_TO_COMPLEX) && (dir != GMX_FFT_COMPLEX_TO_REAL)) )
    {
        gmx_fatal(FARGS, "FFT plan mismatch - bad plan or direction.");
        return EINVAL;
    }

    if (dir == GMX_FFT_REAL_TO_COMPLEX)
    {
        if (inplace)
        {
            /* real-to-complex in Y dimension, in-place */
            status = DftiComputeForward(fft->inplace[1], in_data);

            /* complex-to-complex in X dimension, in-place */
            if (status == 0)
            {
                status = DftiComputeForward(fft->inplace[0], in_data);
            }
        }
        else
        {
            /* real-to-complex in Y dimension, in_data to out_data */
            status = DftiComputeForward(fft->ooplace[1], in_data, out_data);

            /* complex-to-complex in X dimension, in-place to out_data */
            if (status == 0)
            {
                status = DftiComputeForward(fft->inplace[0], out_data);
            }
        }
    }
    else
    {
        /* prior implementation was incorrect. See fft.cpp unit test */
        gmx_incons("Complex -> Real is not supported by MKL.");
    }

    if (status != 0)
    {
        gmx_fatal(FARGS, "Error executing Intel MKL FFT.");
        status = -1;
    }

    return status;
}

void
gmx_fft_destroy(gmx_fft_t    fft)
{
    int d;

    if (fft != NULL)
    {
        for (d = 0; d < 3; d++)
        {
            if (fft->inplace[d] != NULL)
            {
                DftiFreeDescriptor(&fft->inplace[d]);
            }
            if (fft->ooplace[d] != NULL)
            {
                DftiFreeDescriptor(&fft->ooplace[d]);
            }
        }
        if (fft->ooplace[3] != NULL)
        {
            DftiFreeDescriptor(&fft->ooplace[3]);
        }
        if (fft->work != NULL)
        {
            free(fft->work);
        }
        free(fft);
    }
}

void gmx_fft_cleanup()
{
    mkl_free_buffers();
}

const char *gmx_fft_get_version_info()
{
    return "Intel MKL";
}
