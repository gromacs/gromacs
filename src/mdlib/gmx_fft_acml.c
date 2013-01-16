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

#ifdef GMX_FFT_ACML

#include <errno.h>
#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include <assert.h>

#include "acml.h"


#include "gmx_fft.h"
#include "gmx_fatal.h"

/*	Since c has no support of templates, we get around the type restrictions
    with macros.  Our custom macro names obscure the float vs double interfaces
 */
#ifdef GMX_DOUBLE
#define ACML_FFT1DX zfft1dx
#define ACML_FFT2DX zfft2dx
#define ACML_FFT3DX zfft3dy
#define ACML_FFT1DMX zfft1mx

#define ACML_RCFFT1D dzfft
#define ACML_CRFFT1D zdfft
#define ACML_RCFFT1DM dzfftm
#define ACML_CRFFT1DM zdfftm

#define acmlComplex doublecomplex
#else

#define ACML_FFT1DX cfft1dx
#define ACML_FFT2DX cfft2dx
#define ACML_FFT3DX cfft3dy
#define ACML_FFT1DMX cfft1mx

#define ACML_RCFFT1D scfft
#define ACML_CRFFT1D csfft
#define ACML_RCFFT1DM scfftm
#define ACML_CRFFT1DM csfftm

#define acmlComplex complex
#endif

void cPrintArray( t_complex* arr, int outer, int mid, int inner )
{
    int i, j, k;

    printf("\n");
    for (i = 0; i < outer; ++i)
    {
        for (j = 0; j < mid; ++j)
        {
            for (k = 0; k < inner; ++k)
            {
                printf("%f,%f ", arr[i*inner*mid+j*inner+k].re, arr[i*inner*mid+j*inner+k].im );
            }
            printf("\n");
        }
        printf("\n");
    }

    return;
}

void rPrintArray( real* arr, int outer, int mid, int inner )
{
    int i, j, k;

    printf("\n");
    for (i = 0; i < outer; ++i)
    {
        for (j = 0; j < mid; ++j)
        {
            for (k = 0; k < inner; ++k)
            {
                printf("%f ", arr[i*inner*mid+j*inner+k] );
            }
            printf("\n");
        }
        printf("\n");
    }

    return;
}

/* Contents of the AMD ACML FFT fft datatype.
 *
 * Note that this is one of several possible implementations of gmx_fft_t.
 *
 *  The ACML _API_ supports 1D,2D, and 3D transforms, including real-to-complex.
 *  Unfortunately the actual library implementation does not support 2D&3D real
 *  transforms as of version 4.2  In addition, ACML outputs their hermitian data
 *  in such a form that it can't be fed straight back into their complex API's.
 *
 *  To work around this we roll our own 2D and 3D real-to-complex transforms,
 *  using separate X/Y/Z handles defined to perform (ny*nz), (nx*nz), and
 *  (nx*ny) transforms at once when necessary. To perform strided multiple
 *  transforms out-of-place (i.e., without padding in the last dimension)
 *  on the fly we also need to separate the forward and backward
 *  handles for real-to-complex/complex-to-real data permutation.
 *
 *  So, the handles are enumerated as follows:
 *
 *  1D FFT (real too):    Index 0 is the handle for the entire FFT
 *  2D complex FFT:       Index 0 is the handle for the entire FFT
 *  3D complex FFT:       Index 0 is the handle for the entire FFT
 *  2D, real FFT:			0=FFTx, 1=FFTrc, 2=FFTcr handle
 *  3D, real FFT:			0=FFTx, 1=FFTy, 2=FFTrc, 3=FFTcr handle
 *
 *  AMD people reading this: Learn from FFTW what a good interface looks like :-)
 *
 */
struct gmx_fft
{
    /* Work arrays;
     * c2c = index0
     * r2c = index0, c2r = index1
     */
    void       *       comm[4];
    void       *       realScratch;

    int                ndim;              /**< Number of dimensions in FFT  */
    int                nx;                /**< Length of X transform        */
    int                ny;                /**< Length of Y transform        */
    int                nz;                /**< Length of Z transform        */
    int                real_fft;          /**< 1 if real FFT, otherwise 0   */
};

/*	ACML's real FFT support leaves the data in a swizzled 'hermitian' format.  The
    problem with this is that you can't just feed this data back into ACML :(  A
    pre-processing step is required to transform the hermitian swizzled complex
    values into unswizzled complex values
    -This function assumes that the complex conjugates are not expanded; the
    calling function should not assume that they exist.
    -This function pads all real numbers with 0's for imaginary components.
 */
void    hermitianUnPacking( float scale, int row, int col, void* src, int srcPitch, void* dst )
{
    int mid         = col/2;
    int loopCount   = (col-1)/2;
    int hermLength  = mid + 1;

    int nFFT;

    /* Currently this function is not written to handle aliased input/output pointers */
    assert( src != dst );

    /* These two pointers could be aliased on top of each other; be careful */
    real    *   realData        = (real*)src;
    t_complex*  hermData        = (t_complex*)dst;

    /*	We have to expand the data from real to complex, which will clobber data if you start
     *	from the beginning.  For this reason, we start at the last array and work backwards,
     *	however, we skip the very first row (nFFT == 0) because we still have data collision
     *	on that row.  We copy the first row into a temporary buffer.
     */
    for (nFFT = row-1; nFFT >= 0; --nFFT)
    {
        realData        = (real*)src + (nFFT*srcPitch);
        hermData        = (t_complex*)dst + (nFFT*hermLength);

        /* The first complex number is real valued */
        t_complex   tmp = { realData[0]*scale, 0 };
        hermData[0] = tmp;

        int i;
        for (i = 1; i <= loopCount; ++i)
        {
            t_complex   tmp = { realData[i]*scale, realData[col-i]*scale };
            hermData[i] = tmp;
        }

        /* A little cleanup for even lengths */
        if (loopCount != mid)
        {
            /* The n/2 number is real valued */
            t_complex   tmp = { realData[mid]*scale, 0 };
            hermData[mid] = tmp;
        }
    }

    return;
}

/*	ACML's real FFT support requires the data in a swizzled 'hermitian' format.
    This is a non-standard format, and requires us to manually swizzle the data :(
    A pre-processing step is required to transform the complex values into the
    swizzled hermitian format
 */
void    hermitianPacking( float scale, int row, int col, void* src, void* dst, int dstPitch )
{
    int         mid         = col/2;
    int         loopCount   = (col-1)/2;
    int         hermLength  = mid + 1;

    int         nFFT;

    real    *   realData;
    t_complex*  hermData;

    /* Currently this function is not written to handle aliased input/output pointers */
    assert( src != dst );

    for (nFFT = 0; nFFT < row; ++nFFT)
    {
        realData        = (real*)dst + (nFFT*dstPitch);
        hermData        = (t_complex*)src + (nFFT*hermLength);

        /* The first complex number is real valued */
        realData[0] = hermData[0].re * scale;

        /* ACML's complex to real function is documented as only handling 'forward'
         * transforms, which mean we have to manually conjugate the data before doing
         * the backwards transform. */
        int i;
        for (i = 1; i <= loopCount; ++i)
        {
            realData[i]     = hermData[i].re * scale;
            realData[col-i] = -hermData[i].im * scale;
        }

        /* A little cleanup for even lengths */
        if (loopCount != mid)
        {
            /* The n/2 number is real valued */
            realData[mid] = hermData[mid].re * scale;
        }

    }

    return;
}

/*	Gromacs expects results that are returned to be in a packed form; No
    complex conjugate's.  This routine takes a 2D array of complex data,
    and compresses it to fit into 'real' form by removing all the complex
    conjugates and any known imaginary 0's.  Data is expected in row-major
    form, with ny describing the length of a row.
 */
void    compressConj2D( int nX, int nY, void* src, void* dst )
{
    /* These two pointers are aliased on top of each other; be careful */
    real    *   realData    = (real*)dst;
    t_complex*  complexData;

    int         halfRows    = nX/2 + 1;
    int         halfCols    = nY/2 + 1;
    int         rowOdd      = nX & 0x1;
    int         colOdd      = nY & 0x1;


    int nRow;
    for (nRow = 0; nRow < nX; ++nRow)
    {
        complexData = (t_complex*)src + (nRow*halfCols);

        int nCol;
        for (nCol = 0; nCol < halfCols; ++nCol)
        {
            if (( nRow == 0) && (nCol == 0 ))
            {
                /* The complex number is real valued */
                realData[0] = complexData[nCol].re;
                realData   += 1;
            }
            /* Last column is real if even column length */
            else if ( (nRow == 0) && ( !colOdd && ( nCol == (halfCols-1) ) ) )
            {
                /* The complex number is real valued */
                realData[0] = complexData[nCol].re;
                realData   += 1;
            }
            /* The middle value of the first column is real if we have an even number of rows and
             * the last column of middle row is real if we have an even number of rows and columns */
            else if (!rowOdd && ( ( nRow == (halfRows-1) ) && ( ( nCol == 0 ) || ( !colOdd && ( nCol == (halfCols-1) ) ) ) ) )
            {
                /* The complex number is real valued */
                realData[0] = complexData[nCol].re;
                realData   += 1;
            }
            else if ( (nCol == 0) || ( !colOdd && ( nCol == (halfCols-1) ) ) )
            {
                /* Last half of real columns are conjugates */
                if (nRow < halfRows)
                {
                    /* Copy a complex number, which is two reals */
                    realData[0] = complexData[nCol].re;
                    realData[1] = complexData[nCol].im;
                    realData   += 2;
                }
            }
            else
            {
                /* Copy a complex number, which is two reals */
                realData[0] = complexData[nCol].re;
                realData[1] = complexData[nCol].im;
                realData   += 2;
            }
        }
    }

    return;
}

/*	Gromacs expects results that are returned to be in a packed form; No
    complex conjugate's.  This routine takes a 2D array of real data,
    and uncompresses it to fit into 'complex' form by creating all the complex
    conjugates and any known imaginary 0's.  Data is expected in row-major
    form, with ny describing the length of a row.
 */
void    unCompressConj2D( int nX, int nY, void* src, void* dst )
{
    /* These two pointers are aliased on top of each other; be careful */
    real*       realData    = (real*)src;
    t_complex*  complexData = (t_complex*)dst;

    int         halfRows    = nX/2 + 1;
    int         halfCols    = nY/2 + 1;
    int         rowOdd      = nX & 0x1;
    int         colOdd      = nY & 0x1;

    int         nRow;
    for (nRow = 0; nRow < nX; ++nRow)
    {
        int nCol;
        for (nCol = 0; nCol < halfCols; ++nCol)
        {
            int iIndex  = nRow*halfCols + nCol;

            if (( nRow == 0) && (nCol == 0 ))
            {
                /* The complex number is real valued */
                complexData[iIndex].re = realData[0];
                complexData[iIndex].im =  0;
                realData              += 1;
            }
            /* Last column is real if even column length */
            else if ( (nRow == 0) && ( !colOdd && ( nCol == (halfCols-1) ) ) )
            {
                /* The complex number is real valued */
                complexData[iIndex].re = realData[0];
                complexData[iIndex].im =  0;
                realData              += 1;
            }
            /* The middle value of the first column is real if we have an even number of rows and
             * the last column of middle row is real if we have an even number of rows and columns */
            else if (!rowOdd && ( ( nRow == (halfRows-1) ) && ( ( nCol == 0 ) || ( !colOdd && ( nCol == (halfCols-1) ) ) ) ) )
            {
                /* The complex number is real valued */
                complexData[iIndex].re = realData[0];
                complexData[iIndex].im =  0;
                realData              += 1;
            }
            else if ( (nCol == 0) || (!colOdd && ( nCol == (halfCols-1) ) ) )
            {
                /* Last half of real columns are conjugates */
                if (nRow < halfRows)
                {
                    /* Copy a complex number, which is two reals */
                    complexData[iIndex].re = realData[0];
                    complexData[iIndex].im = realData[1];
                    realData              += 2;
                }
                else
                {
                    int oIndex  = (nX-nRow)*halfCols + nCol;
                    complexData[iIndex].re = complexData[oIndex].re;
                    complexData[iIndex].im = -complexData[oIndex].im;
                }
            }
            else
            {
                /* Copy a complex number, which is two reals */
                complexData[iIndex].re = realData[0];
                complexData[iIndex].im = realData[1];
                realData              += 2;
            }
        }
    }

    return;
}

/* Support routine for the 1D array tranforms.  ACML does not support a MODE
 * flag on the real/complex interface.  It assumes a 'forward' transform both
 * directions, so that requires a manual conjugation of the imaginary comps. */
void    negateConj( int len, void* src )
{
    real* imag    = (real*)src;
    int   half    = len/2 + 1;
    int   i;

    for (i = half; i < len; ++i)
    {
        imag[i] = -imag[i];
    }

    return;
}

int
gmx_fft_init_1d(gmx_fft_t *        pfft,
                int                nx,
                enum gmx_fft_flag  flags)
{
    gmx_fft_t          fft;
    int                info     = 0;
    acmlComplex*       comm     = NULL;
    int                commSize = 0;


    if (pfft == NULL)
    {
        gmx_fatal(FARGS, "Invalid opaque FFT datatype pointer.");
        return EINVAL;
    }
    *pfft = NULL;

    if ( (fft = malloc(sizeof(struct gmx_fft))) == NULL)
    {
        return ENOMEM;
    }

    //	Single precision requires 5*nx+100
    //	Double precision requires 3*nx+100
    if (sizeof( acmlComplex ) == 16)
    {
        commSize    = (3*nx+100)*sizeof( acmlComplex );
    }
    else
    {
        commSize    = (5*nx+100)*sizeof( acmlComplex );
    }

    // Allocate communication work array
    if ( (comm = (acmlComplex*)malloc( commSize ) ) == NULL)
    {
        return ENOMEM;
    }

    // Initialize communication work array
    ACML_FFT1DX( 100, 1.0f, TRUE, nx, NULL, 1, NULL, 1, comm, &info );

    if (info != 0)
    {
        gmx_fatal(FARGS, "Error initializing ACML FFT; status=%d", info);
        gmx_fft_destroy( fft );
        return info;
    }

    fft->ndim     = 1;
    fft->nx       = nx;
    fft->real_fft = 0;
    fft->comm[0]  = comm;

    *pfft = fft;
    return 0;
}

int
gmx_fft_init_1d_real(gmx_fft_t *        pfft,
                     int                nx,
                     enum gmx_fft_flag  flags)
{
    gmx_fft_t      fft;
    int            info     = 0;
    real    *      commRC   = NULL;
    real    *      commCR   = NULL;
    int            commSize = 0;

    if (pfft == NULL)
    {
        gmx_fatal(FARGS, "Invalid opaque FFT datatype pointer.");
        return EINVAL;
    }
    *pfft = NULL;

    if ( (fft = malloc(sizeof(struct gmx_fft))) == NULL)
    {
        return ENOMEM;
    }


    commSize    = (3*nx+100)*sizeof( float );

    // Allocate communication work array, r2c
    if ( (commRC = (real*)malloc( commSize ) ) == NULL)
    {
        return ENOMEM;
    }

    // Allocate communication work array, c2r
    if ( (commCR = (real*)malloc( commSize ) ) == NULL)
    {
        return ENOMEM;
    }

    // Initialize communication work array
    ACML_RCFFT1D( 100, nx, NULL, (real*)commRC, &info );

    if (info != 0)
    {
        gmx_fatal(FARGS, "Error initializing ACML FFT; status=%d", info);
        gmx_fft_destroy( fft );
        return info;
    }

    // Initialize communication work array
    ACML_CRFFT1D( 100, nx, NULL, (real*)commCR, &info );

    if (info != 0)
    {
        gmx_fatal(FARGS, "Error initializing ACML FFT; status=%d", info);
        gmx_fft_destroy( fft );
        return info;
    }

    /* Allocate scratch work array that ACML uses to splat from hermitian complex format to
     * full complex format */
    if ( (fft->realScratch = (acmlComplex*)malloc( nx*sizeof( acmlComplex ) ) ) == NULL)
    {
        return ENOMEM;
    }

    fft->ndim     = 1;
    fft->nx       = nx;
    fft->real_fft = 1;
    fft->comm[0]  = commRC;
    fft->comm[1]  = commCR;

    *pfft = fft;
    return 0;
}

int
gmx_fft_init_2d(gmx_fft_t *        pfft,
                int                nx,
                int                ny,
                enum gmx_fft_flag  flags)
{
    gmx_fft_t          fft;
    int                info     = 0;
    acmlComplex*       comm     = NULL;
    int                commSize = 0;


    if (pfft == NULL)
    {
        gmx_fatal(FARGS, "Invalid opaque FFT datatype pointer.");
        return EINVAL;
    }
    *pfft = NULL;

    if ( (fft = malloc(sizeof(struct gmx_fft))) == NULL)
    {
        return ENOMEM;
    }

    //	Single precision requires nx*ny+5*(nx+ny)
    //	Double precision requires nx*ny+3*(nx+ny)
    if (sizeof( acmlComplex ) == 16)
    {
        commSize    = (nx*ny+3*(nx+ny)+200)*sizeof( acmlComplex );
    }
    else
    {
        commSize    = (nx*ny+5*(nx+ny)+200)*sizeof( acmlComplex );
    }

    // Allocate communication work array
    if ( (comm = (acmlComplex*)malloc( commSize ) ) == NULL)
    {
        return ENOMEM;
    }

    // Initialize communication work array
    ACML_FFT2DX( 100, 1.0f, TRUE, TRUE, nx, ny, NULL, 1, nx, NULL, 1, nx, (acmlComplex*)comm, &info );

    if (info != 0)
    {
        gmx_fatal(FARGS, "Error initializing ACML FFT; status=%d", info);
        gmx_fft_destroy( fft );
        return info;
    }

    fft->ndim     = 2;
    fft->nx       = nx;
    fft->ny       = ny;
    fft->real_fft = 0;
    fft->comm[0]  = comm;

    *pfft = fft;
    return 0;
}

int
gmx_fft_init_2d_real(gmx_fft_t *        pfft,
                     int                nx,
                     int                ny,
                     enum gmx_fft_flag  flags)
{
    gmx_fft_t      fft;
    int            info     = 0;
    acmlComplex*   comm     = NULL;
    real*          commRC   = NULL;
    real*          commCR   = NULL;
    int            commSize = 0;
    int            nyc      = 0;

    if (pfft == NULL)
    {
        gmx_fatal(FARGS, "Invalid opaque FFT datatype pointer.");
        return EINVAL;
    }
    *pfft = NULL;

    if ( (fft = malloc(sizeof(struct gmx_fft))) == NULL)
    {
        return ENOMEM;
    }

    nyc = (ny/2 + 1);

    /* Roll our own 2D real transform using multiple transforms in ACML,
     * since the current ACML versions does not support our storage format,
     * and all but the most recent don't even have 2D real FFTs.
     */

    //	Single precision requires 5*nx+100
    //	Double precision requires 3*nx+100
    if (sizeof( acmlComplex ) == 16)
    {
        commSize    = (3*nx+100)*sizeof( acmlComplex );
    }
    else
    {
        commSize    = (5*nx+100)*sizeof( acmlComplex );
    }

    // Allocate communication work array
    if ( (comm = (acmlComplex*)malloc( commSize ) ) == NULL)
    {
        return ENOMEM;
    }

    // Initialize communication work array
    ACML_FFT1DMX( 100, 1.0f, FALSE, nyc, nx, NULL, nyc, 1, NULL, nyc, 1, comm, &info );

    if (info != 0)
    {
        gmx_fatal(FARGS, "Error initializing ACML FFT; status=%d", info);
        gmx_fft_destroy( fft );
        return info;
    }

    commSize    = (3*ny+100)*sizeof( real );
    // Allocate communication work array
    if ( (commRC = (real*)malloc( commSize ) ) == NULL)
    {
        return ENOMEM;
    }

    //	TODO:  Is there no MODE or PLAN for multiple hermetian sequences?
    // Initialize communication work array
    ACML_RCFFT1D( 100, ny, NULL, commRC, &info );

    if (info != 0)
    {
        gmx_fatal(FARGS, "Error initializing ACML FFT; status=%d", info);
        gmx_fft_destroy( fft );
        return info;
    }

    commSize    = (3*ny+100)*sizeof( real );
    // Allocate communication work array
    if ( (commCR = (real*)malloc( commSize ) ) == NULL)
    {
        return ENOMEM;
    }

    //	TODO:  Is there no MODE or PLAN for multiple hermetian sequences?
    // Initialize communication work array
    ACML_CRFFT1D( 100, ny, NULL, commCR, &info );

    if (info != 0)
    {
        gmx_fatal(FARGS, "Error initializing ACML FFT; status=%d", info);
        gmx_fft_destroy( fft );
        return info;
    }

    /* Allocate scratch work array that ACML uses to splat from hermitian complex format to
     * full complex format */
    if ( (fft->realScratch = (acmlComplex*)malloc( (nx*ny)*sizeof( acmlComplex ) ) ) == NULL)
    {
        return ENOMEM;
    }

    fft->ndim     = 2;
    fft->nx       = nx;
    fft->ny       = ny;
    fft->real_fft = 1;
    fft->comm[0]  = comm;
    fft->comm[1]  = commRC;
    fft->comm[2]  = commCR;

    *pfft = fft;
    return 0;
}

int
gmx_fft_init_3d(gmx_fft_t *        pfft,
                int                nx,
                int                ny,
                int                nz,
                enum gmx_fft_flag  flags)
{
    gmx_fft_t          fft;
    int                info     = 0;
    acmlComplex*       comm     = NULL;
    int                commSize = 0;


    if (pfft == NULL)
    {
        gmx_fatal(FARGS, "Invalid opaque FFT datatype pointer.");
        return EINVAL;
    }
    *pfft = NULL;

    if ( (fft = malloc(sizeof(struct gmx_fft))) == NULL)
    {
        return ENOMEM;
    }

    commSize    = (nx*ny*nz+4*(nx+ny+nz)+300)*sizeof( acmlComplex );

    // Allocate communication work array
    if ( (comm = (acmlComplex*)malloc( commSize ) ) == NULL)
    {
        return ENOMEM;
    }

    ACML_FFT3DX( 100, 1.0f, TRUE, nx, ny, nz, NULL, 1, nx, nx*ny, NULL, 1, nx, nx*ny, comm, commSize, &info );

    fft->ndim     = 3;
    fft->nx       = nx;
    fft->ny       = ny;
    fft->nz       = nz;
    fft->real_fft = 0;
    fft->comm[0]  = comm;

    *pfft = fft;
    return 0;
}

int
gmx_fft_init_3d_real(gmx_fft_t *        pfft,
                     int                nx,
                     int                ny,
                     int                nz,
                     enum gmx_fft_flag  flags)
{
    gmx_fft_t      fft;
    int            info     = 0;
    acmlComplex*   commX    = NULL;
    acmlComplex*   commY    = NULL;
    real*          commRC   = NULL;
    real*          commCR   = NULL;
    int            commSize = 0;

    if (pfft == NULL)
    {
        gmx_fatal(FARGS, "Invalid opaque FFT datatype pointer.");
        return EINVAL;
    }
    *pfft = NULL;

    /* nzc = (nz/2 + 1); */

    if ( (fft = malloc(sizeof(struct gmx_fft))) == NULL)
    {
        return ENOMEM;
    }

    /* Roll our own 3D real transform using multiple transforms in ACML,
     * since the current ACML versions does not support 2D
     * or 3D real transforms.
     */

    /* In-place X FFT.
     * ny*nz complex-to-complex transforms, length nx
     * transform distance: 1
     * element strides: ny*nz
     */

    /*	Single precision requires 5*nx+100
        Double precision requires 3*nx+100
     */
    if (sizeof( acmlComplex ) == 16)
    {
        commSize    = (3*nx+100)*sizeof( acmlComplex );
    }
    else
    {
        commSize    = (5*nx+100)*sizeof( acmlComplex );
    }

    /* Allocate communication work array */
    if ( (commX = (acmlComplex*)malloc( commSize ) ) == NULL)
    {
        return ENOMEM;
    }

    /* Initialize communication work array */
    ACML_FFT1DMX( 100, 1.0f, TRUE, ny*nz, nx, NULL, ny*nz, 1, NULL, ny*nz, 1, commX, &info );

    if (info != 0)
    {
        gmx_fatal(FARGS, "Error initializing ACML FFT; status=%d", info);
        gmx_fft_destroy( fft );
        return info;
    }

    /* In-place Y FFT.
     * We cannot do all NX*NZ transforms at once, so define a handle to do
     * NZ transforms, and then execute it NX times.
     * nz complex-to-complex transforms, length ny
     * transform distance: 1
     * element strides: nz
     */
    /*	Single precision requires 5*nx+100
        Double precision requires 3*nx+100
     */
    if (sizeof( acmlComplex ) == 16)
    {
        commSize    = (3*ny+100)*sizeof( acmlComplex );
    }
    else
    {
        commSize    = (5*ny+100)*sizeof( acmlComplex );
    }

    /* Allocate communication work array */
    if ( (commY = (acmlComplex*)malloc( commSize ) ) == NULL)
    {
        return ENOMEM;
    }

    /* Initialize communication work array */
    /* We want to do multiple 1D FFT's in z-y plane, so we have to loop over x
     * dimension recalculating z-y plane for each slice.
     */
    ACML_FFT1DMX( 100, 1.0f, TRUE, nz, ny, NULL, nz, 1, NULL, nz, 1, commY, &info );

    if (info != 0)
    {
        gmx_fatal(FARGS, "Error initializing ACML FFT; status=%d", info);
        gmx_fft_destroy( fft );
        return info;
    }

    /* In-place Z FFT:
     * nx*ny real-to-complex transforms, length nz
     * transform distance: nzc*2 -> nzc*2
     * element strides: 1
     */

    commSize    = (3*nz+100)*sizeof( real );
    /* Allocate communication work array */
    if ( (commRC = (real*)malloc( commSize ) ) == NULL)
    {
        return ENOMEM;
    }

    /*	TODO:  Is there no MODE or PLAN for multiple hermetian sequences? */
    // Initialize communication work array
    ACML_RCFFT1D( 100, nz, NULL, commRC, &info );


    if (info != 0)
    {
        gmx_fatal(FARGS, "Error initializing ACML FFT; status=%d", info);
        gmx_fft_destroy( fft );
        return info;
    }

    /* Out-of-place complex-to-real (affects distance) Z FFT:
     * nx*ny real-to-complex transforms, length nz
     * transform distance: nzc*2 -> nz
     * element STRIDES: 1
     */
    commSize    = (3*nz+100)*sizeof( real );
    /* Allocate communication work array */
    if ( (commCR = (real*)malloc( commSize ) ) == NULL)
    {
        return ENOMEM;
    }

    // Initialize communication work array
    ACML_CRFFT1D( 100, nz, NULL, commCR, &info );

    if (info != 0)
    {
        gmx_fatal(FARGS, "Error initializing ACML FFT; status=%d", info);
        gmx_fft_destroy( fft );
        return info;
    }

    /* Allocate scratch work array that ACML uses to splat from hermitian complex format to
     * full complex format */
    if ( (fft->realScratch = (acmlComplex*)malloc( (nx*ny*nz)*sizeof( acmlComplex ) ) ) == NULL)
    {
        return ENOMEM;
    }

    fft->ndim     = 3;
    fft->nx       = nx;
    fft->ny       = ny;
    fft->nz       = nz;
    fft->real_fft = 1;
    fft->comm[0]  = commX;
    fft->comm[1]  = commY;
    fft->comm[2]  = commRC;
    fft->comm[3]  = commCR;

    *pfft = fft;
    return 0;
}

int
gmx_fft_1d(gmx_fft_t                  fft,
           enum gmx_fft_direction     dir,
           void *                     in_data,
           void *                     out_data)
{
    int inpl = (in_data == out_data);
    int info = 0;
    int mode = (dir == GMX_FFT_FORWARD) ? -1 : 1;

    if ( (fft->real_fft == 1) || (fft->ndim != 1) ||
         ((dir != GMX_FFT_FORWARD) && (dir != GMX_FFT_BACKWARD)) )
    {
        gmx_fatal(FARGS, "FFT plan mismatch - bad plan or direction.");
        return EINVAL;
    }

    ACML_FFT1DX( mode, 1.0f, inpl, fft->nx, in_data, 1, out_data, 1, fft->comm[0], &info );

    if (info != 0)
    {
        gmx_fatal(FARGS, "Error executing AMD ACML FFT.");
        info = -1;
    }

    return info;
}

int
gmx_fft_1d_real(gmx_fft_t                  fft,
                enum gmx_fft_direction     dir,
                void *                     inPtr,
                void *                     outPtr)
{
    int info = 0;
    int nx   = fft->nx;

    if ( (fft->real_fft != 1) || (fft->ndim != 1) ||
         ((dir != GMX_FFT_REAL_TO_COMPLEX) && (dir != GMX_FFT_COMPLEX_TO_REAL)) )
    {
        gmx_fatal(FARGS, "FFT plan mismatch - bad plan or direction.");
        return EINVAL;
    }

    /* I apply a correction scale to the automatic scaling that was done in the real-complex step */
    float   recipCorrection = sqrtf( nx );

    /*	ACML needs to do complex2complex intermediate transforms, which will not fit in the amount
     *  of memory allocated by the gromacs program, which assumes real.  The real2complex transforms
     *  are also only in-place, so we manually do a memcpy() first */
    if (dir == GMX_FFT_REAL_TO_COMPLEX)
    {
        memcpy( fft->realScratch, inPtr, nx*sizeof( real ) );
        ACML_RCFFT1D( 1, nx, fft->realScratch, fft->comm[0], &info );

        hermitianUnPacking( recipCorrection, 1, nx, fft->realScratch, 0, outPtr );
    }
    else
    {
        memcpy( fft->realScratch, inPtr, nx*sizeof( t_complex ) );
        hermitianPacking( recipCorrection, 1, nx, fft->realScratch, outPtr, 0 );

        ACML_CRFFT1D( 1, nx, outPtr, fft->comm[1], &info );
    }

    if (info != 0)
    {
        gmx_fatal(FARGS, "Error executing AMD ACML FFT.");
        info = -1;
    }

    return info;
}

int
gmx_fft_2d(gmx_fft_t                  fft,
           enum gmx_fft_direction     dir,
           void *                     in_data,
           void *                     out_data)
{
    int inpl = (in_data == out_data);
    int info = 0;
    int mode = (dir == GMX_FFT_FORWARD) ? -1 : 1;

    if ( (fft->real_fft == 1) || (fft->ndim != 2) ||
         ((dir != GMX_FFT_FORWARD) && (dir != GMX_FFT_BACKWARD)) )
    {
        gmx_fatal(FARGS, "FFT plan mismatch - bad plan or direction.");
        return EINVAL;
    }

    ACML_FFT2DX( mode, 1.0f, TRUE, inpl, fft->nx, fft->ny,
                 in_data, 1, fft->nx, out_data, 1, fft->nx, fft->comm[0], &info );

    if (info != 0)
    {
        gmx_fatal(FARGS, "Error executing AMD ACML FFT.");
        info = -1;
    }

    return info;
}


int
gmx_fft_2d_real(gmx_fft_t                  fft,
                enum gmx_fft_direction     dir,
                void *                     inPtr,
                void *                     outPtr )
{
    int inpl = (inPtr == outPtr);
    int info = 0;
    int i    = 0;

    int nx      = fft->nx;
    int ny      = fft->ny;
    int nyc     = (fft->ny/2 + 1);
    int nyLen   = 0;

    /* Depending on whether we are calculating the FFT in place or not, the rows of the
     * 2D FFT will be packed or not.
     */
    if (inpl)
    {
        /* Not packed; 1 or 2 extra reals padding at end of each row */
        nyLen   = nyc*2;
    }
    else
    {
        nyLen   = ny;
    }

    if ( (fft->real_fft != 1) || (fft->ndim != 2) ||
         ((dir != GMX_FFT_REAL_TO_COMPLEX) && (dir != GMX_FFT_COMPLEX_TO_REAL)) )
    {
        gmx_fatal(FARGS, "FFT plan mismatch - bad plan or direction.");
        return EINVAL;
    }

    /* I apply a correction scale to the automatic scaling that was done in the real-complex step */
    float   recipCorrection = sqrtf( ny );

    if (dir == GMX_FFT_REAL_TO_COMPLEX)
    {
        /*	ACML needs to do complex2complex intermediate transforms, which will not fit in the amount
         *  of memory allocated by the gromacs program, which assumes real.  The real2complex transforms
         *  are also only in-place, so we manually do a memcpy() first */
        if (!inpl)
        {
            memcpy( outPtr, inPtr, nx*ny*sizeof( real ) );
        }

        /* real-to-complex in Y dimension, in-place
         * SCFFTM is not valid to call here.  SCFFTM does not take any stride information, and assumes that
         * the rows are tightly packed.  GROMACS pads rows with either 1 or 2 extra reals depending
         * on even or odd lengths.
         */
        for (i = 0; i < nx; ++i)
        {
            if (info == 0)
            {
                ACML_RCFFT1D( 1, ny, (real*)outPtr+i*nyLen, fft->comm[1], &info );
            }
        }

        hermitianUnPacking( 1.0f, nx, ny, outPtr, nyLen, fft->realScratch );

        /* complex-to-complex in X dimension */
        if (info == 0)
        {
            ACML_FFT1DMX( -1, recipCorrection, FALSE, nyc, nx, fft->realScratch, nyc, 1,
                          outPtr, nyc, 1, fft->comm[0], &info );
        }
    }
    else
    {
        /* complex-to-complex in X dimension, in-place */
        ACML_FFT1DMX( 1, recipCorrection, FALSE, nyc, nx, inPtr, nyc, 1,
                      fft->realScratch, nyc, 1, fft->comm[0], &info );

        hermitianPacking( 1.0f, nx, ny, fft->realScratch, outPtr, nyLen );

        /* complex-to-real in Y dimension
         * CSFFTM is not valid to call here.  SCFFTM does not take any stride information, and assumes that
         * the rows are tightly packed.  GROMACS pads rows with either 1 or 2 extra reals depending
         * on even or odd lengths.
         */
        if (info == 0)
        {
            for (i = 0; i < nx; ++i)
            {
                if (info == 0)
                {
                    ACML_CRFFT1D( 1, ny, (real*)outPtr+i*nyLen, fft->comm[2], &info );
                }
            }
        }
    }

    if (info != 0)
    {
        gmx_fatal(FARGS, "Error executing AMD ACML FFT.");
        info = -1;
    }

    return info;
}


int
gmx_fft_3d(gmx_fft_t                  fft,
           enum gmx_fft_direction     dir,
           void *                     in_data,
           void *                     out_data)
{
    int mode = (dir == GMX_FFT_FORWARD) ? -1 : 1;
    int inpl = ( in_data == out_data );

    int commSize    = 0;
    int nx          = fft->nx;
    int ny          = fft->ny;
    int nz          = fft->nz;
    int info        = 0;

    if ( (fft->real_fft == 1) || (fft->ndim != 3) ||
         ((dir != GMX_FFT_FORWARD) && (dir != GMX_FFT_BACKWARD)) )
    {
        gmx_fatal(FARGS, "FFT plan mismatch - bad plan or direction.");
        return EINVAL;
    }

    commSize    = (nx*ny*nz+4*(nx+ny+nz)+300)*sizeof( acmlComplex );

    ACML_FFT3DX( mode, 1.0f, inpl, nx, ny, nz, in_data, 1, nx, nx*ny,
                 out_data, 1, nx, nx*ny, fft->comm[0], commSize, &info );

    if (info != 0)
    {
        gmx_fatal(FARGS, "Error executing AMD ACML FFT.");
        info = -1;
    }

    return info;
}

int
gmx_fft_3d_real(gmx_fft_t                  fft,
                enum gmx_fft_direction     dir,
                void *                     inPtr,
                void *                     outPtr)
{
    int inpl = (inPtr == outPtr);
    int info = 0;
    int i;
    int nx, ny, nz, nzLen = 0;

    nx  = fft->nx;
    ny  = fft->ny;
    nz  = fft->nz;
    int nzc = (nz/2 + 1);

    /* Depending on whether we are calculating the FFT in place or not, the rows of the
     * 3D FFT will be packed or not.
     */
    if (inpl)
    {
        /* Not packed; 1 or 2 extra reals padding at end of each row */
        nzLen   = nzc*2;
    }
    else
    {
        nzLen   = nz;
    }

    if ( (fft->real_fft != 1) || (fft->ndim != 3) ||
         ((dir != GMX_FFT_REAL_TO_COMPLEX) && (dir != GMX_FFT_COMPLEX_TO_REAL)) )
    {
        gmx_fatal(FARGS, "FFT plan mismatch - bad plan or direction.");
        return EINVAL;
    }

    /* I apply a correction scale to the automatic scaling that was done in the real-complex step */
    float   recipCorrection = sqrtf( nz );

    if (dir == GMX_FFT_REAL_TO_COMPLEX)
    {
        /*	ACML needs to do complex2complex intermediate transforms, which will not fit in the amount
         *  of memory allocated by the gromacs program, which assumes real.  The real2complex transforms
         *  are also only in-place, so we manually do a memcpy() first */
        if (!inpl)
        {
            memcpy( outPtr, inPtr, nx*ny*nz*sizeof( real ) );
        }

        /* real-to-complex in Z dimension, in-place
         * SCFFTM is not valid to call here.  SCFFTM does not take any stride information, and assumes that
         * the rows are tightly packed.  GROMACS pads rows with either 1 or 2 extra reals depending
         * on even or odd lengths.
         */
        for (i = 0; i < nx*ny; ++i)
        {
            if (info == 0)
            {
                ACML_RCFFT1D( 1, nz, (real*)outPtr+i*nzLen, fft->comm[2], &info );
            }
        }

        hermitianUnPacking( 1.0f, nx*ny, nz, outPtr, nzLen, fft->realScratch );

        /* complex-to-complex in Y dimension, in-place */
        for (i = 0; i < nx; i++)
        {
            if (info == 0)
            {
                ACML_FFT1DMX( -1, 1.0f, TRUE, nzc, ny, (acmlComplex*)fft->realScratch+i*ny*nzc, nzc, 1,
                              (acmlComplex*)fft->realScratch+i*ny*nzc, nzc, 1, fft->comm[1], &info );
            }
        }

        /* complex-to-complex in X dimension, in-place */
        if (info == 0)
        {
            ACML_FFT1DMX( -1, recipCorrection, FALSE, ny*nzc, nx, fft->realScratch, ny*nzc, 1, outPtr, ny*nzc, 1, fft->comm[0], &info );
        }
    }
    else
    {
        /* complex-to-complex in X dimension, from inPtr to work */
        ACML_FFT1DMX( 1, recipCorrection, FALSE, ny*nzc, nx, inPtr, ny*nzc, 1, fft->realScratch, ny*nzc, 1, fft->comm[0], &info );

        /* complex-to-complex in Y dimension, in-place */
        for (i = 0; i < nx; i++)
        {
            if (info == 0)
            {
                ACML_FFT1DMX( 1, 1.0f, TRUE, nzc, ny, (acmlComplex*)fft->realScratch+i*ny*nzc, nzc, 1,
                              (acmlComplex*)fft->realScratch+i*ny*nzc, nzc, 1, fft->comm[1], &info );
            }
        }

        hermitianPacking( 1.0f, nx*ny, nz, fft->realScratch, outPtr, nzLen );

        /* complex-to-real in Z dimension, in-place
         * CSFFTM is not valid to call here.  CSFFTM does not take any stride information, and assumes that
         * the rows are tightly packed.  GROMACS pads rows with either 1 or 2 extra reals depending
         * on even or odd lengths.
         */
        for (i = 0; i < nx*ny; ++i)
        {
            if (info == 0)
            {
                ACML_CRFFT1D( 1, nz, (real*)outPtr+i*nzLen, fft->comm[3], &info );
            }
        }

    }

    if (info != 0)
    {
        gmx_fatal(FARGS, "Error executing AMD ACML FFT.");
        info = -1;
    }

    return info;
}

void
gmx_fft_destroy(gmx_fft_t    fft)
{
    int d;

    if (fft != NULL)
    {
        for (d = 0; d < 4; d++)
        {
            if (fft->comm[d] != NULL)
            {
                free(fft->comm[d]);
            }
        }

        if (fft->realScratch != NULL)
        {
            free( fft->realScratch );
        }

        free( fft );
    }
}

#else
int
    gmx_fft_acml_empty;
#endif /* GMX_FFT_ACML */
