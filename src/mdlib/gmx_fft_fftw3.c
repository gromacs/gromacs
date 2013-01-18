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

#ifdef GMX_FFT_FFTW3

#include <errno.h>
#include <stdlib.h>

#include <fftw3.h>

#include "gmx_fft.h"
#include "gmx_fatal.h"


#ifdef GMX_DOUBLE
#define FFTWPREFIX(name) fftw_ ## name
#else
#define FFTWPREFIX(name) fftwf_ ## name
#endif

#ifdef GMX_THREAD_MPI
#include "thread_mpi/threads.h"
#endif


#ifdef GMX_THREAD_MPI
/* none of the fftw3 calls, except execute(), are thread-safe, so
   we need to serialize them with this mutex. */
static tMPI_Thread_mutex_t big_fftw_mutex              = TMPI_THREAD_MUTEX_INITIALIZER;
static gmx_bool            gmx_fft_threads_initialized = FALSE;
#define FFTW_LOCK tMPI_Thread_mutex_lock(&big_fftw_mutex)
#define FFTW_UNLOCK tMPI_Thread_mutex_unlock(&big_fftw_mutex)
#else /* GMX_THREAD_MPI */
#define FFTW_LOCK
#define FFTW_UNLOCK
#endif /* GMX_THREAD_MPI */

/* We assume here that aligned memory starts at multiple of 16 bytes and unaligned memory starts at multiple of 8 bytes. The later is guranteed for all malloc implementation.
   Consequesences:
   - It is not allowed to use these FFT plans from memory which doesn't have a starting address as a multiple of 8 bytes.
     This is OK as long as the memory directly comes from malloc and is not some subarray within alloated memory.
   - This has to be fixed if any future architecute requires memory to be aligned to multiples of 32 bytes.
 */

struct gmx_fft
{
    /* Three alternatives (unaligned/aligned, out-of-place/in-place, forward/backward)
     * results in 8 different FFTW plans. Keep track of them with 3 array indices:
     * first index:   0=unaligned, 1=aligned
     * second index:  0=out-of-place, 1=in-place
     * third index:   0=backward, 1=forward
     */
    FFTWPREFIX(plan)         plan[2][2][2];
    /* Catch user mistakes */
    int                      real_transform;
    int                      ndim;
};



int
gmx_fft_init_1d(gmx_fft_t *        pfft,
                int                nx,
                gmx_fft_flag       flags)
{
    return gmx_fft_init_many_1d(pfft, nx, 1, flags);
}


int
gmx_fft_init_many_1d(gmx_fft_t *        pfft,
                     int                nx,
                     int                howmany,
                     gmx_fft_flag       flags)
{
    gmx_fft_t              fft;
    FFTWPREFIX(complex)   *p1, *p2, *up1, *up2;
    size_t                 pc;
    int                    i, j, k;
    int                    fftw_flags;

#ifdef GMX_DISABLE_FFTW_MEASURE
    flags |= GMX_FFT_FLAG_CONSERVATIVE;
#endif

    fftw_flags = (flags & GMX_FFT_FLAG_CONSERVATIVE) ? FFTW_ESTIMATE : FFTW_MEASURE;

    if (pfft == NULL)
    {
        gmx_fatal(FARGS, "Invalid opaque FFT datatype pointer.");
        return EINVAL;
    }
    *pfft = NULL;

    FFTW_LOCK;
    if ( (fft = (gmx_fft_t)FFTWPREFIX(malloc)(sizeof(struct gmx_fft))) == NULL)
    {
        FFTW_UNLOCK;
        return ENOMEM;
    }

    /* allocate aligned, and extra memory to make it unaligned */
    p1  = (FFTWPREFIX(complex) *) FFTWPREFIX(malloc)(sizeof(FFTWPREFIX(complex))*(nx+2)*howmany);
    if (p1 == NULL)
    {
        FFTWPREFIX(free)(fft);
        FFTW_UNLOCK;
        return ENOMEM;
    }

    p2  = (FFTWPREFIX(complex) *) FFTWPREFIX(malloc)(sizeof(FFTWPREFIX(complex))*(nx+2)*howmany);
    if (p2 == NULL)
    {
        FFTWPREFIX(free)(p1);
        FFTWPREFIX(free)(fft);
        FFTW_UNLOCK;
        return ENOMEM;
    }

    /* make unaligned pointers.
     * In double precision the actual complex datatype will be 16 bytes,
     * so go to a char pointer and force an offset of 8 bytes instead.
     */
    pc  = (size_t)p1;
    pc += 8;
    up1 = (FFTWPREFIX(complex) *)pc;

    pc  = (size_t)p2;
    pc += 8;
    up2 = (FFTWPREFIX(complex) *)pc;

    /*                            int rank, const int *n, int howmany,
                                  fftw_complex *in, const int *inembed,
                                  int istride, int idist,
                                  fftw_complex *out, const int *onembed,
                                  int ostride, int odist,
                                  int sign, unsigned flags */
    fft->plan[0][0][0] = FFTWPREFIX(plan_many_dft)(1, &nx, howmany, up1, &nx, 1, nx, up2, &nx, 1, nx, FFTW_BACKWARD, fftw_flags);
    fft->plan[0][0][1] = FFTWPREFIX(plan_many_dft)(1, &nx, howmany, up1, &nx, 1, nx, up2, &nx, 1, nx, FFTW_FORWARD, fftw_flags);
    fft->plan[0][1][0] = FFTWPREFIX(plan_many_dft)(1, &nx, howmany, up1, &nx, 1, nx, up1, &nx, 1, nx, FFTW_BACKWARD, fftw_flags);
    fft->plan[0][1][1] = FFTWPREFIX(plan_many_dft)(1, &nx, howmany, up1, &nx, 1, nx, up1, &nx, 1, nx, FFTW_FORWARD, fftw_flags);
    fft->plan[1][0][0] = FFTWPREFIX(plan_many_dft)(1, &nx, howmany, p1, &nx, 1, nx, p2, &nx, 1, nx, FFTW_BACKWARD, fftw_flags);
    fft->plan[1][0][1] = FFTWPREFIX(plan_many_dft)(1, &nx, howmany, p1, &nx, 1, nx, p2, &nx, 1, nx, FFTW_FORWARD, fftw_flags);
    fft->plan[1][1][0] = FFTWPREFIX(plan_many_dft)(1, &nx, howmany, p1, &nx, 1, nx, p1, &nx, 1, nx, FFTW_BACKWARD, fftw_flags);
    fft->plan[1][1][1] = FFTWPREFIX(plan_many_dft)(1, &nx, howmany, p1, &nx, 1, nx, p1, &nx, 1, nx, FFTW_FORWARD, fftw_flags);

    for (i = 0; i < 2; i++)
    {
        for (j = 0; j < 2; j++)
        {
            for (k = 0; k < 2; k++)
            {
                if (fft->plan[i][j][k] == NULL)
                {
                    gmx_fatal(FARGS, "Error initializing FFTW3 plan.");
                    FFTW_UNLOCK;
                    gmx_fft_destroy(fft);
                    FFTW_LOCK;
                    FFTWPREFIX(free)(p1);
                    FFTWPREFIX(free)(p2);
                    FFTW_UNLOCK;
                    return -1;
                }
            }
        }
    }

    FFTWPREFIX(free)(p1);
    FFTWPREFIX(free)(p2);

    fft->real_transform = 0;
    fft->ndim           = 1;

    *pfft = fft;
    FFTW_UNLOCK;
    return 0;
}


int
gmx_fft_init_1d_real(gmx_fft_t *        pfft,
                     int                nx,
                     gmx_fft_flag       flags)
{
    return gmx_fft_init_many_1d_real(pfft, nx, 1, flags);
}

int
gmx_fft_init_many_1d_real(gmx_fft_t *        pfft,
                          int                nx,
                          int                howmany,
                          gmx_fft_flag       flags)
{
    gmx_fft_t              fft;
    real                  *p1, *p2, *up1, *up2;
    size_t                 pc;
    int                    i, j, k;
    int                    fftw_flags;

#ifdef GMX_DISABLE_FFTW_MEASURE
    flags |= GMX_FFT_FLAG_CONSERVATIVE;
#endif

    fftw_flags = (flags & GMX_FFT_FLAG_CONSERVATIVE) ? FFTW_ESTIMATE : FFTW_MEASURE;

    if (pfft == NULL)
    {
        gmx_fatal(FARGS, "Invalid opaque FFT datatype pointer.");
        return EINVAL;
    }
    *pfft = NULL;

    FFTW_LOCK;
    if ( (fft = (gmx_fft_t) FFTWPREFIX(malloc)(sizeof(struct gmx_fft))) == NULL)
    {
        FFTW_UNLOCK;
        return ENOMEM;
    }

    /* allocate aligned, and extra memory to make it unaligned */
    p1  = (real *) FFTWPREFIX(malloc)(sizeof(real)*(nx/2+1)*2*howmany + 8);
    if (p1 == NULL)
    {
        FFTWPREFIX(free)(fft);
        FFTW_UNLOCK;
        return ENOMEM;
    }

    p2  = (real *) FFTWPREFIX(malloc)(sizeof(real)*(nx/2+1)*2*howmany + 8);
    if (p2 == NULL)
    {
        FFTWPREFIX(free)(p1);
        FFTWPREFIX(free)(fft);
        FFTW_UNLOCK;
        return ENOMEM;
    }

    /* make unaligned pointers.
     * In double precision the actual complex datatype will be 16 bytes,
     * so go to a char pointer and force an offset of 8 bytes instead.
     */
    pc  = (size_t)p1;
    pc += 8;
    up1 = (real *)pc;

    pc  = (size_t)p2;
    pc += 8;
    up2 = (real *)pc;

    /*                                int rank, const int *n, int howmany,
                                      double *in, const int *inembed,
                                      int istride, int idist,
                                      fftw_complex *out, const int *onembed,
                                      int ostride, int odist,
                                      unsigned flag    */
    fft->plan[0][0][1] = FFTWPREFIX(plan_many_dft_r2c)(1, &nx, howmany, up1, 0, 1, (nx/2+1) *2, (FFTWPREFIX(complex) *) up2, 0, 1, (nx/2+1), fftw_flags);
    fft->plan[0][1][1] = FFTWPREFIX(plan_many_dft_r2c)(1, &nx, howmany, up1, 0, 1, (nx/2+1) *2, (FFTWPREFIX(complex) *) up1, 0, 1, (nx/2+1), fftw_flags);
    fft->plan[1][0][1] = FFTWPREFIX(plan_many_dft_r2c)(1, &nx, howmany, p1, 0, 1, (nx/2+1) *2, (FFTWPREFIX(complex) *) p2, 0, 1, (nx/2+1), fftw_flags);
    fft->plan[1][1][1] = FFTWPREFIX(plan_many_dft_r2c)(1, &nx, howmany, p1, 0, 1, (nx/2+1) *2, (FFTWPREFIX(complex) *) p1, 0, 1, (nx/2+1), fftw_flags);

    fft->plan[0][0][0] = FFTWPREFIX(plan_many_dft_c2r)(1, &nx, howmany, (FFTWPREFIX(complex) *) up1, 0, 1, (nx/2+1), up2, 0, 1, (nx/2+1) *2, fftw_flags);
    fft->plan[0][1][0] = FFTWPREFIX(plan_many_dft_c2r)(1, &nx, howmany, (FFTWPREFIX(complex) *) up1, 0, 1, (nx/2+1), up1, 0, 1, (nx/2+1) *2, fftw_flags);
    fft->plan[1][0][0] = FFTWPREFIX(plan_many_dft_c2r)(1, &nx, howmany, (FFTWPREFIX(complex) *) p1, 0, 1, (nx/2+1), p2, 0, 1, (nx/2+1) *2, fftw_flags);
    fft->plan[1][1][0] = FFTWPREFIX(plan_many_dft_c2r)(1, &nx, howmany, (FFTWPREFIX(complex) *) p1, 0, 1, (nx/2+1), p1, 0, 1, (nx/2+1) *2, fftw_flags);

    for (i = 0; i < 2; i++)
    {
        for (j = 0; j < 2; j++)
        {
            for (k = 0; k < 2; k++)
            {
                if (fft->plan[i][j][k] == NULL)
                {
                    gmx_fatal(FARGS, "Error initializing FFTW3 plan.");
                    FFTW_UNLOCK;
                    gmx_fft_destroy(fft);
                    FFTW_LOCK;
                    FFTWPREFIX(free)(p1);
                    FFTWPREFIX(free)(p2);
                    FFTW_UNLOCK;
                    return -1;
                }
            }
        }
    }

    FFTWPREFIX(free)(p1);
    FFTWPREFIX(free)(p2);

    fft->real_transform = 1;
    fft->ndim           = 1;

    *pfft = fft;
    FFTW_UNLOCK;
    return 0;
}



int
gmx_fft_init_2d(gmx_fft_t *        pfft,
                int                nx,
                int                ny,
                gmx_fft_flag       flags)
{
    gmx_fft_t              fft;
    FFTWPREFIX(complex)   *p1, *p2, *up1, *up2;
    size_t                 pc;
    int                    i, j, k;
    int                    fftw_flags;

#ifdef GMX_DISABLE_FFTW_MEASURE
    flags |= GMX_FFT_FLAG_CONSERVATIVE;
#endif

    fftw_flags = (flags & GMX_FFT_FLAG_CONSERVATIVE) ? FFTW_ESTIMATE : FFTW_MEASURE;

    if (pfft == NULL)
    {
        gmx_fatal(FARGS, "Invalid opaque FFT datatype pointer.");
        return EINVAL;
    }
    *pfft = NULL;

    FFTW_LOCK;
    if ( (fft = (gmx_fft_t) FFTWPREFIX(malloc)(sizeof(struct gmx_fft))) == NULL)
    {
        FFTW_UNLOCK;
        return ENOMEM;
    }

    /* allocate aligned, and extra memory to make it unaligned */
    p1  = (FFTWPREFIX(complex) *) FFTWPREFIX(malloc)(sizeof(FFTWPREFIX(complex)) *(nx*ny+2));
    if (p1 == NULL)
    {
        FFTWPREFIX(free)(fft);
        FFTW_UNLOCK;
        return ENOMEM;
    }

    p2  = (FFTWPREFIX(complex) *) FFTWPREFIX(malloc)(sizeof(FFTWPREFIX(complex)) *(nx*ny+2));
    if (p2 == NULL)
    {
        FFTWPREFIX(free)(p1);
        FFTWPREFIX(free)(fft);
        FFTW_UNLOCK;
        return ENOMEM;
    }

    /* make unaligned pointers.
     * In double precision the actual complex datatype will be 16 bytes,
     * so go to a char pointer and force an offset of 8 bytes instead.
     */
    pc  = (size_t)p1;
    pc += 8;
    up1 = (FFTWPREFIX(complex) *)pc;

    pc  = (size_t)p2;
    pc += 8;
    up2 = (FFTWPREFIX(complex) *)pc;


    fft->plan[0][0][0] = FFTWPREFIX(plan_dft_2d)(nx, ny, up1, up2, FFTW_BACKWARD, fftw_flags);
    fft->plan[0][0][1] = FFTWPREFIX(plan_dft_2d)(nx, ny, up1, up2, FFTW_FORWARD, fftw_flags);
    fft->plan[0][1][0] = FFTWPREFIX(plan_dft_2d)(nx, ny, up1, up1, FFTW_BACKWARD, fftw_flags);
    fft->plan[0][1][1] = FFTWPREFIX(plan_dft_2d)(nx, ny, up1, up1, FFTW_FORWARD, fftw_flags);

    fft->plan[1][0][0] = FFTWPREFIX(plan_dft_2d)(nx, ny, p1, p2, FFTW_BACKWARD, fftw_flags);
    fft->plan[1][0][1] = FFTWPREFIX(plan_dft_2d)(nx, ny, p1, p2, FFTW_FORWARD, fftw_flags);
    fft->plan[1][1][0] = FFTWPREFIX(plan_dft_2d)(nx, ny, p1, p1, FFTW_BACKWARD, fftw_flags);
    fft->plan[1][1][1] = FFTWPREFIX(plan_dft_2d)(nx, ny, p1, p1, FFTW_FORWARD, fftw_flags);


    for (i = 0; i < 2; i++)
    {
        for (j = 0; j < 2; j++)
        {
            for (k = 0; k < 2; k++)
            {
                if (fft->plan[i][j][k] == NULL)
                {
                    gmx_fatal(FARGS, "Error initializing FFTW3 plan.");
                    FFTW_UNLOCK;
                    gmx_fft_destroy(fft);
                    FFTW_LOCK;
                    FFTWPREFIX(free)(p1);
                    FFTWPREFIX(free)(p2);
                    FFTW_UNLOCK;
                    return -1;
                }
            }
        }
    }

    FFTWPREFIX(free)(p1);
    FFTWPREFIX(free)(p2);

    fft->real_transform = 0;
    fft->ndim           = 2;

    *pfft = fft;
    FFTW_UNLOCK;
    return 0;
}



int
gmx_fft_init_2d_real(gmx_fft_t *        pfft,
                     int                nx,
                     int                ny,
                     gmx_fft_flag       flags)
{
    gmx_fft_t              fft;
    real                  *p1, *p2, *up1, *up2;
    size_t                 pc;
    int                    i, j, k;
    int                    fftw_flags;

#ifdef GMX_DISABLE_FFTW_MEASURE
    flags |= GMX_FFT_FLAG_CONSERVATIVE;
#endif

    fftw_flags = (flags & GMX_FFT_FLAG_CONSERVATIVE) ? FFTW_ESTIMATE : FFTW_MEASURE;

    if (pfft == NULL)
    {
        gmx_fatal(FARGS, "Invalid opaque FFT datatype pointer.");
        return EINVAL;
    }
    *pfft = NULL;

    FFTW_LOCK;
    if ( (fft = (gmx_fft_t) FFTWPREFIX(malloc)(sizeof(struct gmx_fft))) == NULL)
    {
        FFTW_UNLOCK;
        return ENOMEM;
    }

    /* allocate aligned, and extra memory to make it unaligned */
    p1  = (real *) FFTWPREFIX(malloc)(sizeof(real) *( nx*(ny/2+1)*2 + 2) );
    if (p1 == NULL)
    {
        FFTWPREFIX(free)(fft);
        FFTW_UNLOCK;
        return ENOMEM;
    }

    p2  = (real *) FFTWPREFIX(malloc)(sizeof(real) *( nx*(ny/2+1)*2 + 2) );
    if (p2 == NULL)
    {
        FFTWPREFIX(free)(p1);
        FFTWPREFIX(free)(fft);
        FFTW_UNLOCK;
        return ENOMEM;
    }

    /* make unaligned pointers.
     * In double precision the actual complex datatype will be 16 bytes,
     * so go to a char pointer and force an offset of 8 bytes instead.
     */
    pc  = (size_t)p1;
    pc += 8;
    up1 = (real *)pc;

    pc  = (size_t)p2;
    pc += 8;
    up2 = (real *)pc;


    fft->plan[0][0][0] = FFTWPREFIX(plan_dft_c2r_2d)(nx, ny, (FFTWPREFIX(complex) *) up1, up2, fftw_flags);
    fft->plan[0][0][1] = FFTWPREFIX(plan_dft_r2c_2d)(nx, ny, up1, (FFTWPREFIX(complex) *) up2, fftw_flags);
    fft->plan[0][1][0] = FFTWPREFIX(plan_dft_c2r_2d)(nx, ny, (FFTWPREFIX(complex) *) up1, up1, fftw_flags);
    fft->plan[0][1][1] = FFTWPREFIX(plan_dft_r2c_2d)(nx, ny, up1, (FFTWPREFIX(complex) *) up1, fftw_flags);

    fft->plan[1][0][0] = FFTWPREFIX(plan_dft_c2r_2d)(nx, ny, (FFTWPREFIX(complex) *) p1, p2, fftw_flags);
    fft->plan[1][0][1] = FFTWPREFIX(plan_dft_r2c_2d)(nx, ny, p1, (FFTWPREFIX(complex) *) p2, fftw_flags);
    fft->plan[1][1][0] = FFTWPREFIX(plan_dft_c2r_2d)(nx, ny, (FFTWPREFIX(complex) *) p1, p1, fftw_flags);
    fft->plan[1][1][1] = FFTWPREFIX(plan_dft_r2c_2d)(nx, ny, p1, (FFTWPREFIX(complex) *) p1, fftw_flags);


    for (i = 0; i < 2; i++)
    {
        for (j = 0; j < 2; j++)
        {
            for (k = 0; k < 2; k++)
            {
                if (fft->plan[i][j][k] == NULL)
                {
                    gmx_fatal(FARGS, "Error initializing FFTW3 plan.");
                    FFTW_UNLOCK;
                    gmx_fft_destroy(fft);
                    FFTW_LOCK;
                    FFTWPREFIX(free)(p1);
                    FFTWPREFIX(free)(p2);
                    FFTW_UNLOCK;
                    return -1;
                }
            }
        }
    }

    FFTWPREFIX(free)(p1);
    FFTWPREFIX(free)(p2);

    fft->real_transform = 1;
    fft->ndim           = 2;

    *pfft = fft;
    FFTW_UNLOCK;
    return 0;
}



int
gmx_fft_init_3d(gmx_fft_t *        pfft,
                int                nx,
                int                ny,
                int                nz,
                gmx_fft_flag       flags)
{
    gmx_fft_t              fft;
    FFTWPREFIX(complex)   *p1, *p2, *up1, *up2;
    size_t                 pc;
    int                    i, j, k;
    int                    fftw_flags;

#ifdef GMX_DISABLE_FFTW_MEASURE
    flags |= GMX_FFT_FLAG_CONSERVATIVE;
#endif

    fftw_flags = (flags & GMX_FFT_FLAG_CONSERVATIVE) ? FFTW_ESTIMATE : FFTW_MEASURE;

    if (pfft == NULL)
    {
        gmx_fatal(FARGS, "Invalid opaque FFT datatype pointer.");
        return EINVAL;
    }
    *pfft = NULL;

    FFTW_LOCK;
    if ( (fft = (gmx_fft_t) FFTWPREFIX(malloc)(sizeof(struct gmx_fft))) == NULL)
    {
        FFTW_UNLOCK;
        return ENOMEM;
    }

    /* allocate aligned, and extra memory to make it unaligned */
    p1  = (FFTWPREFIX(complex) *) FFTWPREFIX(malloc)(sizeof(FFTWPREFIX(complex)) *(nx*ny*nz+2));
    if (p1 == NULL)
    {
        FFTWPREFIX(free)(fft);
        FFTW_UNLOCK;
        return ENOMEM;
    }

    p2  = (FFTWPREFIX(complex) *) FFTWPREFIX(malloc)(sizeof(FFTWPREFIX(complex)) *(nx*ny*nz+2));
    if (p2 == NULL)
    {
        FFTWPREFIX(free)(p1);
        FFTWPREFIX(free)(fft);
        FFTW_UNLOCK;
        return ENOMEM;
    }

    /* make unaligned pointers.
     * In double precision the actual complex datatype will be 16 bytes,
     * so go to a char pointer and force an offset of 8 bytes instead.
     */
    pc  = (size_t)p1;
    pc += 8;
    up1 = (FFTWPREFIX(complex) *)pc;

    pc  = (size_t)p2;
    pc += 8;
    up2 = (FFTWPREFIX(complex) *)pc;


    fft->plan[0][0][0] = FFTWPREFIX(plan_dft_3d)(nx, ny, nz, up1, up2, FFTW_BACKWARD, fftw_flags);
    fft->plan[0][0][1] = FFTWPREFIX(plan_dft_3d)(nx, ny, nz, up1, up2, FFTW_FORWARD, fftw_flags);
    fft->plan[0][1][0] = FFTWPREFIX(plan_dft_3d)(nx, ny, nz, up1, up1, FFTW_BACKWARD, fftw_flags);
    fft->plan[0][1][1] = FFTWPREFIX(plan_dft_3d)(nx, ny, nz, up1, up1, FFTW_FORWARD, fftw_flags);

    fft->plan[1][0][0] = FFTWPREFIX(plan_dft_3d)(nx, ny, nz, p1, p2, FFTW_BACKWARD, fftw_flags);
    fft->plan[1][0][1] = FFTWPREFIX(plan_dft_3d)(nx, ny, nz, p1, p2, FFTW_FORWARD, fftw_flags);
    fft->plan[1][1][0] = FFTWPREFIX(plan_dft_3d)(nx, ny, nz, p1, p1, FFTW_BACKWARD, fftw_flags);
    fft->plan[1][1][1] = FFTWPREFIX(plan_dft_3d)(nx, ny, nz, p1, p1, FFTW_FORWARD, fftw_flags);


    for (i = 0; i < 2; i++)
    {
        for (j = 0; j < 2; j++)
        {
            for (k = 0; k < 2; k++)
            {
                if (fft->plan[i][j][k] == NULL)
                {
                    gmx_fatal(FARGS, "Error initializing FFTW3 plan.");
                    FFTW_UNLOCK;
                    gmx_fft_destroy(fft);
                    FFTW_LOCK;
                    FFTWPREFIX(free)(p1);
                    FFTWPREFIX(free)(p2);
                    FFTW_UNLOCK;
                    return -1;
                }
            }
        }
    }

    FFTWPREFIX(free)(p1);
    FFTWPREFIX(free)(p2);

    fft->real_transform = 0;
    fft->ndim           = 3;

    *pfft = fft;
    FFTW_UNLOCK;
    return 0;
}



int
gmx_fft_init_3d_real(gmx_fft_t *        pfft,
                     int                nx,
                     int                ny,
                     int                nz,
                     gmx_fft_flag       flags)
{
    gmx_fft_t              fft;
    real                  *p1, *p2, *up1, *up2;
    size_t                 pc;
    int                    i, j, k;
    int                    fftw_flags;

#ifdef GMX_DISABLE_FFTW_MEASURE
    flags |= GMX_FFT_FLAG_CONSERVATIVE;
#endif

    fftw_flags = (flags & GMX_FFT_FLAG_CONSERVATIVE) ? FFTW_ESTIMATE : FFTW_MEASURE;

    if (pfft == NULL)
    {
        gmx_fatal(FARGS, "Invalid opaque FFT datatype pointer.");
        return EINVAL;
    }
    *pfft = NULL;

    FFTW_LOCK;
    if ( (fft = (gmx_fft_t) FFTWPREFIX(malloc)(sizeof(struct gmx_fft))) == NULL)
    {
        FFTW_UNLOCK;
        return ENOMEM;
    }

    /* allocate aligned, and extra memory to make it unaligned */
    p1  = (real *) FFTWPREFIX(malloc)(sizeof(real) *( nx*ny*(nz/2+1)*2 + 2) );
    if (p1 == NULL)
    {
        FFTWPREFIX(free)(fft);
        FFTW_UNLOCK;
        return ENOMEM;
    }

    p2  = (real *) FFTWPREFIX(malloc)(sizeof(real) *( nx*ny*(nz/2+1)*2 + 2) );
    if (p2 == NULL)
    {
        FFTWPREFIX(free)(p1);
        FFTWPREFIX(free)(fft);
        FFTW_UNLOCK;
        return ENOMEM;
    }

    /* make unaligned pointers.
     * In double precision the actual complex datatype will be 16 bytes,
     * so go to a void pointer and force an offset of 8 bytes instead.
     */
    pc  = (size_t)p1;
    pc += 8;
    up1 = (real *)pc;

    pc  = (size_t)p2;
    pc += 8;
    up2 = (real *)pc;


    fft->plan[0][0][0] = FFTWPREFIX(plan_dft_c2r_3d)(nx, ny, nz, (FFTWPREFIX(complex) *) up1, up2, fftw_flags);
    fft->plan[0][0][1] = FFTWPREFIX(plan_dft_r2c_3d)(nx, ny, nz, up1, (FFTWPREFIX(complex) *) up2, fftw_flags);
    fft->plan[0][1][0] = FFTWPREFIX(plan_dft_c2r_3d)(nx, ny, nz, (FFTWPREFIX(complex) *) up1, up1, fftw_flags);
    fft->plan[0][1][1] = FFTWPREFIX(plan_dft_r2c_3d)(nx, ny, nz, up1, (FFTWPREFIX(complex) *) up1, fftw_flags);

    fft->plan[1][0][0] = FFTWPREFIX(plan_dft_c2r_3d)(nx, ny, nz, (FFTWPREFIX(complex) *) p1, p2, fftw_flags);
    fft->plan[1][0][1] = FFTWPREFIX(plan_dft_r2c_3d)(nx, ny, nz, p1, (FFTWPREFIX(complex) *) p2, fftw_flags);
    fft->plan[1][1][0] = FFTWPREFIX(plan_dft_c2r_3d)(nx, ny, nz, (FFTWPREFIX(complex) *) p1, p1, fftw_flags);
    fft->plan[1][1][1] = FFTWPREFIX(plan_dft_r2c_3d)(nx, ny, nz, p1, (FFTWPREFIX(complex) *) p1, fftw_flags);


    for (i = 0; i < 2; i++)
    {
        for (j = 0; j < 2; j++)
        {
            for (k = 0; k < 2; k++)
            {
                if (fft->plan[i][j][k] == NULL)
                {
                    gmx_fatal(FARGS, "Error initializing FFTW3 plan.");
                    FFTW_UNLOCK;
                    gmx_fft_destroy(fft);
                    FFTW_LOCK;
                    FFTWPREFIX(free)(p1);
                    FFTWPREFIX(free)(p2);
                    FFTW_UNLOCK;
                    return -1;
                }
            }
        }
    }

    FFTWPREFIX(free)(p1);
    FFTWPREFIX(free)(p2);

    fft->real_transform = 1;
    fft->ndim           = 3;

    *pfft = fft;
    FFTW_UNLOCK;
    return 0;
}


int
gmx_fft_1d               (gmx_fft_t                  fft,
                          enum gmx_fft_direction     dir,
                          void *                     in_data,
                          void *                     out_data)
{
    int           aligned   = ((((size_t)in_data | (size_t)out_data) & 0xf) == 0);
    int           inplace   = (in_data == out_data);
    int           isforward = (dir == GMX_FFT_FORWARD);

    /* Some checks */
    if ( (fft->real_transform == 1) || (fft->ndim != 1) ||
         ((dir != GMX_FFT_FORWARD) && (dir != GMX_FFT_BACKWARD)) )
    {
        gmx_fatal(FARGS, "FFT plan mismatch - bad plan or direction.");
        return EINVAL;
    }

    FFTWPREFIX(execute_dft)(fft->plan[aligned][inplace][isforward],
                            (FFTWPREFIX(complex) *) in_data,
                            (FFTWPREFIX(complex) *) out_data);

    return 0;
}

int
gmx_fft_many_1d               (gmx_fft_t                  fft,
                               enum gmx_fft_direction     dir,
                               void *                     in_data,
                               void *                     out_data)
{
    return gmx_fft_1d(fft, dir, in_data, out_data);
}

int
gmx_fft_1d_real          (gmx_fft_t                  fft,
                          enum gmx_fft_direction     dir,
                          void *                     in_data,
                          void *                     out_data)
{
    int           aligned   = ((((size_t)in_data | (size_t)out_data) & 0xf) == 0);
    int           inplace   = (in_data == out_data);
    int           isforward = (dir == GMX_FFT_REAL_TO_COMPLEX);

    /* Some checks */
    if ( (fft->real_transform != 1) || (fft->ndim != 1) ||
         ((dir != GMX_FFT_REAL_TO_COMPLEX) && (dir != GMX_FFT_COMPLEX_TO_REAL)) )
    {
        gmx_fatal(FARGS, "FFT plan mismatch - bad plan or direction.");
        return EINVAL;
    }

    if (isforward)
    {
        FFTWPREFIX(execute_dft_r2c)(fft->plan[aligned][inplace][isforward],
                                    (real *)in_data, (FFTWPREFIX(complex) *) out_data);
    }
    else
    {
        FFTWPREFIX(execute_dft_c2r)(fft->plan[aligned][inplace][isforward],
                                    (FFTWPREFIX(complex) *) in_data, (real *)out_data);
    }

    return 0;
}

int
gmx_fft_many_1d_real     (gmx_fft_t                  fft,
                          enum gmx_fft_direction     dir,
                          void *                     in_data,
                          void *                     out_data)
{
    return gmx_fft_1d_real(fft, dir, in_data, out_data);
}

int
gmx_fft_2d               (gmx_fft_t                  fft,
                          enum gmx_fft_direction     dir,
                          void *                     in_data,
                          void *                     out_data)
{
    int           aligned   = ((((size_t)in_data | (size_t)out_data) & 0xf) == 0);
    int           inplace   = (in_data == out_data);
    int           isforward = (dir == GMX_FFT_FORWARD);

    /* Some checks */
    if ( (fft->real_transform == 1) || (fft->ndim != 2) ||
         ((dir != GMX_FFT_FORWARD) && (dir != GMX_FFT_BACKWARD)) )
    {
        gmx_fatal(FARGS, "FFT plan mismatch - bad plan or direction.");
        return EINVAL;
    }

    FFTWPREFIX(execute_dft)(fft->plan[aligned][inplace][isforward],
                            (FFTWPREFIX(complex) *) in_data,
                            (FFTWPREFIX(complex) *) out_data);

    return 0;
}


int
gmx_fft_2d_real          (gmx_fft_t                  fft,
                          enum gmx_fft_direction     dir,
                          void *                     in_data,
                          void *                     out_data)
{
    int           aligned   = ((((size_t)in_data | (size_t)out_data) & 0xf) == 0);
    int           inplace   = (in_data == out_data);
    int           isforward = (dir == GMX_FFT_REAL_TO_COMPLEX);

    /* Some checks */
    if ( (fft->real_transform != 1) || (fft->ndim != 2) ||
         ((dir != GMX_FFT_REAL_TO_COMPLEX) && (dir != GMX_FFT_COMPLEX_TO_REAL)) )
    {
        gmx_fatal(FARGS, "FFT plan mismatch - bad plan or direction.");
        return EINVAL;
    }

    if (isforward)
    {
        FFTWPREFIX(execute_dft_r2c)(fft->plan[aligned][inplace][isforward],
                                    (real *)in_data,
                                    (FFTWPREFIX(complex) *) out_data);
    }
    else
    {
        FFTWPREFIX(execute_dft_c2r)(fft->plan[aligned][inplace][isforward],
                                    (FFTWPREFIX(complex) *) in_data,
                                    (real *)out_data);
    }


    return 0;
}


int
gmx_fft_3d               (gmx_fft_t                  fft,
                          enum gmx_fft_direction     dir,
                          void *                     in_data,
                          void *                     out_data)
{
    int           aligned   = ((((size_t)in_data | (size_t)out_data) & 0xf) == 0);
    int           inplace   = (in_data == out_data);
    int           isforward = (dir == GMX_FFT_FORWARD);

    /* Some checks */
    if ( (fft->real_transform == 1) || (fft->ndim != 3) ||
         ((dir != GMX_FFT_FORWARD) && (dir != GMX_FFT_BACKWARD)) )
    {
        gmx_fatal(FARGS, "FFT plan mismatch - bad plan or direction.");
        return EINVAL;
    }

    FFTWPREFIX(execute_dft)(fft->plan[aligned][inplace][isforward],
                            (FFTWPREFIX(complex) *) in_data,
                            (FFTWPREFIX(complex) *) out_data);

    return 0;
}


int
gmx_fft_3d_real          (gmx_fft_t                  fft,
                          enum gmx_fft_direction     dir,
                          void *                     in_data,
                          void *                     out_data)
{
    int           aligned   = ((((size_t)in_data | (size_t)out_data) & 0xf) == 0);
    int           inplace   = (in_data == out_data);
    int           isforward = (dir == GMX_FFT_REAL_TO_COMPLEX);

    /* Some checks */
    if ( (fft->real_transform != 1) || (fft->ndim != 3) ||
         ((dir != GMX_FFT_REAL_TO_COMPLEX) && (dir != GMX_FFT_COMPLEX_TO_REAL)) )
    {
        gmx_fatal(FARGS, "FFT plan mismatch - bad plan or direction.");
        return EINVAL;
    }

    if (isforward)
    {
        FFTWPREFIX(execute_dft_r2c)(fft->plan[aligned][inplace][isforward],
                                    (real *)in_data,
                                    (FFTWPREFIX(complex) *) out_data);
    }
    else
    {
        FFTWPREFIX(execute_dft_c2r)(fft->plan[aligned][inplace][isforward],
                                    (FFTWPREFIX(complex) *) in_data,
                                    (real *)out_data);
    }


    return 0;
}


void
gmx_fft_destroy(gmx_fft_t      fft)
{
    int                   i, j, k;

    if (fft != NULL)
    {
        for (i = 0; i < 2; i++)
        {
            for (j = 0; j < 2; j++)
            {
                for (k = 0; k < 2; k++)
                {
                    if (fft->plan[i][j][k] != NULL)
                    {
                        FFTW_LOCK;
                        FFTWPREFIX(destroy_plan)(fft->plan[i][j][k]);
                        FFTW_UNLOCK;
                        fft->plan[i][j][k] = NULL;
                    }
                }
            }
        }
        FFTW_LOCK;
        FFTWPREFIX(free)(fft);
        FFTW_UNLOCK;
    }

}

void
gmx_many_fft_destroy(gmx_fft_t    fft)
{
    gmx_fft_destroy(fft);
}

#else
int
    gmx_fft_fftw3_empty;
#endif /* GMX_FFT_FFTW3 */
