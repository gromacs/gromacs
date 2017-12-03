/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2009,2010,2012,2013,2014,2015,2016,2017, by the GROMACS development team, led by
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

#include "fft5d.h"

#include "config.h"

#include <assert.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cmath>

#include <algorithm>

#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/gpu_utils/hostallocator.h"
#include "gromacs/gpu_utils/pinning.h"
#include "gromacs/utility/alignedallocator.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/smalloc.h"

#ifdef NOGMX
#define GMX_PARALLEL_ENV_INITIALIZED 1
#else
#if GMX_MPI
#define GMX_PARALLEL_ENV_INITIALIZED 1
#else
#define GMX_PARALLEL_ENV_INITIALIZED 0
#endif
#endif

#if GMX_OPENMP
/* TODO: Do we still need this? Are we still planning ot use fftw + OpenMP? */
#define FFT5D_THREADS
/* requires fftw compiled with openmp */
/* #define FFT5D_FFTW_THREADS (now set by cmake) */
#endif

#ifndef __FLT_EPSILON__
#define __FLT_EPSILON__ FLT_EPSILON
#define __DBL_EPSILON__ DBL_EPSILON
#endif

#ifdef NOGMX
FILE* debug = 0;
#endif

#if GMX_FFT_FFTW3

#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/mutex.h"
/* none of the fftw3 calls, except execute(), are thread-safe, so
   we need to serialize them with this mutex. */
static gmx::Mutex big_fftw_mutex;
#define FFTW_LOCK try { big_fftw_mutex.lock(); } GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
#define FFTW_UNLOCK try { big_fftw_mutex.unlock(); } GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
#endif /* GMX_FFT_FFTW3 */

#if GMX_MPI
/* largest factor smaller than sqrt */
static int lfactor(int z)
{
    int i = static_cast<int>(sqrt(static_cast<double>(z)));
    while (z%i != 0)
    {
        i--;
    }
    return i;
}
#endif

#if !GMX_MPI
#if HAVE_GETTIMEOFDAY
#include <sys/time.h>
double MPI_Wtime()
{
    struct timeval tv;
    gettimeofday(&tv, 0);
    return tv.tv_sec+tv.tv_usec*1e-6;
}
#else
double MPI_Wtime()
{
    return 0.0;
}
#endif
#endif

static int vmax(int* a, int s)
{
    int i, max = 0;
    for (i = 0; i < s; i++)
    {
        if (a[i] > max)
        {
            max = a[i];
        }
    }
    return max;
}


/* NxMxK the size of the data
 * comm communicator to use for fft5d
 * P0 number of processor in 1st axes (can be null for automatic)
 * lin is allocated by fft5d because size of array is only known after planning phase
 * rlout2 is only used as intermediate buffer - only returned after allocation to reuse for back transform - should not be used by caller
 */
fft5d_plan fft5d_plan_3d(int NG, int MG, int KG, MPI_Comm comm[2], int flags, t_complex** rlin, t_complex** rlout, t_complex** rlout2, t_complex** rlout3, int nthreads, gmx::PinningPolicy realGridAllocationPinningPolicy)
{

    int        P[2], bMaster, prank[2], i, t;
    int        rNG, rMG, rKG;
    int       *N0 = nullptr, *N1 = nullptr, *M0 = nullptr, *M1 = nullptr, *K0 = nullptr, *K1 = nullptr, *oN0 = nullptr, *oN1 = nullptr, *oM0 = nullptr, *oM1 = nullptr, *oK0 = nullptr, *oK1 = nullptr;
    int        N[3], M[3], K[3], pN[3], pM[3], pK[3], oM[3], oK[3], *iNin[3] = {nullptr}, *oNin[3] = {nullptr}, *iNout[3] = {nullptr}, *oNout[3] = {nullptr};
    int        C[3], rC[3], nP[2];
    int        lsize;
    t_complex *lin = nullptr, *lout = nullptr, *lout2 = nullptr, *lout3 = nullptr;
    fft5d_plan plan;
    int        s;

    /* comm, prank and P are in the order of the decomposition (plan->cart is in the order of transposes) */
#if GMX_MPI
    if (GMX_PARALLEL_ENV_INITIALIZED && comm[0] != MPI_COMM_NULL)
    {
        MPI_Comm_size(comm[0], &P[0]);
        MPI_Comm_rank(comm[0], &prank[0]);
    }
    else
#endif
    {
        P[0]     = 1;
        prank[0] = 0;
    }
#if GMX_MPI
    if (GMX_PARALLEL_ENV_INITIALIZED && comm[1] != MPI_COMM_NULL)
    {
        MPI_Comm_size(comm[1], &P[1]);
        MPI_Comm_rank(comm[1], &prank[1]);
    }
    else
#endif
    {
        P[1]     = 1;
        prank[1] = 0;
    }

    bMaster = (prank[0] == 0 && prank[1] == 0);


    if (debug)
    {
        fprintf(debug, "FFT5D: Using %dx%d rank grid, rank %d,%d\n",
                P[0], P[1], prank[0], prank[1]);
    }

    if (bMaster)
    {
        if (debug)
        {
            fprintf(debug, "FFT5D: N: %d, M: %d, K: %d, P: %dx%d, real2complex: %d, backward: %d, order yz: %d, debug %d\n",
                    NG, MG, KG, P[0], P[1], (flags&FFT5D_REALCOMPLEX) > 0, (flags&FFT5D_BACKWARD) > 0, (flags&FFT5D_ORDER_YZ) > 0, (flags&FFT5D_DEBUG) > 0);
        }
        /* The check below is not correct, one prime factor 11 or 13 is ok.
           if (fft5d_fmax(fft5d_fmax(lpfactor(NG),lpfactor(MG)),lpfactor(KG))>7) {
            printf("WARNING: FFT very slow with prime factors larger 7\n");
            printf("Change FFT size or in case you cannot change it look at\n");
            printf("http://www.fftw.org/fftw3_doc/Generating-your-own-code.html\n");
           }
         */
    }

    if (NG == 0 || MG == 0 || KG == 0)
    {
        if (bMaster)
        {
            printf("FFT5D: FATAL: Datasize cannot be zero in any dimension\n");
        }
        return nullptr;
    }

    rNG = NG; rMG = MG; rKG = KG;

    if (flags&FFT5D_REALCOMPLEX)
    {
        if (!(flags&FFT5D_BACKWARD))
        {
            NG = NG/2+1;
        }
        else
        {
            if (!(flags&FFT5D_ORDER_YZ))
            {
                MG = MG/2+1;
            }
            else
            {
                KG = KG/2+1;
            }
        }
    }


    /*for transpose we need to know the size for each processor not only our own size*/

    N0  = (int*)malloc(P[0]*sizeof(int)); N1 = (int*)malloc(P[1]*sizeof(int));
    M0  = (int*)malloc(P[0]*sizeof(int)); M1 = (int*)malloc(P[1]*sizeof(int));
    K0  = (int*)malloc(P[0]*sizeof(int)); K1 = (int*)malloc(P[1]*sizeof(int));
    oN0 = (int*)malloc(P[0]*sizeof(int)); oN1 = (int*)malloc(P[1]*sizeof(int));
    oM0 = (int*)malloc(P[0]*sizeof(int)); oM1 = (int*)malloc(P[1]*sizeof(int));
    oK0 = (int*)malloc(P[0]*sizeof(int)); oK1 = (int*)malloc(P[1]*sizeof(int));

    for (i = 0; i < P[0]; i++)
    {
        #define EVENDIST
        #ifndef EVENDIST
        oN0[i] = i*ceil((double)NG/P[0]);
        oM0[i] = i*ceil((double)MG/P[0]);
        oK0[i] = i*ceil((double)KG/P[0]);
        #else
        oN0[i] = (NG*i)/P[0];
        oM0[i] = (MG*i)/P[0];
        oK0[i] = (KG*i)/P[0];
        #endif
    }
    for (i = 0; i < P[1]; i++)
    {
        #ifndef EVENDIST
        oN1[i] = i*ceil((double)NG/P[1]);
        oM1[i] = i*ceil((double)MG/P[1]);
        oK1[i] = i*ceil((double)KG/P[1]);
        #else
        oN1[i] = (NG*i)/P[1];
        oM1[i] = (MG*i)/P[1];
        oK1[i] = (KG*i)/P[1];
        #endif
    }
    for (i = 0; i < P[0]-1; i++)
    {
        N0[i] = oN0[i+1]-oN0[i];
        M0[i] = oM0[i+1]-oM0[i];
        K0[i] = oK0[i+1]-oK0[i];
    }
    N0[P[0]-1] = NG-oN0[P[0]-1];
    M0[P[0]-1] = MG-oM0[P[0]-1];
    K0[P[0]-1] = KG-oK0[P[0]-1];
    for (i = 0; i < P[1]-1; i++)
    {
        N1[i] = oN1[i+1]-oN1[i];
        M1[i] = oM1[i+1]-oM1[i];
        K1[i] = oK1[i+1]-oK1[i];
    }
    N1[P[1]-1] = NG-oN1[P[1]-1];
    M1[P[1]-1] = MG-oM1[P[1]-1];
    K1[P[1]-1] = KG-oK1[P[1]-1];

    /* for step 1-3 the local N,M,K sizes of the transposed system
       C: contiguous dimension, and nP: number of processor in subcommunicator
       for that step */


    pM[0] = M0[prank[0]];
    oM[0] = oM0[prank[0]];
    pK[0] = K1[prank[1]];
    oK[0] = oK1[prank[1]];
    C[0]  = NG;
    rC[0] = rNG;
    if (!(flags&FFT5D_ORDER_YZ))
    {
        N[0]     = vmax(N1, P[1]);
        M[0]     = M0[prank[0]];
        K[0]     = vmax(K1, P[1]);
        pN[0]    = N1[prank[1]];
        iNout[0] = N1;
        oNout[0] = oN1;
        nP[0]    = P[1];
        C[1]     = KG;
        rC[1]    = rKG;
        N[1]     = vmax(K0, P[0]);
        pN[1]    = K0[prank[0]];
        iNin[1]  = K1;
        oNin[1]  = oK1;
        iNout[1] = K0;
        oNout[1] = oK0;
        M[1]     = vmax(M0, P[0]);
        pM[1]    = M0[prank[0]];
        oM[1]    = oM0[prank[0]];
        K[1]     = N1[prank[1]];
        pK[1]    = N1[prank[1]];
        oK[1]    = oN1[prank[1]];
        nP[1]    = P[0];
        C[2]     = MG;
        rC[2]    = rMG;
        iNin[2]  = M0;
        oNin[2]  = oM0;
        M[2]     = vmax(K0, P[0]);
        pM[2]    = K0[prank[0]];
        oM[2]    = oK0[prank[0]];
        K[2]     = vmax(N1, P[1]);
        pK[2]    = N1[prank[1]];
        oK[2]    = oN1[prank[1]];
        free(N0); free(oN0); /*these are not used for this order*/
        free(M1); free(oM1); /*the rest is freed in destroy*/
    }
    else
    {
        N[0]     = vmax(N0, P[0]);
        M[0]     = vmax(M0, P[0]);
        K[0]     = K1[prank[1]];
        pN[0]    = N0[prank[0]];
        iNout[0] = N0;
        oNout[0] = oN0;
        nP[0]    = P[0];
        C[1]     = MG;
        rC[1]    = rMG;
        N[1]     = vmax(M1, P[1]);
        pN[1]    = M1[prank[1]];
        iNin[1]  = M0;
        oNin[1]  = oM0;
        iNout[1] = M1;
        oNout[1] = oM1;
        M[1]     = N0[prank[0]];
        pM[1]    = N0[prank[0]];
        oM[1]    = oN0[prank[0]];
        K[1]     = vmax(K1, P[1]);
        pK[1]    = K1[prank[1]];
        oK[1]    = oK1[prank[1]];
        nP[1]    = P[1];
        C[2]     = KG;
        rC[2]    = rKG;
        iNin[2]  = K1;
        oNin[2]  = oK1;
        M[2]     = vmax(N0, P[0]);
        pM[2]    = N0[prank[0]];
        oM[2]    = oN0[prank[0]];
        K[2]     = vmax(M1, P[1]);
        pK[2]    = M1[prank[1]];
        oK[2]    = oM1[prank[1]];
        free(N1); free(oN1); /*these are not used for this order*/
        free(K0); free(oK0); /*the rest is freed in destroy*/
    }
    N[2] = pN[2] = -1;       /*not used*/

    /*
       Difference between x-y-z regarding 2d decomposition is whether they are
       distributed along axis 1, 2 or both
     */

    /* int lsize = fmax(N[0]*M[0]*K[0]*nP[0],N[1]*M[1]*K[1]*nP[1]); */
    lsize = std::max(N[0]*M[0]*K[0]*nP[0], std::max(N[1]*M[1]*K[1]*nP[1], C[2]*M[2]*K[2]));
    /* int lsize = fmax(C[0]*M[0]*K[0],fmax(C[1]*M[1]*K[1],C[2]*M[2]*K[2])); */
    if (!(flags&FFT5D_NOMALLOC))
    {
        // only needed for PME GPU mixed mode
        if (realGridAllocationPinningPolicy == gmx::PinningPolicy::CanBePinned)
        {
            const std::size_t numBytes = lsize * sizeof(t_complex);
            lin = static_cast<t_complex *>(gmx::PageAlignedAllocationPolicy::malloc(numBytes));
            gmx::pinBuffer(lin, numBytes);
        }
        else
        {
            snew_aligned(lin, lsize, 32);
        }
        snew_aligned(lout, lsize, 32);
        if (nthreads > 1)
        {
            /* We need extra transpose buffers to avoid OpenMP barriers */
            snew_aligned(lout2, lsize, 32);
            snew_aligned(lout3, lsize, 32);
        }
        else
        {
            /* We can reuse the buffers to avoid cache misses */
            lout2 = lin;
            lout3 = lout;
        }
    }
    else
    {
        lin  = *rlin;
        lout = *rlout;
        if (nthreads > 1)
        {
            lout2 = *rlout2;
            lout3 = *rlout3;
        }
        else
        {
            lout2 = lin;
            lout3 = lout;
        }
    }

    plan = (fft5d_plan)calloc(1, sizeof(struct fft5d_plan_t));


    if (debug)
    {
        fprintf(debug, "Running on %d threads\n", nthreads);
    }

#if GMX_FFT_FFTW3
    /* Don't add more stuff here! We have already had at least one bug because we are reimplementing
     * the low-level FFT interface instead of using the Gromacs FFT module. If we need more
     * generic functionality it is far better to extend the interface so we can use it for
     * all FFT libraries instead of writing FFTW-specific code here.
     */

    /*if not FFTW - then we don't do a 3d plan but instead use only 1D plans */
    /* It is possible to use the 3d plan with OMP threads - but in that case it is not allowed to be called from
     * within a parallel region. For now deactivated. If it should be supported it has to made sure that
     * that the execute of the 3d plan is in a master/serial block (since it contains it own parallel region)
     * and that the 3d plan is faster than the 1d plan.
     */
    if ((!(flags&FFT5D_INPLACE)) && (!(P[0] > 1 || P[1] > 1)) && nthreads == 1) /*don't do 3d plan in parallel or if in_place requested */
    {
        int fftwflags = FFTW_DESTROY_INPUT;
        FFTW(iodim) dims[3];
        int inNG = NG, outMG = MG, outKG = KG;

        FFTW_LOCK;

        fftwflags |= (flags & FFT5D_NOMEASURE) ? FFTW_ESTIMATE : FFTW_MEASURE;

        if (flags&FFT5D_REALCOMPLEX)
        {
            if (!(flags&FFT5D_BACKWARD))        /*input pointer is not complex*/
            {
                inNG *= 2;
            }
            else                                /*output pointer is not complex*/
            {
                if (!(flags&FFT5D_ORDER_YZ))
                {
                    outMG *= 2;
                }
                else
                {
                    outKG *= 2;
                }
            }
        }

        if (!(flags&FFT5D_BACKWARD))
        {
            dims[0].n  = KG;
            dims[1].n  = MG;
            dims[2].n  = rNG;

            dims[0].is = inNG*MG;         /*N M K*/
            dims[1].is = inNG;
            dims[2].is = 1;
            if (!(flags&FFT5D_ORDER_YZ))
            {
                dims[0].os = MG;           /*M K N*/
                dims[1].os = 1;
                dims[2].os = MG*KG;
            }
            else
            {
                dims[0].os = 1;           /*K N M*/
                dims[1].os = KG*NG;
                dims[2].os = KG;
            }
        }
        else
        {
            if (!(flags&FFT5D_ORDER_YZ))
            {
                dims[0].n  = NG;
                dims[1].n  = KG;
                dims[2].n  = rMG;

                dims[0].is = 1;
                dims[1].is = NG*MG;
                dims[2].is = NG;

                dims[0].os = outMG*KG;
                dims[1].os = outMG;
                dims[2].os = 1;
            }
            else
            {
                dims[0].n  = MG;
                dims[1].n  = NG;
                dims[2].n  = rKG;

                dims[0].is = NG;
                dims[1].is = 1;
                dims[2].is = NG*MG;

                dims[0].os = outKG*NG;
                dims[1].os = outKG;
                dims[2].os = 1;
            }
        }
#ifdef FFT5D_THREADS
#ifdef FFT5D_FFTW_THREADS
        FFTW(plan_with_nthreads)(nthreads);
#endif
#endif
        if ((flags&FFT5D_REALCOMPLEX) && !(flags&FFT5D_BACKWARD))
        {
            plan->p3d = FFTW(plan_guru_dft_r2c)(/*rank*/ 3, dims,
                                                         /*howmany*/ 0, /*howmany_dims*/ nullptr,
                                                         (real*)lin, (FFTW(complex) *) lout,
                                                         /*flags*/ fftwflags);
        }
        else if ((flags&FFT5D_REALCOMPLEX) && (flags&FFT5D_BACKWARD))
        {
            plan->p3d = FFTW(plan_guru_dft_c2r)(/*rank*/ 3, dims,
                                                         /*howmany*/ 0, /*howmany_dims*/ nullptr,
                                                         (FFTW(complex) *) lin, (real*)lout,
                                                         /*flags*/ fftwflags);
        }
        else
        {
            plan->p3d = FFTW(plan_guru_dft)(/*rank*/ 3, dims,
                                                     /*howmany*/ 0, /*howmany_dims*/ nullptr,
                                                     (FFTW(complex) *) lin, (FFTW(complex) *) lout,
                                                     /*sign*/ (flags&FFT5D_BACKWARD) ? 1 : -1, /*flags*/ fftwflags);
        }
#ifdef FFT5D_THREADS
#ifdef FFT5D_FFTW_THREADS
        FFTW(plan_with_nthreads)(1);
#endif
#endif
        FFTW_UNLOCK;
    }
    if (!plan->p3d) /* for decomposition and if 3d plan did not work */
    {
#endif              /* GMX_FFT_FFTW3 */
    for (s = 0; s < 3; s++)
    {
        if (debug)
        {
            fprintf(debug, "FFT5D: Plan s %d rC %d M %d pK %d C %d lsize %d\n",
                    s, rC[s], M[s], pK[s], C[s], lsize);
        }
        plan->p1d[s] = (gmx_fft_t*)malloc(sizeof(gmx_fft_t)*nthreads);

        /* Make sure that the init routines are only called by one thread at a time and in order
           (later is only important to not confuse valgrind)
         */
#pragma omp parallel for num_threads(nthreads) schedule(static) ordered
        for (t = 0; t < nthreads; t++)
        {
#pragma omp ordered
            {
                try
                {
                    int tsize = ((t+1)*pM[s]*pK[s]/nthreads)-(t*pM[s]*pK[s]/nthreads);

                    if ((flags&FFT5D_REALCOMPLEX) && ((!(flags&FFT5D_BACKWARD) && s == 0) || ((flags&FFT5D_BACKWARD) && s == 2)))
                    {
                        gmx_fft_init_many_1d_real( &plan->p1d[s][t], rC[s], tsize, (flags&FFT5D_NOMEASURE) ? GMX_FFT_FLAG_CONSERVATIVE : 0 );
                    }
                    else
                    {
                        gmx_fft_init_many_1d     ( &plan->p1d[s][t],  C[s], tsize, (flags&FFT5D_NOMEASURE) ? GMX_FFT_FLAG_CONSERVATIVE : 0 );
                    }
                }
                GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
            }
        }
    }

#if GMX_FFT_FFTW3
}
#endif
    if ((flags&FFT5D_ORDER_YZ))   /*plan->cart is in the order of transposes */
    {
        plan->cart[0] = comm[0]; plan->cart[1] = comm[1];
    }
    else
    {
        plan->cart[1] = comm[0]; plan->cart[0] = comm[1];
    }
#ifdef FFT5D_MPI_TRANSPOSE
    FFTW_LOCK;
    for (s = 0; s < 2; s++)
    {
        if ((s == 0 && !(flags&FFT5D_ORDER_YZ)) || (s == 1 && (flags&FFT5D_ORDER_YZ)))
        {
            plan->mpip[s] = FFTW(mpi_plan_many_transpose)(nP[s], nP[s], N[s]*K[s]*pM[s]*2, 1, 1, (real*)lout2, (real*)lout3, plan->cart[s], FFTW_PATIENT);
        }
        else
        {
            plan->mpip[s] = FFTW(mpi_plan_many_transpose)(nP[s], nP[s], N[s]*pK[s]*M[s]*2, 1, 1, (real*)lout2, (real*)lout3, plan->cart[s], FFTW_PATIENT);
        }
    }
    FFTW_UNLOCK;
#endif


    plan->lin   = lin;
    plan->lout  = lout;
    plan->lout2 = lout2;
    plan->lout3 = lout3;

    plan->NG = NG; plan->MG = MG; plan->KG = KG;

    for (s = 0; s < 3; s++)
    {
        plan->N[s]    = N[s]; plan->M[s] = M[s]; plan->K[s] = K[s]; plan->pN[s] = pN[s]; plan->pM[s] = pM[s]; plan->pK[s] = pK[s];
        plan->oM[s]   = oM[s]; plan->oK[s] = oK[s];
        plan->C[s]    = C[s]; plan->rC[s] = rC[s];
        plan->iNin[s] = iNin[s]; plan->oNin[s] = oNin[s]; plan->iNout[s] = iNout[s]; plan->oNout[s] = oNout[s];
    }
    for (s = 0; s < 2; s++)
    {
        plan->P[s] = nP[s]; plan->coor[s] = prank[s];
    }

/*    plan->fftorder=fftorder;
    plan->direction=direction;
    plan->realcomplex=realcomplex;
 */
    plan->flags         = flags;
    plan->nthreads      = nthreads;
    plan->pinningPolicy = realGridAllocationPinningPolicy;
    *rlin               = lin;
    *rlout              = lout;
    *rlout2             = lout2;
    *rlout3             = lout3;
    return plan;
}


enum order {
    XYZ,
    XZY,
    YXZ,
    YZX,
    ZXY,
    ZYX
};



/*here x,y,z and N,M,K is in rotated coordinate system!!
   x (and N) is mayor (consecutive) dimension, y (M) middle and z (K) major
   maxN,maxM,maxK is max size of local data
   pN, pM, pK is local size specific to current processor (only different to max if not divisible)
   NG, MG, KG is size of global data*/
static void splitaxes(t_complex* lout, const t_complex* lin,
                      int maxN, int maxM, int maxK, int pM,
                      int P, int NG, int *N, int* oN, int starty, int startz, int endy, int endz)
{
    int x, y, z, i;
    int in_i, out_i, in_z, out_z, in_y, out_y;
    int s_y, e_y;

    for (z = startz; z < endz+1; z++) /*3. z l*/
    {
        if (z == startz)
        {
            s_y = starty;
        }
        else
        {
            s_y = 0;
        }
        if (z == endz)
        {
            e_y = endy;
        }
        else
        {
            e_y = pM;
        }
        out_z  = z*maxN*maxM;
        in_z   = z*NG*pM;

        for (i = 0; i < P; i++) /*index cube along long axis*/
        {
            out_i  = out_z  + i*maxN*maxM*maxK;
            in_i   = in_z + oN[i];
            for (y = s_y; y < e_y; y++)   /*2. y k*/
            {
                out_y  = out_i  + y*maxN;
                in_y   = in_i + y*NG;
                for (x = 0; x < N[i]; x++)       /*1. x j*/
                {
                    lout[out_y+x] = lin[in_y+x]; /*in=z*NG*pM+oN[i]+y*NG+x*/
                    /*after split important that each processor chunk i has size maxN*maxM*maxK and thus being the same size*/
                    /*before split data contiguos - thus if different processor get different amount oN is different*/
                }
            }
        }
    }
}

/*make axis contiguous again (after AllToAll) and also do local transpose*/
/*transpose mayor and major dimension
   variables see above
   the major, middle, minor order is only correct for x,y,z (N,M,K) for the input
   N,M,K local dimensions
   KG global size*/
static void joinAxesTrans13(t_complex* lout, const t_complex* lin,
                            int maxN, int maxM, int maxK, int pM,
                            int P, int KG, int* K, int* oK, int starty, int startx, int endy, int endx)
{
    int i, x, y, z;
    int out_i, in_i, out_x, in_x, out_z, in_z;
    int s_y, e_y;

    for (x = startx; x < endx+1; x++) /*1.j*/
    {
        if (x == startx)
        {
            s_y = starty;
        }
        else
        {
            s_y = 0;
        }
        if (x == endx)
        {
            e_y = endy;
        }
        else
        {
            e_y = pM;
        }

        out_x  = x*KG*pM;
        in_x   = x;

        for (i = 0; i < P; i++) /*index cube along long axis*/
        {
            out_i  = out_x  + oK[i];
            in_i   = in_x + i*maxM*maxN*maxK;
            for (z = 0; z < K[i]; z++) /*3.l*/
            {
                out_z  = out_i  + z;
                in_z   = in_i + z*maxM*maxN;
                for (y = s_y; y < e_y; y++)              /*2.k*/
                {
                    lout[out_z+y*KG] = lin[in_z+y*maxN]; /*out=x*KG*pM+oK[i]+z+y*KG*/
                }
            }
        }
    }
}

/*make axis contiguous again (after AllToAll) and also do local transpose
   tranpose mayor and middle dimension
   variables see above
   the minor, middle, major order is only correct for x,y,z (N,M,K) for the input
   N,M,K local size
   MG, global size*/
static void joinAxesTrans12(t_complex* lout, const t_complex* lin, int maxN, int maxM, int maxK, int pN,
                            int P, int MG, int* M, int* oM, int startx, int startz, int endx, int endz)
{
    int i, z, y, x;
    int out_i, in_i, out_z, in_z, out_x, in_x;
    int s_x, e_x;

    for (z = startz; z < endz+1; z++)
    {
        if (z == startz)
        {
            s_x = startx;
        }
        else
        {
            s_x = 0;
        }
        if (z == endz)
        {
            e_x = endx;
        }
        else
        {
            e_x = pN;
        }
        out_z  = z*MG*pN;
        in_z   = z*maxM*maxN;

        for (i = 0; i < P; i++) /*index cube along long axis*/
        {
            out_i  = out_z  + oM[i];
            in_i   = in_z + i*maxM*maxN*maxK;
            for (x = s_x; x < e_x; x++)
            {
                out_x  = out_i  + x*MG;
                in_x   = in_i + x;
                for (y = 0; y < M[i]; y++)
                {
                    lout[out_x+y] = lin[in_x+y*maxN]; /*out=z*MG*pN+oM[i]+x*MG+y*/
                }
            }
        }
    }
}


static void rotate_offsets(int x[])
{
    int t = x[0];
/*    x[0]=x[2];
    x[2]=x[1];
    x[1]=t;*/
    x[0] = x[1];
    x[1] = x[2];
    x[2] = t;
}

/*compute the offset to compare or print transposed local data in original input coordinates
   xs matrix dimension size, xl dimension length, xc decomposition offset
   s: step in computation = number of transposes*/
static void compute_offsets(fft5d_plan plan, int xs[], int xl[], int xc[], int NG[], int s)
{
/*    int direction = plan->direction;
    int fftorder = plan->fftorder;*/

    int  o = 0;
    int  pos[3], i;
    int *pM = plan->pM, *pK = plan->pK, *oM = plan->oM, *oK = plan->oK,
    *C      = plan->C, *rC = plan->rC;

    NG[0] = plan->NG; NG[1] = plan->MG; NG[2] = plan->KG;

    if (!(plan->flags&FFT5D_ORDER_YZ))
    {
        switch (s)
        {
            case 0: o = XYZ; break;
            case 1: o = ZYX; break;
            case 2: o = YZX; break;
            default: assert(0);
        }
    }
    else
    {
        switch (s)
        {
            case 0: o = XYZ; break;
            case 1: o = YXZ; break;
            case 2: o = ZXY; break;
            default: assert(0);
        }
    }

    switch (o)
    {
        case XYZ: pos[0] = 1; pos[1] = 2; pos[2] = 3; break;
        case XZY: pos[0] = 1; pos[1] = 3; pos[2] = 2; break;
        case YXZ: pos[0] = 2; pos[1] = 1; pos[2] = 3; break;
        case YZX: pos[0] = 3; pos[1] = 1; pos[2] = 2; break;
        case ZXY: pos[0] = 2; pos[1] = 3; pos[2] = 1; break;
        case ZYX: pos[0] = 3; pos[1] = 2; pos[2] = 1; break;
    }
    /*if (debug) printf("pos: %d %d %d\n",pos[0],pos[1],pos[2]);*/

    /*xs, xl give dimension size and data length in local transposed coordinate system
       for 0(/1/2): x(/y/z) in original coordinate system*/
    for (i = 0; i < 3; i++)
    {
        switch (pos[i])
        {
            case 1: xs[i] = 1;         xc[i] = 0;     xl[i] = C[s]; break;
            case 2: xs[i] = C[s];      xc[i] = oM[s]; xl[i] = pM[s]; break;
            case 3: xs[i] = C[s]*pM[s]; xc[i] = oK[s]; xl[i] = pK[s]; break;
        }
    }
    /*input order is different for test program to match FFTW order
       (important for complex to real)*/
    if (plan->flags&FFT5D_BACKWARD)
    {
        rotate_offsets(xs);
        rotate_offsets(xl);
        rotate_offsets(xc);
        rotate_offsets(NG);
        if (plan->flags&FFT5D_ORDER_YZ)
        {
            rotate_offsets(xs);
            rotate_offsets(xl);
            rotate_offsets(xc);
            rotate_offsets(NG);
        }
    }
    if ((plan->flags&FFT5D_REALCOMPLEX) && ((!(plan->flags&FFT5D_BACKWARD) && s == 0) || ((plan->flags&FFT5D_BACKWARD) && s == 2)))
    {
        xl[0] = rC[s];
    }
}

static void print_localdata(const t_complex* lin, const char* txt, int s, fft5d_plan plan)
{
    int  x, y, z, l;
    int *coor = plan->coor;
    int  xs[3], xl[3], xc[3], NG[3];
    int  ll = (plan->flags&FFT5D_REALCOMPLEX) ? 1 : 2;
    compute_offsets(plan, xs, xl, xc, NG, s);
    fprintf(debug, txt, coor[0], coor[1]);
    /*printf("xs: %d %d %d, xl: %d %d %d\n",xs[0],xs[1],xs[2],xl[0],xl[1],xl[2]);*/
    for (z = 0; z < xl[2]; z++)
    {
        for (y = 0; y < xl[1]; y++)
        {
            fprintf(debug, "%d %d: ", coor[0], coor[1]);
            for (x = 0; x < xl[0]; x++)
            {
                for (l = 0; l < ll; l++)
                {
                    fprintf(debug, "%f ", ((real*)lin)[(z*xs[2]+y*xs[1])*2+(x*xs[0])*ll+l]);
                }
                fprintf(debug, ",");
            }
            fprintf(debug, "\n");
        }
    }
}

void fft5d_execute(fft5d_plan plan, int thread, fft5d_time times)
{
    t_complex  *lin   = plan->lin;
    t_complex  *lout  = plan->lout;
    t_complex  *lout2 = plan->lout2;
    t_complex  *lout3 = plan->lout3;
    t_complex  *fftout, *joinin;

    gmx_fft_t **p1d = plan->p1d;
#ifdef FFT5D_MPI_TRANSPOSE
    FFTW(plan) *mpip = plan->mpip;
#endif
#if GMX_MPI
    MPI_Comm *cart = plan->cart;
#endif
#ifdef NOGMX
    double time_fft = 0, time_local = 0, time_mpi[2] = {0}, time = 0;
#endif
    int   *N = plan->N, *M = plan->M, *K = plan->K, *pN = plan->pN, *pM = plan->pM, *pK = plan->pK,
    *C       = plan->C, *P = plan->P, **iNin = plan->iNin, **oNin = plan->oNin, **iNout = plan->iNout, **oNout = plan->oNout;
    int    s = 0, tstart, tend, bParallelDim;


#if GMX_FFT_FFTW3
    if (plan->p3d)
    {
        if (thread == 0)
        {
#ifdef NOGMX
            if (times != 0)
            {
                time = MPI_Wtime();
            }
#endif
            FFTW(execute)(plan->p3d);
#ifdef NOGMX
            if (times != 0)
            {
                times->fft += MPI_Wtime()-time;
            }
#endif
        }
        return;
    }
#endif

    s = 0;

    /*lin: x,y,z*/
    if ((plan->flags&FFT5D_DEBUG) && thread == 0)
    {
        print_localdata(lin, "%d %d: copy in lin\n", s, plan);
    }

    for (s = 0; s < 2; s++)  /*loop over first two FFT steps (corner rotations)*/

    {
#if GMX_MPI
        if (GMX_PARALLEL_ENV_INITIALIZED && cart[s] != MPI_COMM_NULL && P[s] > 1)
        {
            bParallelDim = 1;
        }
        else
#endif
        {
            bParallelDim = 0;
        }

        /* ---------- START FFT ------------ */
#ifdef NOGMX
        if (times != 0 && thread == 0)
        {
            time = MPI_Wtime();
        }
#endif

        if (bParallelDim || plan->nthreads == 1)
        {
            fftout = lout;
        }
        else
        {
            if (s == 0)
            {
                fftout = lout3;
            }
            else
            {
                fftout = lout2;
            }
        }

        tstart = (thread*pM[s]*pK[s]/plan->nthreads)*C[s];
        if ((plan->flags&FFT5D_REALCOMPLEX) && !(plan->flags&FFT5D_BACKWARD) && s == 0)
        {
            gmx_fft_many_1d_real(p1d[s][thread], (plan->flags&FFT5D_BACKWARD) ? GMX_FFT_COMPLEX_TO_REAL : GMX_FFT_REAL_TO_COMPLEX, lin+tstart, fftout+tstart);
        }
        else
        {
            gmx_fft_many_1d(     p1d[s][thread], (plan->flags&FFT5D_BACKWARD) ? GMX_FFT_BACKWARD : GMX_FFT_FORWARD,               lin+tstart, fftout+tstart);

        }

#ifdef NOGMX
        if (times != NULL && thread == 0)
        {
            time_fft += MPI_Wtime()-time;
        }
#endif
        if ((plan->flags&FFT5D_DEBUG) && thread == 0)
        {
            print_localdata(lout, "%d %d: FFT %d\n", s, plan);
        }
        /* ---------- END FFT ------------ */

        /* ---------- START SPLIT + TRANSPOSE------------ (if parallel in in this dimension)*/
        if (bParallelDim)
        {
#ifdef NOGMX
            if (times != NULL && thread == 0)
            {
                time = MPI_Wtime();
            }
#endif
            /*prepare for A
               llToAll
               1. (most outer) axes (x) is split into P[s] parts of size N[s]
               for sending*/
            if (pM[s] > 0)
            {
                tend    = ((thread+1)*pM[s]*pK[s]/plan->nthreads);
                tstart /= C[s];
                splitaxes(lout2, lout, N[s], M[s], K[s], pM[s], P[s], C[s], iNout[s], oNout[s], tstart%pM[s], tstart/pM[s], tend%pM[s], tend/pM[s]);
            }
#pragma omp barrier /*barrier required before AllToAll (all input has to be their) - before timing to make timing more acurate*/
#ifdef NOGMX
            if (times != NULL && thread == 0)
            {
                time_local += MPI_Wtime()-time;
            }
#endif

            /* ---------- END SPLIT , START TRANSPOSE------------ */

            if (thread == 0)
            {
#ifdef NOGMX
                if (times != 0)
                {
                    time = MPI_Wtime();
                }
#else
                wallcycle_start(times, ewcPME_FFTCOMM);
#endif
#ifdef FFT5D_MPI_TRANSPOSE
                FFTW(execute)(mpip[s]);
#else
#if GMX_MPI
                if ((s == 0 && !(plan->flags&FFT5D_ORDER_YZ)) || (s == 1 && (plan->flags&FFT5D_ORDER_YZ)))
                {
                    MPI_Alltoall((real *)lout2, N[s]*pM[s]*K[s]*sizeof(t_complex)/sizeof(real), GMX_MPI_REAL, (real *)lout3, N[s]*pM[s]*K[s]*sizeof(t_complex)/sizeof(real), GMX_MPI_REAL, cart[s]);
                }
                else
                {
                    MPI_Alltoall((real *)lout2, N[s]*M[s]*pK[s]*sizeof(t_complex)/sizeof(real), GMX_MPI_REAL, (real *)lout3, N[s]*M[s]*pK[s]*sizeof(t_complex)/sizeof(real), GMX_MPI_REAL, cart[s]);
                }
#else
                gmx_incons("fft5d MPI call without MPI configuration");
#endif /*GMX_MPI*/
#endif /*FFT5D_MPI_TRANSPOSE*/
#ifdef NOGMX
                if (times != 0)
                {
                    time_mpi[s] = MPI_Wtime()-time;
                }
#else
                wallcycle_stop(times, ewcPME_FFTCOMM);
#endif
            } /*master*/
        }     /* bPrallelDim */
#pragma omp barrier  /*both needed for parallel and non-parallel dimension (either have to wait on data from AlltoAll or from last FFT*/

        /* ---------- END SPLIT + TRANSPOSE------------ */

        /* ---------- START JOIN ------------ */
#ifdef NOGMX
        if (times != NULL && thread == 0)
        {
            time = MPI_Wtime();
        }
#endif

        if (bParallelDim)
        {
            joinin = lout3;
        }
        else
        {
            joinin = fftout;
        }
        /*bring back in matrix form
           thus make  new 1. axes contiguos
           also local transpose 1 and 2/3
           runs on thread used for following FFT (thus needing a barrier before but not afterwards)
         */
        if ((s == 0 && !(plan->flags&FFT5D_ORDER_YZ)) || (s == 1 && (plan->flags&FFT5D_ORDER_YZ)))
        {
            if (pM[s] > 0)
            {
                tstart = ( thread   *pM[s]*pN[s]/plan->nthreads);
                tend   = ((thread+1)*pM[s]*pN[s]/plan->nthreads);
                joinAxesTrans13(lin, joinin, N[s], pM[s], K[s], pM[s], P[s], C[s+1], iNin[s+1], oNin[s+1], tstart%pM[s], tstart/pM[s], tend%pM[s], tend/pM[s]);
            }
        }
        else
        {
            if (pN[s] > 0)
            {
                tstart = ( thread   *pK[s]*pN[s]/plan->nthreads);
                tend   = ((thread+1)*pK[s]*pN[s]/plan->nthreads);
                joinAxesTrans12(lin, joinin, N[s], M[s], pK[s], pN[s], P[s], C[s+1], iNin[s+1], oNin[s+1], tstart%pN[s], tstart/pN[s], tend%pN[s], tend/pN[s]);
            }
        }

#ifdef NOGMX
        if (times != NULL && thread == 0)
        {
            time_local += MPI_Wtime()-time;
        }
#endif
        if ((plan->flags&FFT5D_DEBUG) && thread == 0)
        {
            print_localdata(lin, "%d %d: tranposed %d\n", s+1, plan);
        }
        /* ---------- END JOIN ------------ */

        /*if (debug) print_localdata(lin, "%d %d: transposed x-z\n", N1, M0, K, ZYX, coor);*/
    }  /* for(s=0;s<2;s++) */
#ifdef NOGMX
    if (times != NULL && thread == 0)
    {
        time = MPI_Wtime();
    }
#endif

    if (plan->flags&FFT5D_INPLACE)
    {
        lout = lin;                          /*in place currently not supported*/

    }
    /*  ----------- FFT ----------- */
    tstart = (thread*pM[s]*pK[s]/plan->nthreads)*C[s];
    if ((plan->flags&FFT5D_REALCOMPLEX) && (plan->flags&FFT5D_BACKWARD))
    {
        gmx_fft_many_1d_real(p1d[s][thread], (plan->flags&FFT5D_BACKWARD) ? GMX_FFT_COMPLEX_TO_REAL : GMX_FFT_REAL_TO_COMPLEX, lin+tstart, lout+tstart);
    }
    else
    {
        gmx_fft_many_1d(     p1d[s][thread], (plan->flags&FFT5D_BACKWARD) ? GMX_FFT_BACKWARD : GMX_FFT_FORWARD,               lin+tstart, lout+tstart);
    }
    /* ------------ END FFT ---------*/

#ifdef NOGMX
    if (times != NULL && thread == 0)
    {
        time_fft += MPI_Wtime()-time;

        times->fft   += time_fft;
        times->local += time_local;
        times->mpi2  += time_mpi[1];
        times->mpi1  += time_mpi[0];
    }
#endif

    if ((plan->flags&FFT5D_DEBUG) && thread == 0)
    {
        print_localdata(lout, "%d %d: FFT %d\n", s, plan);
    }
}

void fft5d_destroy(fft5d_plan plan)
{
    int s, t;

    for (s = 0; s < 3; s++)
    {
        if (plan->p1d[s])
        {
            for (t = 0; t < plan->nthreads; t++)
            {
                gmx_many_fft_destroy(plan->p1d[s][t]);
            }
            free(plan->p1d[s]);
        }
        if (plan->iNin[s])
        {
            free(plan->iNin[s]);
            plan->iNin[s] = nullptr;
        }
        if (plan->oNin[s])
        {
            free(plan->oNin[s]);
            plan->oNin[s] = nullptr;
        }
        if (plan->iNout[s])
        {
            free(plan->iNout[s]);
            plan->iNout[s] = nullptr;
        }
        if (plan->oNout[s])
        {
            free(plan->oNout[s]);
            plan->oNout[s] = nullptr;
        }
    }
#if GMX_FFT_FFTW3
    FFTW_LOCK;
#ifdef FFT5D_MPI_TRANSPOS
    for (s = 0; s < 2; s++)
    {
        FFTW(destroy_plan)(plan->mpip[s]);
    }
#endif /* FFT5D_MPI_TRANSPOS */
    if (plan->p3d)
    {
        FFTW(destroy_plan)(plan->p3d);
    }
    FFTW_UNLOCK;
#endif /* GMX_FFT_FFTW3 */

    if (!(plan->flags&FFT5D_NOMALLOC))
    {
        // only needed for PME GPU mixed mode
        if (plan->pinningPolicy == gmx::PinningPolicy::CanBePinned &&
            isHostMemoryPinned(plan->lin))
        {
            gmx::unpinBuffer(plan->lin);
        }
        sfree_aligned(plan->lin);
        sfree_aligned(plan->lout);
        if (plan->nthreads > 1)
        {
            sfree_aligned(plan->lout2);
            sfree_aligned(plan->lout3);
        }
    }

#ifdef FFT5D_THREADS
#ifdef FFT5D_FFTW_THREADS
    /*FFTW(cleanup_threads)();*/
#endif
#endif

    free(plan);
}

/*Is this better than direct access of plan? enough data?
   here 0,1 reference divided by which processor grid dimension (not FFT step!)*/
void fft5d_local_size(fft5d_plan plan, int* N1, int* M0, int* K0, int* K1, int** coor)
{
    *N1 = plan->N[0];
    *M0 = plan->M[0];
    *K1 = plan->K[0];
    *K0 = plan->N[1];

    *coor = plan->coor;
}


/*same as fft5d_plan_3d but with cartesian coordinator and automatic splitting
   of processor dimensions*/
fft5d_plan fft5d_plan_3d_cart(int NG, int MG, int KG, MPI_Comm comm, int P0, int flags, t_complex** rlin, t_complex** rlout, t_complex** rlout2, t_complex** rlout3, int nthreads)
{
    MPI_Comm cart[2] = {MPI_COMM_NULL, MPI_COMM_NULL};
#if GMX_MPI
    int      size = 1, prank = 0;
    int      P[2];
    int      coor[2];
    int      wrap[] = {0, 0};
    MPI_Comm gcart;
    int      rdim1[] = {0, 1}, rdim2[] = {1, 0};

    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &prank);

    if (P0 == 0)
    {
        P0 = lfactor(size);
    }
    if (size%P0 != 0)
    {
        if (prank == 0)
        {
            printf("FFT5D: WARNING: Number of ranks %d not evenly divisible by %d\n", size, P0);
        }
        P0 = lfactor(size);
    }

    P[0] = P0; P[1] = size/P0; /*number of processors in the two dimensions*/

    /*Difference between x-y-z regarding 2d decomposition is whether they are
       distributed along axis 1, 2 or both*/

    MPI_Cart_create(comm, 2, P, wrap, 1, &gcart); /*parameter 4: value 1: reorder*/
    MPI_Cart_get(gcart, 2, P, wrap, coor);
    MPI_Cart_sub(gcart, rdim1, &cart[0]);
    MPI_Cart_sub(gcart, rdim2, &cart[1]);
#else
    (void)P0;
    (void)comm;
#endif
    return fft5d_plan_3d(NG, MG, KG, cart, flags, rlin, rlout, rlout2, rlout3, nthreads);
}



/*prints in original coordinate system of data (as the input to FFT)*/
void fft5d_compare_data(const t_complex* lin, const t_complex* in, fft5d_plan plan, int bothLocal, int normalize)
{
    int  xs[3], xl[3], xc[3], NG[3];
    int  x, y, z, l;
    int *coor = plan->coor;
    int  ll   = 2; /*compare ll values per element (has to be 2 for complex)*/
    if ((plan->flags&FFT5D_REALCOMPLEX) && (plan->flags&FFT5D_BACKWARD))
    {
        ll = 1;
    }

    compute_offsets(plan, xs, xl, xc, NG, 2);
    if (plan->flags&FFT5D_DEBUG)
    {
        printf("Compare2\n");
    }
    for (z = 0; z < xl[2]; z++)
    {
        for (y = 0; y < xl[1]; y++)
        {
            if (plan->flags&FFT5D_DEBUG)
            {
                printf("%d %d: ", coor[0], coor[1]);
            }
            for (x = 0; x < xl[0]; x++)
            {
                for (l = 0; l < ll; l++)   /*loop over real/complex parts*/
                {
                    real a, b;
                    a = ((real*)lin)[(z*xs[2]+y*xs[1])*2+x*xs[0]*ll+l];
                    if (normalize)
                    {
                        a /= plan->rC[0]*plan->rC[1]*plan->rC[2];
                    }
                    if (!bothLocal)
                    {
                        b = ((real*)in)[((z+xc[2])*NG[0]*NG[1]+(y+xc[1])*NG[0])*2+(x+xc[0])*ll+l];
                    }
                    else
                    {
                        b = ((real*)in)[(z*xs[2]+y*xs[1])*2+x*xs[0]*ll+l];
                    }
                    if (plan->flags&FFT5D_DEBUG)
                    {
                        printf("%f %f, ", a, b);
                    }
                    else
                    {
                        if (fabs(a-b) > 2*NG[0]*NG[1]*NG[2]*GMX_REAL_EPS)
                        {
                            printf("result incorrect on %d,%d at %d,%d,%d: FFT5D:%f reference:%f\n", coor[0], coor[1], x, y, z, a, b);
                        }
/*                        assert(fabs(a-b)<2*NG[0]*NG[1]*NG[2]*GMX_REAL_EPS);*/
                    }
                }
                if (plan->flags&FFT5D_DEBUG)
                {
                    printf(",");
                }
            }
            if (plan->flags&FFT5D_DEBUG)
            {
                printf("\n");
            }
        }
    }

}
