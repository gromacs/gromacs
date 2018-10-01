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

#ifndef GMX_FFT_FFT5D_H
#define GMX_FFT_FFT5D_H

#include "config.h"

#ifdef NOGMX
/*#define GMX_MPI*/
/*#define GMX_FFT_FFTW3*/
FILE* debug;
#endif

#include "gromacs/fft/fft.h"
#include "gromacs/gpu_utils/hostallocator.h"
#include "gromacs/math/gmxcomplex.h"
#include "gromacs/utility/gmxmpi.h"

#if !GMX_MPI
double MPI_Wtime();
#endif

/*currently only special optimization for FFTE*/
#if GMX_FFT_FFTW3
#include <fftw3.h>
#endif

#if !GMX_DOUBLE
#define FFTW(x) fftwf_ ## x
#else
#define FFTW(x) fftw_ ## x
#endif

#ifdef NOGMX
struct fft5d_time_t {
    double fft, local, mpi1, mpi2;
};
typedef struct fft5d_time_t *fft5d_time;
#else
#include "gromacs/timing/wallcycle.h"
typedef gmx_wallcycle_t fft5d_time;
#endif

namespace gmx
{
enum class PinningPolicy : int;
} // namespace

typedef enum fft5d_flags_t {
    FFT5D_ORDER_YZ    = 1,
    FFT5D_BACKWARD    = 2,
    FFT5D_REALCOMPLEX = 4,
    FFT5D_DEBUG       = 8,
    FFT5D_NOMEASURE   = 16,
    FFT5D_INPLACE     = 32,
    FFT5D_NOMALLOC    = 64
} fft5d_flags;

struct fft5d_plan_t {
    t_complex *lin;
    t_complex *lout, *lout2, *lout3;
    gmx_fft_t* p1d[3]; /*1D plans*/
#if GMX_FFT_FFTW3
    FFTW(plan) p2d;    /*2D plan: used for 1D decomposition if FFT supports transposed output*/
    FFTW(plan) p3d;    /*3D plan: used for 0D decomposition if FFT supports transposed output*/
    FFTW(plan) mpip[2];
#endif
    MPI_Comm cart[2];

    int      N[3], M[3], K[3];                        /*local length in transposed coordinate system (if not divisisable max)*/
    int      pN[3], pM[3], pK[3];                     /*local length - not max but length for this processor*/
    int      oM[3], oK[3];                            /*offset for current processor*/
    int     *iNin[3], *oNin[3], *iNout[3], *oNout[3]; /*size for each processor (if divisisable=max) for out(=split)
                                                         and in (=join) and offsets in transposed coordinate system*/
    int      C[3], rC[3];                             /*global length (of the one global axes) */
    /* C!=rC for real<->complex. then C=rC/2 but with potential padding*/
    int      P[2];                                    /*size of processor grid*/
/*  int fftorder;*/
/*  int direction;*/
/*  int realcomplex;*/
    int                flags;
    /*int N0,N1,M0,M1,K0,K1;*/
    int                NG, MG, KG;
    /*int P[2];*/
    int                coor[2];
    int                nthreads;
    gmx::PinningPolicy pinningPolicy;
};

typedef struct fft5d_plan_t *fft5d_plan;

void fft5d_execute(fft5d_plan plan, int thread, fft5d_time times);
fft5d_plan fft5d_plan_3d(int N, int M, int K, MPI_Comm comm[2], int flags, t_complex**lin, t_complex**lin2, t_complex**lout2, t_complex**lout3, int nthreads, gmx::PinningPolicy realGridAllocationPinningPolicy = gmx::PinningPolicy::CannotBePinned);
void fft5d_local_size(fft5d_plan plan, int* N1, int* M0, int* K0, int* K1, int** coor);
void fft5d_destroy(fft5d_plan plan);
fft5d_plan fft5d_plan_3d_cart(int N, int M, int K, MPI_Comm comm, int P0, int flags, t_complex** lin, t_complex** lin2, t_complex** lout2, t_complex** lout3, int nthreads);
void fft5d_compare_data(const t_complex* lin, const t_complex* in, fft5d_plan plan, int bothLocal, int normarlize);

#endif
