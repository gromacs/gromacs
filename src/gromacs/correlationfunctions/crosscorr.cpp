/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2014- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
/*! \internal
 * \file
 * \brief
 * Implements routine for computing a cross correlation between two data sets
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_correlationfunctions
 */
#include "gmxpre.h"

#include "crosscorr.h"

#include "gromacs/fft/fft.h"
#include "gromacs/math/gmxcomplex.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/smalloc.h"

/*! \brief
 * return the size witch the array should be after zero padding
 *
 * \param[in] n Factor to multiply by 2
 * \return zeroPaddingSize
 */
static int zeroPaddingSize(int n)
{
    return 2 * n;
}

/*! \brief
 * Compute complex conjugate. Output in the first input variable.
 *
 * \param[in] in1 first complex number
 * \param[in] in2 second complex number
 */
static void complexConjugatMult(t_complex* in1, t_complex* in2)
{
    t_complex res;
    res.re  = in1->re * in2->re + in1->im * in2->im;
    res.im  = in1->re * -in2->im + in1->im * in2->re;
    in1->re = res.re;
    in1->im = res.im;
}

/*! \brief
 * Compute one cross correlation corr = f x g using FFT.
 *
 * \param[in] n number of data point
 * \param[in] f first function
 * \param[in] g second function
 * \param[out] corr output correlation
 * \param[in] fft FFT data structure
 */
static void cross_corr_low(int n, const real f[], const real g[], real corr[], gmx_fft_t fft)
{
    int        i;
    const int  size = zeroPaddingSize(n);
    t_complex *in1, *in2;

    snew(in1, size);
    snew(in2, size);

    for (i = 0; i < n; i++)
    {
        in1[i].re = f[i];
        in1[i].im = 0;
        in2[i].re = g[i];
        in2[i].im = 0;
    }
    for (; i < size; i++)
    {
        in1[i].re = 0;
        in1[i].im = 0;
        in2[i].re = 0;
        in2[i].im = 0;
    }
    gmx_fft_1d(fft, GMX_FFT_FORWARD, in1, in1);
    gmx_fft_1d(fft, GMX_FFT_FORWARD, in2, in2);

    for (i = 0; i < size; i++)
    {
        complexConjugatMult(&in1[i], &in2[i]);
        in1[i].re /= size;
    }
    gmx_fft_1d(fft, GMX_FFT_BACKWARD, in1, in1);

    for (i = 0; i < n; i++)
    {
        corr[i] = in1[i].re;
    }

    sfree(in1);
    sfree(in2);
}

void cross_corr(int n, real f[], real g[], real corr[])
{
    gmx_fft_t fft;
    gmx_fft_init_1d(&fft, zeroPaddingSize(n), GMX_FFT_FLAG_CONSERVATIVE);
    cross_corr_low(n, f, g, corr, fft);
    gmx_fft_destroy(fft);
    gmx_fft_cleanup();
}

void many_cross_corr(int nFunc, int* nData, real** f, real** g, real** corr)
{
#pragma omp parallel
    // gmx_fft_t is not thread safe, so structure are allocated per thread.
    {
        int i;

#pragma omp for
        for (i = 0; i < nFunc; i++)
        {
            try
            {
                gmx_fft_t fft;
                gmx_fft_init_1d(&fft, zeroPaddingSize(nData[i]), GMX_FFT_FLAG_CONSERVATIVE);
                cross_corr_low(nData[i], f[i], g[i], corr[i], fft);
                gmx_fft_destroy(fft);
            }
            GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
        }
    }
    gmx_fft_cleanup();
}
