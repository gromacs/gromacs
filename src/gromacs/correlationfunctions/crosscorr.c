/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014, by the GROMACS development team, led by
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
#include "crosscorr.h"
#include <gromacs/fft/fft.h>
#include "gromacs/fft/fft.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/legacyheaders/macros.h"


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

typedef struct {
    double real;
    double img;
}  complex;

//return the size witch the array should be after zero padding
int zeroPadingSize(int n)
{
    return 2*n;
}

complex complexConjugatMult(complex in1, complex in2)
{
    complex res;
    res.real = in1.real * in2.real + in1.img * in2.img;
    res.img  = in1.real * -in2.img + in1.img * in2.real;
    return res;
}



void cross_corr_low(int n, real f[], real g[], real corr[], gmx_fft_t fft)
{
    int           i;
    const int     size = zeroPadingSize(n);
    complex       in1[size];
    complex       in2[size];

    for (i = 0; i < n; i++)
    {
        in1[i].real = (double)f[i];
        in1[i].img  = 0;
        in2[i].real = (double)g[i];
        in2[i].img  = 0;
    }
    for (; i < size; i++)
    {
        in1[i].real = 0;
        in1[i].img  = 0;
        in2[i].real = 0;
        in2[i].img  = 0;
    }


    gmx_fft_1d(fft, GMX_FFT_FORWARD, in1, in1);
    gmx_fft_1d(fft, GMX_FFT_FORWARD, in2, in2);

    for (i = 0; i < size; i++)
    {
        in1[i]       = complexConjugatMult(in1[i], in2[i]);
        in1[i].real /= size;
    }
    gmx_fft_1d(fft, GMX_FFT_BACKWARD, in1, in1);

    for (i = 0; i < n; i++)
    {
        corr[i] = (real)(in1[i].real);
    }
}



void cross_corr(int n, real f[], real g[], real corr[])
{
    gmx_fft_t fft;
    gmx_fft_init_1d(&fft, zeroPadingSize(n), GMX_FFT_FLAG_CONSERVATIVE);
    cross_corr_low( n,  f,  g, corr, fft);
    gmx_fft_destroy(fft);
    gmx_fft_cleanup();
}

void multi_cross_corr(int nFunc, int * nData, real ** f, real ** g, real ** corr)
{
#pragma omp parallel //gmx_fft_t is not thread safe
    {
        int i;

#pragma omp for
        for (i = 0; i < nFunc; i++)
        {
            gmx_fft_t fft;
            gmx_fft_init_1d(&fft, zeroPadingSize(nData[i]), GMX_FFT_FLAG_CONSERVATIVE);
            cross_corr_low( nData[i],  f[i],  g[i], corr[i], fft);
            gmx_fft_destroy(fft);
        }
    }
    gmx_fft_cleanup();

}
