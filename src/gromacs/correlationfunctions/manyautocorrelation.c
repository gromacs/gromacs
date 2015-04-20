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
/*! \internal \file
 * \brief
 * Implements function to compute many autocorrelation functions
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_correlationfunctions
 */
#include "gmxpre.h"

#include "manyautocorrelation.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "gromacs/fft/fft.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/utility/gmxomp.h"
#include "gromacs/utility/smalloc.h"

int many_auto_correl(int nfunc, int ndata, int nfft, real **c)
{
    #pragma omp parallel
    {
        typedef real complex[2];
        int          i, t, j, fftcode;
        gmx_fft_t    fft1;
        complex     *in, *out;
        int          i0, i1;
        int          nthreads, thread_id;

        nthreads  = gmx_omp_get_max_threads();
        thread_id = gmx_omp_get_thread_num();
        if ((0 == thread_id))
        {
            // fprintf(stderr, "There are %d threads for correlation functions\n", nthreads);
        }
        i0 = thread_id*nfunc/nthreads;
        i1 = min(nfunc, (thread_id+1)*nfunc/nthreads);

        fftcode = gmx_fft_init_1d(&fft1, nfft, GMX_FFT_FLAG_CONSERVATIVE);
        /* Allocate temporary arrays */
        snew(in, nfft);
        snew(out, nfft);
        for (i = i0; (i < i1); i++)
        {
            for (j = 0; j < ndata; j++)
            {
                in[j][0] = c[i][j];
                in[j][1] = 0;
            }
            for (; (j < nfft); j++)
            {
                in[j][0] = in[j][1] = 0;
            }

            fftcode = gmx_fft_1d(fft1, GMX_FFT_BACKWARD, (void *)in, (void *)out);
            for (j = 0; j < nfft; j++)
            {
                in[j][0] = (out[j][0]*out[j][0] + out[j][1]*out[j][1])/nfft;
                in[j][1] = 0;
            }
            for (; (j < nfft); j++)
            {
                in[j][0] = in[j][1] = 0;
            }

            fftcode = gmx_fft_1d(fft1, GMX_FFT_FORWARD, (void *)in, (void *)out);
            for (j = 0; (j < nfft); j++)
            {
                c[i][j] = out[j][0]/ndata;
            }
        }
        /* Free the memory */
        gmx_fft_destroy(fft1);
        sfree(in);
        sfree(out);
    }
    // gmx_fft_cleanup();
    return 0;
}
