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
/*! \internal \file
 * \brief
 * Implements function to compute many autocorrelation functions
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_correlationfunctions
 */
#include "gmxpre.h"

#include "manyautocorrelation.h"

#include <cstddef>

#include <algorithm>

#include "gromacs/fft/fft.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxomp.h"

int many_auto_correl(std::vector<std::vector<real>>* c)
{
    size_t nfunc = (*c).size();
    if (nfunc == 0)
    {
        GMX_THROW(gmx::InconsistentInputError("Empty array of vectors supplied"));
    }
    size_t ndata = (*c)[0].size();
    if (ndata == 0)
    {
        GMX_THROW(gmx::InconsistentInputError("Empty vector supplied"));
    }
#ifndef NDEBUG
    for (size_t i = 1; i < nfunc; i++)
    {
        if ((*c)[i].size() != ndata)
        {
            char buf[256];
            snprintf(buf,
                     sizeof(buf),
                     "Vectors of different lengths supplied (%d %d)",
                     static_cast<int>((*c)[i].size()),
                     static_cast<int>(ndata));
            GMX_THROW(gmx::InconsistentInputError(buf));
        }
    }
#endif
    // Add buffer size to the arrays.
    size_t nfft = (3 * ndata / 2) + 1;
    // Pad arrays with zeros
    for (auto& i : *c)
    {
        i.resize(nfft, 0);
    }
#pragma omp parallel
    {
        try
        {
            gmx_fft_t         fft1;
            std::vector<real> in, out;

            int nthreads  = gmx_omp_get_max_threads();
            int thread_id = gmx_omp_get_thread_num();
            int i0        = (thread_id * nfunc) / nthreads;
// nvc++ 24.1+ version has bug due to which it generates incorrect OMP code for this region
// with std::min() so we add a macro for min() function.
#if defined(__NVCOMPILER)
#    define min(l, r) (l < r ? l : r)
            int i1 = min(nfunc, ((thread_id + 1) * nfunc) / nthreads);
#    undef min
#else
            int i1 = std::min(nfunc, ((thread_id + 1) * nfunc) / nthreads);
#endif


            gmx_fft_init_1d(&fft1, nfft, GMX_FFT_FLAG_CONSERVATIVE);
            /* Allocate temporary arrays */
            in.resize(2 * nfft, 0);
            out.resize(2 * nfft, 0);
            for (int i = i0; (i < i1); i++)
            {
                for (size_t j = 0; j < ndata; j++)
                {
                    in[2 * j + 0] = (*c)[i][j];
                    in[2 * j + 1] = 0;
                }
                gmx_fft_1d(fft1, GMX_FFT_BACKWARD, in.data(), out.data());
                for (size_t j = 0; j < nfft; j++)
                {
                    in[2 * j + 0] =
                            (out[2 * j + 0] * out[2 * j + 0] + out[2 * j + 1] * out[2 * j + 1]) / nfft;
                    in[2 * j + 1] = 0;
                }
                gmx_fft_1d(fft1, GMX_FFT_FORWARD, in.data(), out.data());
                for (size_t j = 0; (j < nfft); j++)
                {
                    (*c)[i][j] = out[2 * j + 0];
                }
            }
            /* Free the memory */
            gmx_fft_destroy(fft1);
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
    }
    for (auto& i : *c)
    {
        i.resize(ndata);
    }

    return 0;
}
