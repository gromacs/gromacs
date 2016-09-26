/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
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

/*! \libinternal \file
 *  \brief Defines the CUDA 3D-FFT functions for PME.
 *
 *  \author Aleksei Iupinov <a.yupinov@gmail.com>
 */

#ifndef PME3DFFT_CUH
#define PME3DFFT_CUH

#include "gmxpre.h"

#include <cufft.h>

#include "gromacs/fft/fft.h"

#include "pme-gpu-types.h"

/*! \brief \internal A 3D FFT class for performing R2C/C2R transforms */
class gmx_parallel_3dfft_gpu_t
{
    ivec          nDataReal;
    ivec          sizeReal;
    ivec          sizeComplex;

    cufftHandle   planR2C;
    cufftHandle   planC2R;
    cufftReal    *realGrid;
    cufftComplex *complexGrid;

    /* unused */
    ivec          localOffset;
    public:
        /*! \brief
         * Constructs CUDA FFT plans for performing 3D FFT on a PME grid.
         *
         * \param[in] pmeGPU                  The PME GPU structure.
         */
        gmx_parallel_3dfft_gpu_t(const pme_gpu_t *pmeGPU);
        /*! \brief Destroys CUDA FFT plans. */
        ~gmx_parallel_3dfft_gpu_t();
        /*! \brief
         * Returns the grid dimensions of the local real-space grid.
         *
         * \param[out]   localNData           The numbers of real elements in the local grid.
         * \param[out]   localOffset          The offsets of the local grid.
         * \param[out]   localSize            The numbers of real elements (with padding) in the local grid.
         */
        void get_real_limits(ivec localNData, ivec localOffset, ivec localSize);
        /*! \brief
         * Returns the grid dimensions of the local complex grid.
         *
         * \param[out]   localNData           The numbers of complex elements in the local grid.
         * \param[out]   localOffset          The offsets of the local grid.
         * \param[out]   localSize            The numbers of complex elements (with padding) in the local grid.
         */
        void get_complex_limits(ivec localNData, ivec localOffset, ivec localSize);
        /*! \brief
         * Performs the 3D FFT.
         *
         * \param[in] dir                     The transform direction.
         * \returns                           The cuFFT result code (0 if no error).
         */
        cufftResult_t perform_3dfft(gmx_fft_direction dir);
};

#endif // PME3DFFT_CUH
