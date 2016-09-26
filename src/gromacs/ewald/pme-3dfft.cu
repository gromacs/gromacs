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

/*! \internal \file
 *  \brief Implements CUDA FFT routines for PME GPU.
 *
 *  \author Aleksei Iupinov <a.yupinov@gmail.com>
 */

#include "gmxpre.h"

#include <assert.h>

#include <cufft.h>

#include "gromacs/gpu_utils/cudautils.cuh"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

#include "pme.cuh"
#include "pme-gpu.h"

struct gmx_parallel_3dfft_gpu
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
};

void pme_gpu_init_3dfft_plan(gmx_parallel_3dfft_gpu_t *pfftSetup, const pme_gpu_t *pmeGPU)
{
    cufftResult_t            result;
    gmx_parallel_3dfft_gpu_t setup;
    snew(setup, 1);
    *pfftSetup = setup;

    for (int i = 0; i < DIM; i++)
    {
        setup->nDataReal[i]   = pmeGPU->kernelParams.grid.localGridSize[i];
        setup->sizeComplex[i] = setup->sizeReal[i] = pmeGPU->kernelParams.grid.localGridSizePadded[i];
    }
    if (!pmeGPU->archSpecific->bOutOfPlaceFFT)
    {
        GMX_ASSERT(setup->sizeComplex[ZZ] % 2 == 0, "Odd inplace cuFFT minor dimension");
    }
    setup->sizeComplex[ZZ] /= 2;

    GMX_ASSERT(!pme_gpu_uses_dd(pmeGPU), "FFT decomposition not implemented");

    const int gridSizeComplex = setup->sizeComplex[XX] * setup->sizeComplex[YY] * setup->sizeComplex[ZZ];
    const int gridSizeReal    = setup->sizeReal[XX] * setup->sizeReal[YY] * setup->sizeReal[ZZ];

    memset(setup->localOffset, 0, sizeof(setup->localOffset)); //!

    setup->realGrid = (cufftReal *)pmeGPU->kernelParams.grid.realGrid;
    assert(setup->realGrid);
    setup->complexGrid = (cufftComplex *)pmeGPU->kernelParams.grid.fourierGrid;

    /* Commented code for a simple 3D grid with no padding */
    /*
       result = cufftPlan3d(&setup->planR2C, setup->ndataReal[XX], setup->ndataReal[YY], setup->ndataReal[ZZ], CUFFT_R2C);
       if (result != CUFFT_SUCCESS)
       gmx_fatal(FARGS, "cufftPlan3d R2C error %d\n", result);

       result = cufftPlan3d(&setup->planC2R, setup->ndataReal[XX], setup->ndataReal[YY], setup->ndataReal[ZZ], CUFFT_C2R);
       if (result != CUFFT_SUCCESS)
       gmx_fatal(FARGS, "cufftPlan3d C2R error %d\n", result);
     */

    const int rank = 3, batch = 1;
    result = cufftPlanMany(&setup->planR2C, rank, setup->nDataReal,
                           setup->sizeReal, 1, gridSizeReal,
                           setup->sizeComplex, 1, gridSizeComplex,
                           CUFFT_R2C,
                           batch);
    if (result != CUFFT_SUCCESS)
    {
        gmx_fatal(FARGS, "cufftPlanMany R2C error %d\n", result);
    }

    result = cufftPlanMany(&setup->planC2R, rank, setup->nDataReal,
                           setup->sizeComplex, 1, gridSizeComplex,
                           setup->sizeReal, 1, gridSizeReal,
                           CUFFT_C2R,
                           batch);
    if (result != CUFFT_SUCCESS)
    {
        gmx_fatal(FARGS, "cufftPlanMany C2R error %d\n", result);
    }

    cudaStream_t s = pmeGPU->archSpecific->pmeStream;
    assert(s);
    result = cufftSetStream(setup->planR2C, s);
    if (result != CUFFT_SUCCESS)
    {
        gmx_fatal(FARGS, "cufftSetStream R2C error %d\n", result);
    }

    result = cufftSetStream(setup->planC2R, s);
    if (result != CUFFT_SUCCESS)
    {
        gmx_fatal(FARGS, "cufftSetStream C2R error %d\n", result);
    }
}

void pme_gpu_get_3dfft_real_limits(gmx_parallel_3dfft_gpu_t      setup,
                                   ivec                          localNData,
                                   ivec                          localOffset,
                                   ivec                          localSize)
{
    if (localNData)
    {
        memcpy(localNData, setup->nDataReal, sizeof(setup->nDataReal));
    }
    if (localSize)
    {
        memcpy(localSize, setup->sizeReal, sizeof(setup->sizeReal));
    }
    if (localOffset)
    {
        memcpy(localOffset, setup->localOffset, sizeof(setup->localOffset));
    }
}

void pme_gpu_get_3dfft_complex_limits(const gmx_parallel_3dfft_gpu_t setup,
                                      ivec                           localNData,
                                      ivec                           localOffset,
                                      ivec                           localSize)
{
    if (localNData)
    {
        memcpy(localNData, setup->nDataReal, sizeof(setup->nDataReal));
        localNData[ZZ] = localNData[ZZ] / 2 + 1;
    }
    if (localSize)
    {
        memcpy(localSize, setup->sizeComplex, sizeof(setup->sizeComplex));
    }
    if (localOffset)
    {
        memcpy(localOffset, setup->localOffset, sizeof(setup->localOffset));
    }
}

void pme_gpu_3dfft(const pme_gpu_t        *pmeGPU,
                   gmx_fft_direction       dir,
                   const int               grid_index)
{
    gmx_parallel_3dfft_gpu_t setup = pmeGPU->archSpecific->pfft_setup_gpu[grid_index];

    if (dir == GMX_FFT_REAL_TO_COMPLEX)
    {
        pme_gpu_start_timing(pmeGPU, gtPME_FFT_R2C);
        cufftResult_t result = cufftExecR2C(setup->planR2C, setup->realGrid, setup->complexGrid);
        pme_gpu_stop_timing(pmeGPU, gtPME_FFT_R2C);
        if (result)
        {
            gmx_fatal(FARGS, "cufft R2C error %d\n", result);
        }
    }
    else
    {
        pme_gpu_start_timing(pmeGPU, gtPME_FFT_C2R);
        cufftResult_t result = cufftExecC2R(setup->planC2R, setup->complexGrid, setup->realGrid);
        pme_gpu_stop_timing(pmeGPU, gtPME_FFT_C2R);
        if (result)
        {
            gmx_fatal(FARGS, "cufft C2R error %d\n", result);
        }
    }
}

void pme_gpu_destroy_3dfft_plan(const gmx_parallel_3dfft_gpu_t &pfftSetup)
{
    if (pfftSetup)
    {
        cufftResult_t result;

        result = cufftDestroy(pfftSetup->planR2C);
        if (result != CUFFT_SUCCESS)
        {
            gmx_fatal(FARGS, "cufftDestroy R2C error %d\n", result);
        }
        result = cufftDestroy(pfftSetup->planC2R);
        if (result != CUFFT_SUCCESS)
        {
            gmx_fatal(FARGS, "cufftDestroy C2R error %d\n", result);
        }

        sfree(pfftSetup);
    }
}
