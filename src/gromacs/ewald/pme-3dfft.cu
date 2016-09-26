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

#include "pme-3dfft.cuh"

#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"

#include "pme.cuh"
#include "pme-gpu-types.h"

gmx_parallel_3dfft_gpu_t::gmx_parallel_3dfft_gpu_t(const pme_gpu_t *pmeGPU)
{
    for (int i = 0; i < DIM; i++)
    {
        _nDataReal[i]   = pmeGPU->kernelParams.grid.localGridSize[i];
        _sizeComplex[i] = _sizeReal[i] = pmeGPU->kernelParams.grid.localGridSizePadded[i];
    }
    if (!pmeGPU->archSpecific->bOutOfPlaceFFT)
    {
        GMX_ASSERT(_sizeComplex[ZZ] % 2 == 0, "Odd inplace cuFFT minor dimension");
    }
    _sizeComplex[ZZ] /= 2;

    GMX_ASSERT(!pme_gpu_uses_dd(pmeGPU), "FFT decomposition not implemented");

    const int gridSizeComplex = _sizeComplex[XX] * _sizeComplex[YY] * _sizeComplex[ZZ];
    const int gridSizeReal    = _sizeReal[XX] * _sizeReal[YY] * _sizeReal[ZZ];

    memset(_localOffset, 0, sizeof(_localOffset)); //!

    _realGrid = (cufftReal *)pmeGPU->kernelParams.grid.realGrid;
    assert(_realGrid);
    _complexGrid = (cufftComplex *)pmeGPU->kernelParams.grid.fourierGrid;

    /* Commented code for a simple 3D grid with no padding */
    /*
       result = cufftPlan3d(&_planR2C, _ndataReal[XX], _ndataReal[YY], _ndataReal[ZZ], CUFFT_R2C);
       if (result != CUFFT_SUCCESS)
       gmx_fatal(FARGS, "cufftPlan3d R2C error %d\n", result);

       result = cufftPlan3d(&_planC2R, _ndataReal[XX], _ndataReal[YY], _ndataReal[ZZ], CUFFT_C2R);
       if (result != CUFFT_SUCCESS)
       gmx_fatal(FARGS, "cufftPlan3d C2R error %d\n", result);
     */

    cufftResult_t             result;
    const int                 rank = 3, batch = 1;
    result = cufftPlanMany(&_planR2C, rank, _nDataReal,
                           _sizeReal, 1, gridSizeReal,
                           _sizeComplex, 1, gridSizeComplex,
                           CUFFT_R2C,
                           batch);
    if (result != CUFFT_SUCCESS)
    {
        gmx_fatal(FARGS, "cufftPlanMany R2C error %d\n", result);
    }

    result = cufftPlanMany(&_planC2R, rank, _nDataReal,
                           _sizeComplex, 1, gridSizeComplex,
                           _sizeReal, 1, gridSizeReal,
                           CUFFT_C2R,
                           batch);
    if (result != CUFFT_SUCCESS)
    {
        gmx_fatal(FARGS, "cufftPlanMany C2R error %d\n", result);
    }

    cudaStream_t s = pmeGPU->archSpecific->pmeStream;
    assert(s);
    result = cufftSetStream(_planR2C, s);
    if (result != CUFFT_SUCCESS)
    {
        gmx_fatal(FARGS, "cufftSetStream R2C error %d\n", result);
    }

    result = cufftSetStream(_planC2R, s);
    if (result != CUFFT_SUCCESS)
    {
        gmx_fatal(FARGS, "cufftSetStream C2R error %d\n", result);
    }
}

gmx_parallel_3dfft_gpu_t::~gmx_parallel_3dfft_gpu_t()
{
    cufftResult_t result;
    result = cufftDestroy(_planR2C);
    if (result != CUFFT_SUCCESS)
    {
        gmx_fatal(FARGS, "cufftDestroy R2C error %d\n", result);
    }
    result = cufftDestroy(_planC2R);
    if (result != CUFFT_SUCCESS)
    {
        gmx_fatal(FARGS, "cufftDestroy C2R error %d\n", result);
    }
}

void gmx_parallel_3dfft_gpu_t::get_real_limits(ivec localNData, ivec localOffset, ivec localSize)
{
    if (localNData)
    {
        memcpy(localNData, _nDataReal, sizeof(_nDataReal));
    }
    if (localSize)
    {
        memcpy(localSize, _sizeReal, sizeof(_sizeReal));
    }
    if (localOffset)
    {
        memcpy(localOffset, _localOffset, sizeof(_localOffset));
    }
}

void gmx_parallel_3dfft_gpu_t::get_complex_limits(ivec localNData, ivec localOffset, ivec localSize)
{
    if (localNData)
    {
        memcpy(localNData, _nDataReal, sizeof(_nDataReal));
        localNData[ZZ] = localNData[ZZ] / 2 + 1;
    }
    if (localSize)
    {
        memcpy(localSize, _sizeComplex, sizeof(_sizeComplex));
    }
    if (localOffset)
    {
        memcpy(localOffset, _localOffset, sizeof(_localOffset));
    }
}

cufftResult_t gmx_parallel_3dfft_gpu_t::perform_3dfft(gmx_fft_direction dir)
{
    cufftResult_t result;
    if (dir == GMX_FFT_REAL_TO_COMPLEX)
    {
        result = cufftExecR2C(_planR2C, _realGrid, _complexGrid);
    }
    else
    {
        result = cufftExecC2R(_planC2R, _complexGrid, _realGrid);
    }
    return result;
}

void pme_gpu_3dfft(const pme_gpu_t *pmeGPU, gmx_fft_direction dir, int grid_index)
{
    int           timerId = (dir == GMX_FFT_REAL_TO_COMPLEX) ? gtPME_FFT_R2C : gtPME_FFT_C2R;
    pme_gpu_start_timing(pmeGPU, timerId);
    cufftResult_t result = pmeGPU->archSpecific->pfft_setup_gpu[grid_index]->perform_3dfft(dir);
    pme_gpu_stop_timing(pmeGPU, timerId);
    if (result)
    {
        gmx_fatal(FARGS, "cuFFT %s error %d\n", (dir == GMX_FFT_REAL_TO_COMPLEX) ? "R2C" : "C2R", result);
    }
}
