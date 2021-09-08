/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2017,2018,2019,2020,2021, by the GROMACS development team, led by
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
 *  \brief Implements GPU 3D FFT routines for CUDA.
 *
 *  \author Aleksei Iupinov <a.yupinov@gmail.com>
 *  \author Mark Abraham <mark.j.abraham@gmail.com>
 *  \ingroup module_fft
 */

#include "gmxpre.h"

#include "gpu_3dfft_cufft.h"

#include "gromacs/gpu_utils/device_stream.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"

namespace gmx
{
static void handleCufftError(cufftResult_t status, const char* msg)
{
    if (status != CUFFT_SUCCESS)
    {
        gmx_fatal(FARGS, "%s (error code %d)\n", msg, status);
    }
}

Gpu3dFft::ImplCuFft::ImplCuFft(bool allocateGrids,
                               MPI_Comm /*comm*/,
                               ArrayRef<const int> gridSizesInXForEachRank,
                               ArrayRef<const int> gridSizesInYForEachRank,
                               const int /*nz*/,
                               bool /*performOutOfPlaceFFT*/,
                               const DeviceContext& /*context*/,
                               const DeviceStream&  pmeStream,
                               ivec                 realGridSize,
                               ivec                 realGridSizePadded,
                               ivec                 complexGridSizePadded,
                               DeviceBuffer<float>* realGrid,
                               DeviceBuffer<float>* complexGrid) :
    realGrid_(reinterpret_cast<cufftReal*>(*realGrid)),
    complexGrid_(reinterpret_cast<cufftComplex*>(*complexGrid))
{
    GMX_RELEASE_ASSERT(allocateGrids == false, "Grids needs to be pre-allocated");
    GMX_RELEASE_ASSERT(gridSizesInXForEachRank.size() == 1 && gridSizesInYForEachRank.size() == 1,
                       "FFT decomposition not implemented with cuFFT backend");

    const int complexGridSizePaddedTotal =
            complexGridSizePadded[XX] * complexGridSizePadded[YY] * complexGridSizePadded[ZZ];
    const int realGridSizePaddedTotal =
            realGridSizePadded[XX] * realGridSizePadded[YY] * realGridSizePadded[ZZ];

    GMX_RELEASE_ASSERT(realGrid_, "Bad (null) input real-space grid");
    GMX_RELEASE_ASSERT(complexGrid_, "Bad (null) input complex grid");

    cufftResult_t result;
    /* Commented code for a simple 3D grid with no padding */
    /*
       result = cufftPlan3d(&planR2C_, realGridSize[XX], realGridSize[YY], realGridSize[ZZ],
       CUFFT_R2C); handleCufftError(result, "cufftPlan3d R2C plan failure");

       result = cufftPlan3d(&planC2R_, realGridSize[XX], realGridSize[YY], realGridSize[ZZ],
       CUFFT_C2R); handleCufftError(result, "cufftPlan3d C2R plan failure");
     */

    const int rank = 3, batch = 1;
    result = cufftPlanMany(&planR2C_,
                           rank,
                           realGridSize,
                           realGridSizePadded,
                           1,
                           realGridSizePaddedTotal,
                           complexGridSizePadded,
                           1,
                           complexGridSizePaddedTotal,
                           CUFFT_R2C,
                           batch);
    handleCufftError(result, "cufftPlanMany R2C plan failure");

    result = cufftPlanMany(&planC2R_,
                           rank,
                           realGridSize,
                           complexGridSizePadded,
                           1,
                           complexGridSizePaddedTotal,
                           realGridSizePadded,
                           1,
                           realGridSizePaddedTotal,
                           CUFFT_C2R,
                           batch);
    handleCufftError(result, "cufftPlanMany C2R plan failure");

    cudaStream_t stream = pmeStream.stream();
    GMX_RELEASE_ASSERT(stream, "Can not use the default CUDA stream for PME cuFFT");

    result = cufftSetStream(planR2C_, stream);
    handleCufftError(result, "cufftSetStream R2C failure");

    result = cufftSetStream(planC2R_, stream);
    handleCufftError(result, "cufftSetStream C2R failure");
}

Gpu3dFft::ImplCuFft::~ImplCuFft()
{
    cufftResult_t result;
    result = cufftDestroy(planR2C_);
    handleCufftError(result, "cufftDestroy R2C failure");
    result = cufftDestroy(planC2R_);
    handleCufftError(result, "cufftDestroy C2R failure");
}

void Gpu3dFft::ImplCuFft::perform3dFft(gmx_fft_direction dir, CommandEvent* /*timingEvent*/)
{
    cufftResult_t result;
    if (dir == GMX_FFT_REAL_TO_COMPLEX)
    {
        result = cufftExecR2C(planR2C_, realGrid_, complexGrid_);
        handleCufftError(result, "cuFFT R2C execution failure");
    }
    else
    {
        result = cufftExecC2R(planC2R_, complexGrid_, realGrid_);
        handleCufftError(result, "cuFFT C2R execution failure");
    }
}

} // namespace gmx
