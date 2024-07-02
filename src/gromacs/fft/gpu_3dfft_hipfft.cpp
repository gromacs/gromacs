/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2024- The GROMACS Authors
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
 *  \brief Implements GPU 3D FFT routines for HIP.
 *
 *  \author Aleksei Iupinov <a.yupinov@gmail.com>
 *  \author Mark Abraham <mark.j.abraham@gmail.com>
 *  \author Paul Bauer <paul.bauer.q@gmail.com>
 *  \ingroup module_fft
 */

#include "gmxpre.h"

#include "gpu_3dfft_hipfft.h"

#include "gromacs/gpu_utils/device_stream.h"
#include "gromacs/gpu_utils/devicebuffer.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"

namespace gmx
{
static void handleHipfftError(hipfftResult_t status, const char* msg)
{
    if (status != HIPFFT_SUCCESS)
    {
        gmx_fatal(FARGS, "%s (error code %d)\n", msg, status);
    }
}

Gpu3dFft::ImplHipFft::ImplHipFft(bool allocateRealGrid,
                                 MPI_Comm /*comm*/,
                                 ArrayRef<const int> gridSizesInXForEachRank,
                                 ArrayRef<const int> gridSizesInYForEachRank,
                                 const int /*nz*/,
                                 bool                 performOutOfPlaceFFT,
                                 const DeviceContext& context,
                                 const DeviceStream&  pmeStream,
                                 ivec                 realGridSize,
                                 ivec                 realGridSizePadded,
                                 ivec                 complexGridSizePadded,
                                 DeviceBuffer<float>* realGrid,
                                 DeviceBuffer<float>* complexGrid) :
    Gpu3dFft::Impl::Impl(performOutOfPlaceFFT), realGrid_(reinterpret_cast<hipfftReal*>(*realGrid))
{
    GMX_RELEASE_ASSERT(allocateRealGrid == false, "Grids needs to be pre-allocated");
    GMX_RELEASE_ASSERT(gridSizesInXForEachRank.size() == 1 && gridSizesInYForEachRank.size() == 1,
                       "FFT decomposition not implemented with hipFFT backend");

    allocateComplexGrid(complexGridSizePadded, realGrid, complexGrid, context);

    const int complexGridSizePaddedTotal =
            complexGridSizePadded[XX] * complexGridSizePadded[YY] * complexGridSizePadded[ZZ];
    const int realGridSizePaddedTotal =
            realGridSizePadded[XX] * realGridSizePadded[YY] * realGridSizePadded[ZZ];

    GMX_RELEASE_ASSERT(realGrid_, "Bad (null) input real-space grid");
    GMX_RELEASE_ASSERT(complexGrid_, "Bad (null) input complex grid");

    hipfftResult_t result;

    const int rank = 3, batch = 1;
    result = hipfftPlanMany(&planR2C_,
                            rank,
                            realGridSize,
                            realGridSizePadded,
                            1,
                            realGridSizePaddedTotal,
                            complexGridSizePadded,
                            1,
                            complexGridSizePaddedTotal,
                            HIPFFT_R2C,
                            batch);
    handleHipfftError(result, "hipfftPlanMany R2C plan failure");

    result = hipfftPlanMany(&planC2R_,
                            rank,
                            realGridSize,
                            complexGridSizePadded,
                            1,
                            complexGridSizePaddedTotal,
                            realGridSizePadded,
                            1,
                            realGridSizePaddedTotal,
                            HIPFFT_C2R,
                            batch);
    handleHipfftError(result, "hipfftPlanMany C2R plan failure");

    hipStream_t stream = pmeStream.stream();
    GMX_RELEASE_ASSERT(stream, "Can not use the default HIP stream for PME hipFFT");

    result = hipfftSetStream(planR2C_, stream);
    handleHipfftError(result, "hipfftSetStream R2C failure");

    result = hipfftSetStream(planC2R_, stream);
    handleHipfftError(result, "hipfftSetStream C2R failure");
}

Gpu3dFft::ImplHipFft::~ImplHipFft()
{
    deallocateComplexGrid();

    hipfftResult_t result;
    result = hipfftDestroy(planR2C_);
    handleHipfftError(result, "hipfftDestroy R2C failure");
    result = hipfftDestroy(planC2R_);
    handleHipfftError(result, "hipfftDestroy C2R failure");
}

void Gpu3dFft::ImplHipFft::perform3dFft(gmx_fft_direction dir, CommandEvent* /*timingEvent*/)
{
    hipfftResult_t result;
    if (dir == GMX_FFT_REAL_TO_COMPLEX)
    {
        result = hipfftExecR2C(planR2C_, realGrid_, (hipfftComplex*)complexGrid_);
        handleHipfftError(result, "hipFFT R2C execution failure");
    }
    else
    {
        result = hipfftExecC2R(planC2R_, (hipfftComplex*)complexGrid_, realGrid_);
        handleHipfftError(result, "hipFFT C2R execution failure");
    }
}

} // namespace gmx
