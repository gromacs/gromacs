/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
#include "gmxpre.h"

#include "config.h"

#include "gromacs/math/gmxcomplex.h"
#include "gromacs/utility/gmxassert.h"

#include "gromacs/gpu_utils/gputraits_ocl.h"
//#include "gromacs/ewald/pme-ocl-types-kernel.clh"
//#include "gromacs/ewald/pme-ocl.h"
#include "pme-types-ocl.h"

#include "pme-ocl-definitely-common.h"

//#include "pme-timings.cuh"

void pme_gpu_solve(PmeGpu *pmeGpu, t_complex *h_grid,
                   GridOrdering gridOrdering, bool computeEnergyAndVirial)
{
    const bool    copyInputAndOutputGrid = pme_gpu_is_testing(pmeGpu) || !pme_gpu_performs_FFT(pmeGpu);

    CommandStream stream          = pmeGpu->archSpecific->pmeStream;
    auto         *kernelParamsPtr = pmeGpu->kernelParams.get();

    if (copyInputAndOutputGrid)
    {
        copyToDeviceBuffer(&kernelParamsPtr->grid.d_fourierGrid, (float *)(h_grid), 0, pmeGpu->archSpecific->complexGridSize,
                           stream, pmeGpu->settings.transferKind, nullptr);
    }

    int majorDim = -1, middleDim = -1, minorDim = -1;
    switch (gridOrdering)
    {
        case GridOrdering::YZX:
            majorDim  = YY;
            middleDim = ZZ;
            minorDim  = XX;
            break;

        case GridOrdering::XYZ:
            majorDim  = XX;
            middleDim = YY;
            minorDim  = ZZ;
            break;

        default:
            GMX_ASSERT(false, "Implement grid ordering here and below for the kernel launch");
    }

    const int          maxBlockSize      = c_solveMaxThreadsPerBlock;
    const int          gridLineSize      = kernelParamsPtr->grid.complexGridSize[minorDim];
    const int          gridLinesPerBlock = std::max(maxBlockSize / gridLineSize, 1);
    const int          blocksPerGridLine = (gridLineSize + maxBlockSize - 1) / maxBlockSize;
    const int          cellsPerBlock     = gridLineSize * gridLinesPerBlock;
    const int          blockSize         = (cellsPerBlock + warp_size - 1) / warp_size * warp_size;
    // rounding up to full warps so that shuffle operations produce defined results //FIXME
    KernelLaunchConfig config;
    config.blockSize.x = blockSize;
    config.gridSize.x  = blocksPerGridLine;
    config.gridSize.y  = (kernelParamsPtr->grid.complexGridSize[middleDim] + gridLinesPerBlock - 1) / gridLinesPerBlock;
    config.gridSize.z  = kernelParamsPtr->grid.complexGridSize[majorDim];
    config.stream      = stream;


    cl_kernel kernel;

    if (gridOrdering == GridOrdering::YZX)
    {
        if (computeEnergyAndVirial)
        {
            kernel = pmeGpu->archSpecific->persistent->solveYZXEnergyKernel;
        }
        else
        {
            kernel = pmeGpu->archSpecific->persistent->solveYZXKernel;
        }
    }
    else if (gridOrdering == GridOrdering::XYZ)
    {
        if (computeEnergyAndVirial)
        {
            kernel = pmeGpu->archSpecific->persistent->solveXYZEnergyKernel;
        }
        else
        {
            kernel = pmeGpu->archSpecific->persistent->solveXYZKernel;
        }
    }

    pme_gpu_start_timing(pmeGpu, gtPME_SOLVE);

    launchGpuKernel(config, kernel,
                    kernelParamsPtr,
                    &kernelParamsPtr->grid.d_splineModuli,
                    &kernelParamsPtr->constants.d_virialAndEnergy,
                    &kernelParamsPtr->grid.d_fourierGrid);
    //FIXME CU_LAUNCH_ERR("pme_solve_kernel");
    pme_gpu_stop_timing(pmeGpu, gtPME_SOLVE);

    if (computeEnergyAndVirial)
    {
        copyFromDeviceBuffer(pmeGpu->staging.h_virialAndEnergy, &kernelParamsPtr->constants.d_virialAndEnergy,
                             0, c_virialAndEnergyCount, stream, pmeGpu->settings.transferKind, nullptr);
    }

    if (copyInputAndOutputGrid)
    {
        copyFromDeviceBuffer( (float *)h_grid, &kernelParamsPtr->grid.d_fourierGrid,
                              0, pmeGpu->archSpecific->complexGridSize, stream, pmeGpu->settings.transferKind, nullptr);
    }
}
