/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
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
 * \brief Defines the host-side PME GPU data structure, which is dependent on the GPU types.
 * It's included by pointer in the general PmeGpu host structure in pme_gpu_types_host.h.
 *
 * \author Aleksei Iupinov <a.yupinov@gmail.com>
 * \ingroup module_ewald
 */

#ifndef PMEGPUTYPESHOSTIMPL_H
#define PMEGPUTYPESHOSTIMPL_H

#include "config.h"
#include "gromacs/utility/enumerationhelpers.h"

#include <array>
#include <set>
#include <vector>

#if GMX_GPU_CUDA
#    include "gromacs/gpu_utils/gpuregiontimer.cuh"
#elif GMX_GPU_OPENCL
#    include "gromacs/gpu_utils/gpuregiontimer_ocl.h"
#elif GMX_GPU_SYCL
#    include "gromacs/gpu_utils/gpuregiontimer_sycl.h"
#endif

#include "gromacs/gpu_utils/gpueventsynchronizer.h"

#include "gromacs/fft/gpu_3dfft.h"
#include "gromacs/timing/gpu_timing.h" // for gtPME_EVENT_COUNT

#ifndef NUMFEPSTATES
//! Number of FEP states.
#    define NUMFEPSTATES 2
#endif

namespace gmx
{
class Gpu3dFft;
} // namespace gmx

/*! \internal \brief
 * The main PME CUDA/OpenCL-specific host data structure, included in the PME GPU structure by the archSpecific pointer.
 */
struct PmeGpuSpecific
{
    /*! \brief Constructor
     *
     * \param[in] deviceContext  GPU device context
     * \param[in] pmeStream      GPU pme stream.
     */
    PmeGpuSpecific(const DeviceContext& deviceContext, const DeviceStream& pmeStream) :
        deviceContext_(deviceContext), pmeStream_(pmeStream)
    {
    }

    /*! \brief
     * A handle to the GPU context.
     * TODO: this is currently extracted from the implementation of pmeGpu->programHandle_,
     * but should be a constructor parameter to PmeGpu, as well as PmeGpuProgram,
     * managed by high-level code.
     */
    const DeviceContext& deviceContext_;

    /*! \brief The GPU stream where everything related to the PME happens. */
    const DeviceStream& pmeStream_;

    /* Synchronization events */
    /*! \brief Triggered after the PME Force Calculations have been completed */
    GpuEventSynchronizer pmeForcesReady;
    /*! \brief Triggered after the grid has been copied to the host (after the spreading stage). */
    GpuEventSynchronizer syncSpreadGridD2H;
    /*! \brief Triggered after the grid has been converted from FFT grid to PME grid (before the gather stage). */
    GpuEventSynchronizer syncFftToPmeGrid;
    /*! \brief Triggered after spline/spread computations have been completed. */
    GpuEventSynchronizer spreadCompleted;

    /* Settings which are set at the start of the run */
    /*! \brief A boolean which tells whether the complex and real grids for cu/clFFT are different or same. Currently true. */
    bool performOutOfPlaceFFT = false;
    /*! \brief A boolean which tells if the GPU timing events are enabled.
     *  False by default, can be enabled by setting the environment variable GMX_ENABLE_GPU_TIMING.
     *  Note: will not be reliable when multiple GPU tasks are running concurrently on the same
     * device context, as CUDA events on multiple streams are untrustworthy.
     */
    bool useTiming = false;

    //! Vector of FFT setups
    std::vector<std::unique_ptr<gmx::Gpu3dFft>> fftSetup;

    //! All the timers one might use
    gmx::EnumerationArray<PmeStage, GpuRegionTimer> timingEvents;

    //! Indices of timingEvents actually used
    std::set<PmeStage> activeTimers;

    /* GPU arrays element counts (not the arrays sizes in bytes!).
     * They might be larger than the actual meaningful data sizes.
     * These are paired: the actual element count + the maximum element count that can fit in the current allocated memory.
     * These integer pairs are mostly meaningful for the reallocateDeviceBuffer calls.
     * As such, if DeviceBuffer is refactored into a class, they can be freely changed, too.
     * The only exceptions are realGridSize and complexGridSize which are also used for grid clearing/copying.
     * TODO: these should live in a clean buffered container type, and be refactored in the NB/cudautils as well.
     */
    /*! \brief The kernelParams.atoms.coordinates float element count (actual)*/
    int coordinatesSize = 0;
    /*! \brief The kernelParams.atoms.coordinates float element count (reserved) */
    int coordinatesSizeAlloc = 0;
    /*! \brief The kernelParams.atoms.forces float element count (actual) */
    int forcesSize = 0;
    /*! \brief The kernelParams.atoms.forces float element count (reserved) */
    int forcesSizeAlloc = 0;
    /*! \brief The kernelParams.atoms.gridlineIndices int element count (actual) */
    int gridlineIndicesSize = 0;
    /*! \brief The kernelParams.atoms.gridlineIndices int element count (reserved) */
    int gridlineIndicesSizeAlloc = 0;
    /*! \brief Both the kernelParams.atoms.theta and kernelParams.atoms.dtheta float element count (actual) */
    int splineDataSize = 0;
    /*! \brief Both the kernelParams.atoms.theta and kernelParams.atoms.dtheta float element count (reserved) */
    int splineDataSizeAlloc = 0;
    /*! \brief The kernelParams.atoms.coefficients float element count (actual) */
    int coefficientsSize[NUMFEPSTATES] = { 0, 0 };
    /*! \brief The kernelParams.atoms.coefficients float element count (reserved) */
    int coefficientsCapacity[NUMFEPSTATES] = { 0, 0 };
    /*! \brief The kernelParams.grid.splineValuesArray float element count (actual) */
    int splineValuesSize[NUMFEPSTATES] = { 0, 0 };
    /*! \brief The kernelParams.grid.splineValuesArray float element count (reserved) */
    int splineValuesCapacity[NUMFEPSTATES] = { 0, 0 };
    /*! \brief The kernelParams.grid.realGrid float element count (actual) */
    int realGridSize[NUMFEPSTATES] = { 0, 0 };
    /*! \brief The kernelParams.grid.realGrid float element count (reserved) */
    int realGridCapacity[NUMFEPSTATES] = { 0, 0 };
    /*! \brief The kernelParams.grid.fourierGrid float (not float2!) element count (actual) */
    int complexGridSize[NUMFEPSTATES] = { 0, 0 };
    /*! \brief The kernelParams.grid.fourierGrid float (not float2!) element count (reserved) */
    int complexGridCapacity[NUMFEPSTATES] = { 0, 0 };

    /*! \brief Buffer size used to transfer PME grid overlap region in X-dimension*/
    int overlapXSizeLeft = 0;
    /*! \brief Buffer capacity used to transfer PME grid overlap region in X-dimension*/
    int overlapXCapacityLeft = 0;
    /*! \brief Buffer size used to transfer PME grid overlap region in X-dimension*/
    int overlapXSizeRight = 0;
    /*! \brief Buffer capacity used to transfer PME grid overlap region in X-dimension*/
    int overlapXCapacityRight = 0;
    /*! \brief Buffer capacity used to send PME grid overlap region in Y-dimension*/
    int overlapYSendSizeLeft = 0;
    /*! \brief Buffer capacity used to send PME grid overlap region in Y-dimension*/
    int overlapYSendCapacityLeft = 0;
    /*! \brief Buffer size used to recv PME grid overlap region in Y-dimension*/
    int overlapYRecvSizeLeft = 0;
    /*! \brief Buffer capacity used to recv PME grid overlap region in Y-dimension*/
    int overlapYRecvCapacityLeft = 0;
    /*! \brief Buffer capacity used to send PME grid overlap region in Y-dimension*/
    int overlapYSendSizeRight = 0;
    /*! \brief Buffer capacity used to send PME grid overlap region in Y-dimension*/
    int overlapYSendCapacityRight = 0;
    /*! \brief Buffer size used to recv PME grid overlap region in Y-dimension*/
    int overlapYRecvSizeRight = 0;
    /*! \brief Buffer capacity used to recv PME grid overlap region in Y-dimension*/
    int overlapYRecvCapacityRight = 0;
    /*! \brief Buffer used to transfer PME grid overlap region in X-dimension*/
    DeviceBuffer<float> d_recvGridLeftX = nullptr;
    /*! \brief Buffer used to transfer PME grid overlap region in X-dimension*/
    DeviceBuffer<float> d_recvGridRightX = nullptr;
    /*! \brief Buffer used to send PME grid overlap region in Y-dimension*/
    DeviceBuffer<float> d_sendGridLeftY = nullptr;
    /*! \brief Buffer used to recv PME grid overlap region in Y-dimension*/
    DeviceBuffer<float> d_recvGridLeftY = nullptr;
    /*! \brief Buffer used to send PME grid overlap region in Y-dimension*/
    DeviceBuffer<float> d_sendGridRightY = nullptr;
    /*! \brief Buffer used to recv PME grid overlap region in Y-dimension*/
    DeviceBuffer<float> d_recvGridRightY = nullptr;
};

#endif
