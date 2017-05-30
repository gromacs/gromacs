/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2017, by the GROMACS development team, led by
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
 * \brief This file defines the PME CUDA-specific data structure,
 * various compile-time constants shared among the PME CUDA kernels,
 * and also names some PME CUDA memory management routines.
 * TODO: consider changing defines into variables where possible; have inline getters.
 *
 * \author Aleksei Iupinov <a.yupinov@gmail.com>
 */

#ifndef GMX_EWALD_PME_CUH
#define GMX_EWALD_PME_CUH

#include <cassert>

#include <array>
#include <set>

#include "gromacs/gpu_utils/cuda_arch_utils.cuh" // for warp_size

#include "pme-gpu-internal.h"                    // for the general PME GPU behaviour defines
#include "pme-timings.cuh"

class GpuParallel3dFft;

/* Some defines for PME behaviour follow */

/*
    Here is a current memory layout for the theta/dtheta B-spline float parameter arrays.
    This is the data in global memory used both by spreading and gathering kernels (with same scheduling).
    This example has PME order 4 and 2 particles per warp/data chunk.
    Each particle has 16 threads assigned to it, each thread works on 4 non-sequential global grid contributions.

    ----------------------------------------------------------------------------
    particles 0, 1                                        | particles 2, 3     | ...
    ----------------------------------------------------------------------------
    order index 0           | index 1 | index 2 | index 3 | order index 0 .....
    ----------------------------------------------------------------------------
    tx0 tx1 ty0 ty1 tz0 tz1 | ..........
    ----------------------------------------------------------------------------

    Each data chunk for a single warp is 24 floats. This goes both for theta and dtheta.
    24 = 2 particles per warp * order 4 * 3 dimensions. 48 floats (1.5 warp size) per warp in total.
    I have also tried intertwining theta and theta in a single array (they are used in pairs in gathering stage anyway)
    and it didn't seem to make a performance difference.

    The corresponding defines follow.
 */

/* This is the distance between the neighbour theta elements - would be 2 for the intertwining layout */
#define PME_SPLINE_THETA_STRIDE 1

/*! \brief
 * The number of GPU threads used for computing spread/gather contributions of a single atom as function of the PME order.
 * The assumption is currently that any thread processes only a single atom's contributions.
 */
#define PME_SPREADGATHER_THREADS_PER_ATOM (order * order)

/*! \brief
 * The number of atoms processed by a single warp in spread/gather.
 * This macro depends on the templated order parameter (2 atoms per warp for order 4).
 * It is mostly used for spline data layout tweaked for coalesced access.
 */
#define PME_SPREADGATHER_ATOMS_PER_WARP (warp_size / PME_SPREADGATHER_THREADS_PER_ATOM)

/*! \brief
 * Atom data alignment (in terms of number of atoms).
 * If the GPU atom data buffers are padded (c_usePadding == true),
 * Then the numbers of atoms which would fit in the padded GPU buffers has to be divisible by this.
 * The literal number (16) expresses maximum spread/gather block width in warps.
 * Accordingly, spread and gather block widths in warps should be divisors of this
 * (e.g. in the pme-spread.cu: constexpr int c_spreadMaxThreadsPerBlock = 8 * warp_size;).
 * There are debug asserts for this divisibility.
 */
#define PME_ATOM_DATA_ALIGNMENT (16 * PME_SPREADGATHER_ATOMS_PER_WARP)

/*! \brief \internal
 * An inline CUDA function for checking the global atom data indices against the atom data array sizes.
 *
 * \param[in] atomDataIndexGlobal  The atom data index.
 * \param[in] nAtomData            The atom data array element count.
 * \returns                        Non-0 if index is within bounds (or PME data padding is enabled), 0 otherwise.
 *
 * This is called from the spline_and_spread and gather PME kernels.
 * The goal is to isolate the global range checks, and allow avoiding them with c_usePadding enabled.
 */
int __device__ __forceinline__ pme_gpu_check_atom_data_index(const int atomDataIndex, const int nAtomData)
{
    return c_usePadding ? 1 : (atomDataIndex < nAtomData);
}

/*! \brief \internal
 * An inline CUDA function for skipping the zero-charge atoms.
 *
 * \returns                        Non-0 if atom should be processed, 0 otherwise.
 * \param[in] coefficient          The atom charge.
 *
 * This is called from the spline_and_spread and gather PME kernels.
 */
int __device__ __forceinline__ pme_gpu_check_atom_charge(const float coefficient)
{
    assert(isfinite(coefficient));
    return c_skipNeutralAtoms ? (coefficient != 0.0f) : 1;
}

/*! \brief \internal
 * The main PME CUDA-specific host data structure, included in the PME GPU structure by the archSpecific pointer.
 */
struct pme_gpu_cuda_t
{
    /*! \brief The CUDA stream where everything related to the PME happens. */
    cudaStream_t pmeStream;

    /* Synchronization events */
    /*! \brief Triggered after the energy/virial have been copied to the host (after the solving stage). */
    cudaEvent_t syncEnerVirD2H;
    /*! \brief Triggered after the output forces have been copied to the host (after the gathering stage). */
    cudaEvent_t syncForcesD2H;
    /*! \brief Triggered after the grid has been copied to the host (after the spreading stage). */
    cudaEvent_t syncSpreadGridD2H;
    /*! \brief Triggered after the atom spline data has been copied to the host (after the spline computation). */
    cudaEvent_t syncSplineAtomDataD2H;
    /*! \brief Triggered after the grid hes been copied to the host (after the solving stage) */
    cudaEvent_t syncSolveGridD2H;

    // TODO: consider moving some things below into the non-CUDA struct.

    /* Settings which are set at the start of the run */
    /*! \brief A boolean which tells whether the complex and real grids for cuFFT are different or same. Currenty true. */
    bool performOutOfPlaceFFT;
    /*! \brief A boolean which tells if the CUDA timing events are enabled.
     * True by default, disabled by setting the environment variable GMX_DISABLE_CUDA_TIMING.
     * FIXME: this should also be disabled if any other GPU task is running concurrently on the same device,
     * as CUDA events on multiple streams are untrustworthy.
     */
    bool                                             useTiming;

    std::vector<std::unique_ptr<GpuParallel3dFft > > fftSetup;

    std::array<GpuRegionTimer, gtPME_EVENT_COUNT>    timingEvents;

    std::set<size_t>                                 activeTimers; // indices into timingEvents

    /* GPU arrays element counts (not the arrays sizes in bytes!).
     * They might be larger than the actual meaningful data sizes.
     * These are paired: the actual element count + the maximum element count that can fit in the current allocated memory.
     * These integer pairs are mostly meaningful for the cu_realloc/free_buffered calls.
     * As such, if cu_realloc/free_buffered is refactored, they can be freely changed, too.
     * The only exceptions are realGridSize and complexGridSize which are also used for grid clearing/copying.
     * TODO: these should live in a clean buffered container type, and be refactored in the NB/cudautils as well.
     */
    /*! \brief The kernelParams.atoms.coordinates float element count (actual)*/
    int coordinatesSize;
    /*! \brief The kernelParams.atoms.coordinates float element count (reserved) */
    int coordinatesSizeAlloc;
    /*! \brief The kernelParams.atoms.forces float element count (actual) */
    int forcesSize;
    /*! \brief The kernelParams.atoms.forces float element count (reserved) */
    int forcesSizeAlloc;
    /*! \brief The kernelParams.atoms.gridlineIndices int element count (actual) */
    int gridlineIndicesSize;
    /*! \brief The kernelParams.atoms.gridlineIndices int element count (reserved) */
    int gridlineIndicesSizeAlloc;
    /*! \brief Both the kernelParams.atoms.theta and kernelParams.atoms.dtheta float element count (actual) */
    int splineDataSize;
    /*! \brief Both the kernelParams.atoms.theta and kernelParams.atoms.dtheta float element count (reserved) */
    int splineDataSizeAlloc;
    /*! \brief The kernelParams.atoms.coefficients float element count (actual) */
    int coefficientsSize;
    /*! \brief The kernelParams.atoms.coefficients float element count (reserved) */
    int coefficientsSizeAlloc;
    /*! \brief The kernelParams.grid.splineValuesArray float element count (actual) */
    int splineValuesSize;
    /*! \brief The kernelParams.grid.splineValuesArray float element count (reserved) */
    int splineValuesSizeAlloc;
    /*! \brief The kernelParams.grid.realGrid float element count (actual) */
    int realGridSize;
    /*! \brief The kernelParams.grid.realGrid float element count (reserved) */
    int realGridSizeAlloc;
    /*! \brief The kernelParams.grid.fourierGrid float (not float2!) element count (actual) */
    int complexGridSize;
    /*! \brief The kernelParams.grid.fourierGrid float (not float2!) element count (reserved) */
    int complexGridSizeAlloc;
};


/*! \brief \internal
 * A single structure encompassing all the PME data used in CUDA kernels.
 * This inherits from pme_gpu_kernel_params_base_t and adds a couple cudaTextureObject_t handles,
 * which we would like to avoid in plain C++.
 */
struct pme_gpu_cuda_kernel_params_t : pme_gpu_kernel_params_base_t
{
    /* These are CUDA texture objects, related to the grid size. */
    /*! \brief CUDA texture object for accessing grid.d_fractShiftsTable */
    cudaTextureObject_t fractShiftsTableTexture;
    /*! \brief CUDA texture object for accessing grid.d_gridlineIndicesTable */
    cudaTextureObject_t gridlineIndicesTableTexture;
};

/* CUDA texture reference functions which reside in respective kernel files
 * (due to texture references having scope of a translation unit).
 */
/*! Returns the reference to the gridlineIndices texture. */
const struct texture<int, 1, cudaReadModeElementType>   &pme_gpu_get_gridline_texref();
/*! Returns the reference to the fractShifts texture. */
const struct texture<float, 1, cudaReadModeElementType> &pme_gpu_get_fract_shifts_texref();

#endif
