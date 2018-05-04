/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2017,2018, by the GROMACS development team, led by
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
 *
 * \brief This file contains internal function definitions for performing the PME calculations on GPU.
 * These are not meant to be exposed outside of the PME GPU code.
 * As of now, their bodies are still in the common pme-gpu.cpp files.
 *
 * \author Aleksei Iupinov <a.yupinov@gmail.com>
 * \ingroup module_ewald
 */

#ifndef GMX_EWALD_PME_GPU_INTERNAL_H
#define GMX_EWALD_PME_GPU_INTERNAL_H

#include "pme-ocl-definitely-common.h"

#include "gromacs/fft/fft.h"                   // for the gmx_fft_direction enum
#include "gromacs/gpu_utils/gpu_macros.h"      // for the GPU_FUNC_ macros
#include "gromacs/utility/arrayref.h"

#include "pme-gpu-internal-real.h"
#include "pme-gpu-types.h"                     // for the inline functions accessing PmeGpu members

struct gmx_hw_info_t;
struct gmx_gpu_opt_t;
struct gmx_pme_t;                              // only used in pme_gpu_reinit
struct gmx_wallclock_gpu_pme_t;
struct pme_atomcomm_t;
struct t_complex;

namespace gmx
{
class MDLogger;
}

/* Some general constants for PME GPU behaviour follow. */

/*! \brief \libinternal
 * false: The atom data GPU buffers are sized precisely according to the number of atoms.
 *        (Except GPU spline data layout which is regardless intertwined for 2 atoms per warp).
 *        The atom index checks in the spread/gather code potentially hinder the performance.
 * true:  The atom data GPU buffers are padded with zeroes so that the possible number of atoms
 *        fitting in is divisible by PME_ATOM_DATA_ALIGNMENT.
 *        The atom index checks are not performed. There should be a performance win, but how big is it, remains to be seen.
 *        Additional cudaMemsetAsync calls are done occasionally (only charges/coordinates; spline data is always recalculated now).
 * \todo Estimate performance differences
 */
const bool c_usePadding = true;

/*! \brief \libinternal
 * false: Atoms with zero charges are processed by PME. Could introduce some overhead.
 * true:  Atoms with zero charges are not processed by PME. Adds branching to the spread/gather.
 *        Could be good for performance in specific systems with lots of neutral atoms.
 * \todo Estimate performance differences.
 */
const bool c_skipNeutralAtoms = false;

/*! \brief \libinternal
 * Number of PME solve output floating point numbers.
 * 6 for symmetric virial matrix + 1 for reciprocal energy.
 */
//const int c_virialAndEnergyCount = 7;

//FIXME this is a copy!
//! Spreading max block width in warps picked among powers of 2 (2, 4, 8, 16) for max. occupancy and min. runtime in most cases
//FIXME thsi  should account for max block size (e.g. only 4 with OpenCL "warps" of 64  when max work group size is 256)
constexpr int c_spreadMaxWarpsPerBlock = 8;
/* TODO: it has been observed that the kernel can be faster with smaller block sizes (2 or 4 warps)
 * only on some GPUs (660Ti) with large enough grid (>= 48^3) due to load/store units being overloaded
 * (ldst_fu_utilization metric maxed out in nvprof). Runtime block size choice might be nice to have.
 * This has been tried on architectures up to Maxwell (GTX 750) and it would be good to revisit this.
 */

//! Gathering max block width in warps - picked empirically among 2, 4, 8, 16 for max. occupancy and min. runtime
constexpr int c_gatherMaxWarpsPerBlock = 4;
//! Gathering min blocks per CUDA multiprocessor - for CC2.x, we just take the CUDA limit of 8 to avoid the warning
//constexpr int c_gatherMinBlocksPerMP = (GMX_PTX_ARCH < 300) ? GMX_CUDA_MAX_BLOCKS_PER_MP : (GMX_CUDA_MAX_THREADS_PER_MP / c_gatherMaxThreadsPerBlock);
//FIXME dont forget the launch bounds

/* A block of GPU-only functions that live in pme.cu/pme-ocl.cpp */

/*! \libinternal \brief
 * Returns the number of atoms per chunk in the atom charges/coordinates data layout.
 * Depends on CUDA-specific block sizes, needed for the atom data padding.
 *
 * \param[in] pmeGpu            The PME GPU structure.
 * \returns   Number of atoms in a single GPU atom data chunk.
 */
GPU_FUNC_QUALIFIER int pme_gpu_get_atom_data_alignment(const PmeGpu *GPU_FUNC_ARGUMENT(pmeGpu)) GPU_FUNC_TERM_WITH_RETURN(1)

/*! \libinternal \brief
 * Returns the number of atoms per chunk in the atom spline theta/dtheta data layout.
 *
 * \param[in] pmeGpu            The PME GPU structure.
 * \returns   Number of atoms in a single GPU atom spline data chunk.
 */
GPU_FUNC_QUALIFIER int pme_gpu_get_atoms_per_warp(const PmeGpu *GPU_FUNC_ARGUMENT(pmeGpu)) GPU_FUNC_TERM_WITH_RETURN(1)

/*! \libinternal \brief
 * Synchronizes the current computation, waiting for the GPU kernels/transfers to finish.
 *
 * \param[in] pmeGpu            The PME GPU structure.
 */
GPU_FUNC_QUALIFIER void pme_gpu_synchronize(const PmeGpu *GPU_FUNC_ARGUMENT(pmeGpu)) GPU_FUNC_TERM

/*! \libinternal \brief
 * Allocates the fixed size energy and virial buffer both on GPU and CPU.
 *
 * \param[in] pmeGpu            The PME GPU structure.
 */
GPU_FUNC_QUALIFIER void pme_gpu_alloc_energy_virial(const PmeGpu *GPU_FUNC_ARGUMENT(pmeGpu)) GPU_FUNC_TERM

/*! \libinternal \brief
 * Frees the energy and virial memory both on GPU and CPU.
 *
 * \param[in] pmeGpu            The PME GPU structure.
 */
GPU_FUNC_QUALIFIER void pme_gpu_free_energy_virial(PmeGpu *GPU_FUNC_ARGUMENT(pmeGpu)) GPU_FUNC_TERM

/*! \libinternal \brief
 * Clears the energy and virial memory on GPU with 0.
 * Should be called at the end of PME computation which returned energy/virial.
 *
 * \param[in] pmeGpu            The PME GPU structure.
 */
GPU_FUNC_QUALIFIER void pme_gpu_clear_energy_virial(const PmeGpu *GPU_FUNC_ARGUMENT(pmeGpu)) GPU_FUNC_TERM

/*! \libinternal \brief
 * Reallocates and copies the pre-computed B-spline values to the GPU.
 *
 * \param[in] pmeGpu             The PME GPU structure.
 */
GPU_FUNC_QUALIFIER void pme_gpu_realloc_and_copy_bspline_values(const PmeGpu *GPU_FUNC_ARGUMENT(pmeGpu)) GPU_FUNC_TERM

/*! \libinternal \brief
 * Frees the pre-computed B-spline values on the GPU (and the transfer CPU buffers).
 *
 * \param[in] pmeGpu             The PME GPU structure.
 */
GPU_FUNC_QUALIFIER void pme_gpu_free_bspline_values(const PmeGpu *GPU_FUNC_ARGUMENT(pmeGpu)) GPU_FUNC_TERM

/*! \libinternal \brief
 * Reallocates the GPU buffer for the PME forces.
 *
 * \param[in] pmeGpu             The PME GPU structure.
 */
GPU_FUNC_QUALIFIER void pme_gpu_realloc_forces(PmeGpu *GPU_FUNC_ARGUMENT(pmeGpu)) GPU_FUNC_TERM

/*! \libinternal \brief
 * Frees the GPU buffer for the PME forces.
 *
 * \param[in] pmeGpu             The PME GPU structure.
 */
GPU_FUNC_QUALIFIER void pme_gpu_free_forces(const PmeGpu *GPU_FUNC_ARGUMENT(pmeGpu)) GPU_FUNC_TERM

/*! \libinternal \brief
 * Copies the forces from the CPU buffer to the GPU (to reduce them with the PME GPU gathered forces).
 * To be called e.g. after the bonded calculations.
 *
 * \param[in] pmeGpu             The PME GPU structure.
 */
GPU_FUNC_QUALIFIER void pme_gpu_copy_input_forces(PmeGpu *GPU_FUNC_ARGUMENT(pmeGpu)) GPU_FUNC_TERM

/*! \libinternal \brief
 * Copies the forces from the GPU to the CPU buffer. To be called after the gathering stage.
 *
 * \param[in] pmeGpu             The PME GPU structure.
 */
GPU_FUNC_QUALIFIER void pme_gpu_copy_output_forces(PmeGpu *GPU_FUNC_ARGUMENT(pmeGpu)) GPU_FUNC_TERM

/*! \libinternal \brief
 * Checks whether work in the PME GPU stream has completed.
 *
 * \param[in] pmeGpu            The PME GPU structure.
 *
 * \returns                     True if work in the PME stream has completed.
 */
GPU_FUNC_QUALIFIER bool pme_gpu_stream_query(const PmeGpu *GPU_FUNC_ARGUMENT(pmeGpu)) GPU_FUNC_TERM_WITH_RETURN(0)

/*! \libinternal \brief
 * Reallocates the input coordinates buffer on the GPU (and clears the padded part if needed).
 *
 * \param[in] pmeGpu            The PME GPU structure.
 *
 * Needs to be called on every DD step/in the beginning.
 */
GPU_FUNC_QUALIFIER void pme_gpu_realloc_coordinates(const PmeGpu *GPU_FUNC_ARGUMENT(pmeGpu)) GPU_FUNC_TERM

/*! \libinternal \brief
 * Copies the input coordinates from the CPU buffer onto the GPU.
 *
 * \param[in] pmeGpu            The PME GPU structure.
 * \param[in] h_coordinates     Input coordinates (XYZ rvec array).
 *
 * Needs to be called for every PME computation. The coordinates are then used in the spline calculation.
 */
GPU_FUNC_QUALIFIER void pme_gpu_copy_input_coordinates(const PmeGpu    *GPU_FUNC_ARGUMENT(pmeGpu),
                                                        const rvec      *GPU_FUNC_ARGUMENT(h_coordinates)) GPU_FUNC_TERM

/*! \libinternal \brief
 * Frees the coordinates on the GPU.
 *
 * \param[in] pmeGpu            The PME GPU structure.
 */
GPU_FUNC_QUALIFIER void pme_gpu_free_coordinates(const PmeGpu *GPU_FUNC_ARGUMENT(pmeGpu)) GPU_FUNC_TERM

/*! \libinternal \brief
 * Reallocates the buffer on the GPU and copies the charges/coefficients from the CPU buffer.
 * Clears the padded part if needed.
 *
 * \param[in] pmeGpu            The PME GPU structure.
 * \param[in] h_coefficients    The input atom charges/coefficients.
 *
 * Does not need to be done for every PME computation, only whenever the local charges change.
 * (So, in the beginning of the run, or on DD step).
 */
GPU_FUNC_QUALIFIER void pme_gpu_realloc_and_copy_input_coefficients(const PmeGpu    *GPU_FUNC_ARGUMENT(pmeGpu),
                                                                     const float     *GPU_FUNC_ARGUMENT(h_coefficients)) GPU_FUNC_TERM

/*! \libinternal \brief
 * Frees the charges/coefficients on the GPU.
 *
 * \param[in] pmeGpu             The PME GPU structure.
 */
GPU_FUNC_QUALIFIER void pme_gpu_free_coefficients(const PmeGpu *GPU_FUNC_ARGUMENT(pmeGpu)) GPU_FUNC_TERM

/*! \libinternal \brief
 * Reallocates the buffers on the GPU and the host for the atoms spline data.
 *
 * \param[in] pmeGpu            The PME GPU structure.
 */
GPU_FUNC_QUALIFIER void pme_gpu_realloc_spline_data(const PmeGpu *GPU_FUNC_ARGUMENT(pmeGpu)) GPU_FUNC_TERM

/*! \libinternal \brief
 * Frees the buffers on the GPU for the atoms spline data.
 *
 * \param[in] pmeGpu            The PME GPU structure.
 */
GPU_FUNC_QUALIFIER void pme_gpu_free_spline_data(const PmeGpu *GPU_FUNC_ARGUMENT(pmeGpu)) GPU_FUNC_TERM

/*! \libinternal \brief
 * Reallocates the buffers on the GPU and the host for the particle gridline indices.
 *
 * \param[in] pmeGpu            The PME GPU structure.
 */
GPU_FUNC_QUALIFIER void pme_gpu_realloc_grid_indices(const PmeGpu *GPU_FUNC_ARGUMENT(pmeGpu)) GPU_FUNC_TERM

/*! \libinternal \brief
 * Frees the buffer on the GPU for the particle gridline indices.
 *
 * \param[in] pmeGpu            The PME GPU structure.
 */
GPU_FUNC_QUALIFIER void pme_gpu_free_grid_indices(const PmeGpu *GPU_FUNC_ARGUMENT(pmeGpu)) GPU_FUNC_TERM

/*! \libinternal \brief
 * Reallocates the real space grid and the complex reciprocal grid (if needed) on the GPU.
 *
 * \param[in] pmeGpu            The PME GPU structure.
 */
GPU_FUNC_QUALIFIER void pme_gpu_realloc_grids(PmeGpu *GPU_FUNC_ARGUMENT(pmeGpu)) GPU_FUNC_TERM

/*! \libinternal \brief
 * Frees the real space grid and the complex reciprocal grid (if needed) on the GPU.
 *
 * \param[in] pmeGpu            The PME GPU structure.
 */
GPU_FUNC_QUALIFIER void pme_gpu_free_grids(const PmeGpu *GPU_FUNC_ARGUMENT(pmeGpu)) GPU_FUNC_TERM

/*! \libinternal \brief
 * Clears the real space grid on the GPU.
 * Should be called at the end of each computation.
 *
 * \param[in] pmeGpu            The PME GPU structure.
 */
GPU_FUNC_QUALIFIER void pme_gpu_clear_grids(const PmeGpu *GPU_FUNC_ARGUMENT(pmeGpu)) GPU_FUNC_TERM

/*! \libinternal \brief
 * Reallocates and copies the pre-computed fractional coordinates' shifts to the GPU.
 *
 * \param[in] pmeGpu            The PME GPU structure.
 */
GPU_FUNC_QUALIFIER void pme_gpu_realloc_and_copy_fract_shifts(PmeGpu *GPU_FUNC_ARGUMENT(pmeGpu)) GPU_FUNC_TERM

/*! \libinternal \brief
 * Frees the pre-computed fractional coordinates' shifts on the GPU.
 *
 * \param[in] pmeGpu            The PME GPU structure.
 */
GPU_FUNC_QUALIFIER void pme_gpu_free_fract_shifts(const PmeGpu *GPU_FUNC_ARGUMENT(pmeGpu)) GPU_FUNC_TERM

/*! \libinternal \brief
 * Copies the input real-space grid from the host to the GPU.
 *
 * \param[in] pmeGpu   The PME GPU structure.
 * \param[in] h_grid   The host-side grid buffer.
 */
GPU_FUNC_QUALIFIER void pme_gpu_copy_input_gather_grid(const PmeGpu    *GPU_FUNC_ARGUMENT(pmeGpu),
                                                        float           *GPU_FUNC_ARGUMENT(h_grid)) GPU_FUNC_TERM

/*! \libinternal \brief
 * Copies the output real-space grid from the GPU to the host.
 *
 * \param[in] pmeGpu   The PME GPU structure.
 * \param[out] h_grid  The host-side grid buffer.
 */
GPU_FUNC_QUALIFIER void pme_gpu_copy_output_spread_grid(const PmeGpu    *GPU_FUNC_ARGUMENT(pmeGpu),
                                                         float           *GPU_FUNC_ARGUMENT(h_grid)) GPU_FUNC_TERM

/*! \libinternal \brief
 * Copies the spread output spline data and gridline indices from the GPU to the host.
 *
 * \param[in] pmeGpu   The PME GPU structure.
 */
GPU_FUNC_QUALIFIER void pme_gpu_copy_output_spread_atom_data(PmeGpu *pmeGpu) GPU_FUNC_TERM

/*! \libinternal \brief
 * Copies the gather input spline data and gridline indices from the host to the GPU.
 *
 * \param[in] pmeGpu   The PME GPU structure.
 */
GPU_FUNC_QUALIFIER void pme_gpu_copy_input_gather_atom_data(const PmeGpu *GPU_FUNC_ARGUMENT(pmeGpu)) GPU_FUNC_TERM

/*! \libinternal \brief
 * Waits for the grid copying to the host-side buffer after spreading to finish.
 *
 * \param[in] pmeGpu  The PME GPU structure.
 */
GPU_FUNC_QUALIFIER void pme_gpu_sync_spread_grid(const PmeGpu *GPU_FUNC_ARGUMENT(pmeGpu)) GPU_FUNC_TERM

/*! \libinternal \brief
 * Does the one-time GPU-framework specific PME initialization.
 * For CUDA, the PME stream is created with the highest priority.
 *
 * \param[in] pmeGpu  The PME GPU structure.
 */
GPU_FUNC_QUALIFIER void pme_gpu_init_internal(PmeGpu *GPU_FUNC_ARGUMENT(pmeGpu),
					      PmePersistentDataHandle persistent) GPU_FUNC_TERM

/*! \libinternal \brief
 * Destroys the PME GPU-framework specific data.
 * Should be called last in the PME GPU destructor.
 *
 * \param[in] pmeGpu  The PME GPU structure.
 */
GPU_FUNC_QUALIFIER void pme_gpu_destroy_specific(const PmeGpu *GPU_FUNC_ARGUMENT(pmeGpu)) GPU_FUNC_TERM

/*! \libinternal \brief
 * Initializes the GPU FFT structures.
 *
 * \param[in] pmeGpu  The PME GPU structure.
 */
GPU_FUNC_QUALIFIER void pme_gpu_reinit_3dfft(const PmeGpu *GPU_FUNC_ARGUMENT(pmeGpu)) GPU_FUNC_TERM

/*! \libinternal \brief
 * Destroys the GPU FFT structures.
 *
 * \param[in] pmeGpu  The PME GPU structure.
 */
GPU_FUNC_QUALIFIER void pme_gpu_destroy_3dfft(const PmeGpu *GPU_FUNC_ARGUMENT(pmeGpu)) GPU_FUNC_TERM

/* Several CUDA event-based timing functions that live in pme-timings.cu */

/*! \libinternal \brief
 * Finalizes all the active PME GPU stage timings for the current computation. Should be called at the end of every computation.
 *
 * \param[in] pmeGpu         The PME GPU structure.
 */
GPU_FUNC_QUALIFIER void pme_gpu_update_timings(const PmeGpu *GPU_FUNC_ARGUMENT(pmeGpu)) GPU_FUNC_TERM

/*! \libinternal \brief
 * Updates the internal list of active PME GPU stages (if timings are enabled).
 *
 * \param[in] pmeGpu         The PME GPU data structure.
 */
GPU_FUNC_QUALIFIER void pme_gpu_reinit_timings(const PmeGpu *GPU_FUNC_ARGUMENT(pmeGpu)) GPU_FUNC_TERM

/*! \brief
 * Resets the PME GPU timings. To be called at the reset MD step.
 *
 * \param[in] pmeGpu         The PME GPU structure.
 */
GPU_FUNC_QUALIFIER void pme_gpu_reset_timings(const PmeGpu *GPU_FUNC_ARGUMENT(pmeGpu)) GPU_FUNC_TERM

/*! \libinternal \brief
 * Copies the PME GPU timings to the gmx_wallclock_gpu_t structure (for log output). To be called at the run end.
 *
 * \param[in] pmeGpu         The PME GPU structure.
 * \param[in] timings        The gmx_wallclock_gpu_pme_t structure.
 */
GPU_FUNC_QUALIFIER void pme_gpu_get_timings(const PmeGpu            *GPU_FUNC_ARGUMENT(pmeGpu),
                                             gmx_wallclock_gpu_pme_t *GPU_FUNC_ARGUMENT(timings)) GPU_FUNC_TERM

/* The PME stages themselves */

/*! \libinternal \brief
 * A GPU spline computation and charge spreading function.
 *
 * \param[in]  pmeGpu          The PME GPU structure.
 * \param[in]  gridIndex       Index of the PME grid - unused, assumed to be 0.
 * \param[out] h_grid          The host-side grid buffer (used only if the result of the spread is expected on the host,
 *                             e.g. testing or host-side FFT)
 * \param[in]  computeSplines  Should the computation of spline parameters and gridline indices be performed.
 * \param[in]  spreadCharges   Should the charges/coefficients be spread on the grid.
 */
GPU_FUNC_QUALIFIER void pme_gpu_spread(PmeGpu    *GPU_FUNC_ARGUMENT(pmeGpu),
                                        int              GPU_FUNC_ARGUMENT(gridIndex),
                                        real            *GPU_FUNC_ARGUMENT(h_grid),
                                        bool             GPU_FUNC_ARGUMENT(computeSplines),
                                        bool             GPU_FUNC_ARGUMENT(spreadCharges)) GPU_FUNC_TERM

/*! \libinternal \brief
 * 3D FFT R2C/C2R routine.
 *
 * \param[in]  pmeGpu          The PME GPU structure.
 * \param[in]  direction       Transform direction (real-to-complex or complex-to-real)
 * \param[in]  gridIndex       Index of the PME grid - unused, assumed to be 0.
 */
void pme_gpu_3dfft(const PmeGpu *pmeGpu,
                   enum gmx_fft_direction direction,
                   const int gridIndex);

/* The inlined convenience PME GPU status getters */

/*! \libinternal \brief
 * Tells if PME runs on multiple GPUs with the decomposition.
 *
 * \param[in] pmeGpu         The PME GPU structure.
 * \returns                  True if PME runs on multiple GPUs, false otherwise.
 */
inline bool pme_gpu_uses_dd(const PmeGpu *pmeGpu)
{
    return !pmeGpu->settings.useDecomposition;
}

/*! \libinternal \brief
 * Tells if PME performs the gathering stage on GPU.
 *
 * \param[in] pmeGpu         The PME GPU structure.
 * \returns                  True if the gathering is performed on GPU, false otherwise.
 */
inline bool pme_gpu_performs_gather(const PmeGpu *pmeGpu)
{
    return pmeGpu->settings.performGPUGather;
}

/*! \libinternal \brief
 * Tells if PME performs the FFT stages on GPU.
 *
 * \param[in] pmeGpu         The PME GPU structure.
 * \returns                  True if FFT is performed on GPU, false otherwise.
 */
inline bool pme_gpu_performs_FFT(const PmeGpu *pmeGpu)
{
    return pmeGpu->settings.performGPUFFT;
}

/*! \libinternal \brief
 * Tells if PME performs the grid (un-)wrapping on GPU.
 *
 * \param[in] pmeGpu         The PME GPU structure.
 * \returns                  True if (un-)wrapping is performed on GPU, false otherwise.
 */
inline bool pme_gpu_performs_wrapping(const PmeGpu *pmeGpu)
{
    return pmeGpu->settings.useDecomposition;
}

/*! \libinternal \brief
 * Tells if PME performs the grid solving on GPU.
 *
 * \param[in] pmeGpu         The PME GPU structure.
 * \returns                  True if solving is performed on GPU, false otherwise.
 */
inline bool pme_gpu_performs_solve(const PmeGpu *pmeGpu)
{
    return pmeGpu->settings.performGPUSolve;
}

/*! \libinternal \brief
 * Enables or disables the testing mode.
 * Testing mode only implies copying all the outputs, even the intermediate ones, to the host,
 * and also makes the copies synchronous.
 *
 * \param[in] pmeGpu             The PME GPU structure.
 * \param[in] testing            Should the testing mode be enabled, or disabled.
 */
inline void pme_gpu_set_testing(PmeGpu *pmeGpu, bool testing)
{
    pmeGpu->settings.copyAllOutputs = testing;
    pmeGpu->settings.transferKind   = testing ? GpuApiCallBehavior::Sync : GpuApiCallBehavior::Async;
}

/*! \libinternal \brief
 * Tells if PME is in the testing mode.
 *
 * \param[in] pmeGpu             The PME GPU structure.
 * \returns                      true if testing mode is enabled, false otherwise.
 */
inline bool pme_gpu_is_testing(const PmeGpu *pmeGpu)
{
    return pmeGpu->settings.copyAllOutputs;
}

/* A block of C++ functions that live in pme-gpu-internal.cpp */




//FIXME


#include "config.h"

#include <list>
#include <string>

#include "gromacs/ewald/ewald-utils.h"
#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/math/invertmatrix.h"
#include "gromacs/math/units.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/stringutil.h"

#include "pme-grid.h"
#include "pme-internal.h"

/*! \internal \brief
 * Wrapper for getting a pointer to the plain C++ part of the GPU kernel parameters structure.
 *
 * \param[in] pmeGpu  The PME GPU structure.
 * \returns The pointer to the kernel parameters.
 */
static PmeGpuKernelParamsBase *pme_gpu_get_kernel_params_base_ptr(const PmeGpu *pmeGpu)
{
    // reinterpret_cast is needed because the derived CUDA structure is not known in this file
    auto *kernelParamsPtr = reinterpret_cast<PmeGpuKernelParamsBase *>(pmeGpu->kernelParams.get());
    return kernelParamsPtr;
}

inline gmx::ArrayRef<gmx::RVec> pme_gpu_get_forces(PmeGpu *pmeGpu)
{
    return pmeGpu->staging.h_forces;
}

inline void pme_gpu_get_energy_virial(const PmeGpu *pmeGpu, real *energy, matrix virial)
{
    for (int j = 0; j < c_virialAndEnergyCount; j++)
    {
        GMX_ASSERT(std::isfinite(pmeGpu->staging.h_virialAndEnergy[j]), "PME GPU produces incorrect energy/virial.");
    }

    GMX_ASSERT(energy, "Invalid energy output pointer in PME GPU");
    unsigned int j = 0;
    virial[XX][XX] = 0.25f * pmeGpu->staging.h_virialAndEnergy[j++];
    virial[YY][YY] = 0.25f * pmeGpu->staging.h_virialAndEnergy[j++];
    virial[ZZ][ZZ] = 0.25f * pmeGpu->staging.h_virialAndEnergy[j++];
    virial[XX][YY] = virial[YY][XX] = 0.25f * pmeGpu->staging.h_virialAndEnergy[j++];
    virial[XX][ZZ] = virial[ZZ][XX] = 0.25f * pmeGpu->staging.h_virialAndEnergy[j++];
    virial[YY][ZZ] = virial[ZZ][YY] = 0.25f * pmeGpu->staging.h_virialAndEnergy[j++];
    *energy        = 0.5f * pmeGpu->staging.h_virialAndEnergy[j++];
}

inline void pme_gpu_update_input_box(PmeGpu gmx_unused       *pmeGpu,
                              const matrix gmx_unused  box)
{
#if GMX_DOUBLE
    GMX_THROW(gmx::NotImplementedError("PME is implemented for single-precision only on GPU"));
#else
    matrix  scaledBox;
    pmeGpu->common->boxScaler->scaleBox(box, scaledBox);
    auto   *kernelParamsPtr      = pme_gpu_get_kernel_params_base_ptr(pmeGpu);
    kernelParamsPtr->current.boxVolume = scaledBox[XX][XX] * scaledBox[YY][YY] * scaledBox[ZZ][ZZ];
    GMX_ASSERT(kernelParamsPtr->current.boxVolume != 0.0f, "Zero volume of the unit cell");
    matrix recipBox;
    gmx::invertBoxMatrix(scaledBox, recipBox);

    /* The GPU recipBox is transposed as compared to the CPU recipBox.
     * Spread uses matrix columns (while solve and gather use rows).
     * There is no particular reason for this; it might be further rethought/optimized for better access patterns.
     */
    const real newRecipBox[DIM][DIM] =
    {
        {recipBox[XX][XX], recipBox[YY][XX], recipBox[ZZ][XX]},
        {             0.0, recipBox[YY][YY], recipBox[ZZ][YY]},
        {             0.0,              0.0, recipBox[ZZ][ZZ]}
    };
    memcpy(kernelParamsPtr->current.recipBox, newRecipBox, sizeof(matrix));
#endif
}

/*! \brief \libinternal
 * The PME GPU reinitialization function that is called both at the end of any PME computation and on any load balancing.
 *
 * \param[in] pmeGpu            The PME GPU structure.
 */
inline void pme_gpu_reinit_computation(const PmeGpu *pmeGpu)
{
    pme_gpu_clear_grids(pmeGpu);
    pme_gpu_clear_energy_virial(pmeGpu);
}

/*! \brief \libinternal
 * (Re-)initializes all the PME GPU data related to the grid size and cut-off.
 *
 * \param[in] pmeGpu            The PME GPU structure.
 */
static void pme_gpu_reinit_grids(PmeGpu *pmeGpu)
{
    auto *kernelParamsPtr = pme_gpu_get_kernel_params_base_ptr(pmeGpu);
    kernelParamsPtr->grid.ewaldFactor = (M_PI * M_PI) / (pmeGpu->common->ewaldcoeff_q * pmeGpu->common->ewaldcoeff_q);

    /* The grid size variants */
    for (int i = 0; i < DIM; i++)
    {
        kernelParamsPtr->grid.realGridSize[i]       = pmeGpu->common->nk[i];
        kernelParamsPtr->grid.realGridSizeFP[i]     = (float)kernelParamsPtr->grid.realGridSize[i];
        kernelParamsPtr->grid.realGridSizePadded[i] = kernelParamsPtr->grid.realGridSize[i];

        // The complex grid currently uses no padding;
        // if it starts to do so, then another test should be added for that
        kernelParamsPtr->grid.complexGridSize[i]       = kernelParamsPtr->grid.realGridSize[i];
        kernelParamsPtr->grid.complexGridSizePadded[i] = kernelParamsPtr->grid.realGridSize[i];
    }
    /* FFT: n real elements correspond to (n / 2 + 1) complex elements in minor dimension */
    if (!pme_gpu_performs_FFT(pmeGpu))
    {
        // This allows for GPU spreading grid and CPU fftgrid to have the same layout, so that we can copy the data directly
        kernelParamsPtr->grid.realGridSizePadded[ZZ] = (kernelParamsPtr->grid.realGridSize[ZZ] / 2 + 1) * 2;
    }

    /* GPU FFT: n real elements correspond to (n / 2 + 1) complex elements in minor dimension */
    kernelParamsPtr->grid.complexGridSize[ZZ] /= 2;
    kernelParamsPtr->grid.complexGridSize[ZZ]++;
    kernelParamsPtr->grid.complexGridSizePadded[ZZ] = kernelParamsPtr->grid.complexGridSize[ZZ];

    pme_gpu_realloc_and_copy_fract_shifts(pmeGpu);
    pme_gpu_realloc_and_copy_bspline_values(pmeGpu);
    pme_gpu_realloc_grids(pmeGpu);
    pme_gpu_reinit_3dfft(pmeGpu);
}

/* Several GPU functions that refer to the CPU PME data live here.
 * We would like to keep these away from the GPU-framework specific code for clarity,
 * as well as compilation issues with MPI.
 */

/*! \brief \libinternal
 * Copies everything useful from the PME CPU to the PME GPU structure.
 * The goal is to minimize interaction with the PME CPU structure in the GPU code.
 *
 * \param[in] pme         The PME structure.
 */
static void pme_gpu_copy_common_data_from(const gmx_pme_t *pme)
{
    /* TODO: Consider refactoring the CPU PME code to use the same structure,
     * so that this function becomes 2 lines */
    PmeGpu *pmeGpu             = pme->gpu;
    pmeGpu->common->ngrids        = pme->ngrids;
    pmeGpu->common->epsilon_r     = pme->epsilon_r;
    pmeGpu->common->ewaldcoeff_q  = pme->ewaldcoeff_q;
    pmeGpu->common->nk[XX]        = pme->nkx;
    pmeGpu->common->nk[YY]        = pme->nky;
    pmeGpu->common->nk[ZZ]        = pme->nkz;
    pmeGpu->common->pme_order     = pme->pme_order;
    for (int i = 0; i < DIM; i++)
    {
        pmeGpu->common->bsp_mod[i].assign(pme->bsp_mod[i], pme->bsp_mod[i] + pmeGpu->common->nk[i]);
    }
    const int cellCount = c_pmeNeighborUnitcellCount;
    pmeGpu->common->fsh.resize(0);
    pmeGpu->common->fsh.insert(pmeGpu->common->fsh.end(), pme->fshx, pme->fshx + cellCount * pme->nkx);
    pmeGpu->common->fsh.insert(pmeGpu->common->fsh.end(), pme->fshy, pme->fshy + cellCount * pme->nky);
    pmeGpu->common->fsh.insert(pmeGpu->common->fsh.end(), pme->fshz, pme->fshz + cellCount * pme->nkz);
    pmeGpu->common->nn.resize(0);
    pmeGpu->common->nn.insert(pmeGpu->common->nn.end(), pme->nnx, pme->nnx + cellCount * pme->nkx);
    pmeGpu->common->nn.insert(pmeGpu->common->nn.end(), pme->nny, pme->nny + cellCount * pme->nky);
    pmeGpu->common->nn.insert(pmeGpu->common->nn.end(), pme->nnz, pme->nnz + cellCount * pme->nkz);
    pmeGpu->common->runMode   = pme->runMode;
    pmeGpu->common->boxScaler = pme->boxScaler;
}

/*! \brief \libinternal
 * Finds out if PME with given inputs is possible to run on GPU.
 *
 * \param[in]  pme          The PME structure.
 * \param[out] error        The error message if the input is not supported on GPU.
 * \returns                 True if this PME input is possible to run on GPU, false otherwise.
 */
static bool pme_gpu_check_restrictions(const gmx_pme_t *pme, std::string *error)
{
    std::list<std::string> errorReasons;
    if (pme->nnodes != 1)
    {
        errorReasons.push_back("PME decomposition");
    }
    if (pme->pme_order != 4)
    {
        errorReasons.push_back("interpolation orders other than 4");
    }
    if (pme->bFEP)
    {
        errorReasons.push_back("free energy calculations (multiple grids)");
    }
    if (pme->doLJ)
    {
        errorReasons.push_back("Lennard-Jones PME");
    }
#if GMX_DOUBLE
    {
        errorReasons.push_back("double precision");
    }
#endif
#if GMX_GPU == GMX_GPU_NONE
    {
        errorReasons.push_back("non-GPU build of GROMACS");
    }
#endif

    bool inputSupported = errorReasons.empty();
    if (!inputSupported && error)
    {
        std::string regressionTestMarker = "PME GPU does not support";
        // this prefix is tested for in the regression tests script gmxtest.pl
        *error = regressionTestMarker + ": " + gmx::joinStrings(errorReasons, "; ") + ".";
    }
    return inputSupported;
}

/*! \libinternal \brief
 * Initializes the PME GPU data at the beginning of the run.
 *
 * \param[in,out] pme       The PME structure.
 * \param[in,out] gpuInfo   The GPU information structure.
 */
static void pme_gpu_init(gmx_pme_t *pme, gmx_device_info_t *gpuInfo, PmePersistentDataHandle persistent)
{
    std::string errorString;
    bool        canRunOnGpu = pme_gpu_check_restrictions(pme, &errorString);
    if (!canRunOnGpu)
    {
        GMX_THROW(gmx::NotImplementedError(errorString));
    }

    pme->gpu          = new PmeGpu();
    PmeGpu *pmeGpu = pme->gpu;
    changePinningPolicy(&pmeGpu->staging.h_forces, gmx::PinningPolicy::CanBePinned);
    pmeGpu->common = std::shared_ptr<PmeShared>(new PmeShared());

    /* These settings are set here for the whole run; dynamic ones are set in pme_gpu_reinit() */
    /* A convenience variable. */
    pmeGpu->settings.useDecomposition = (pme->nnodes == 1);
    /* TODO: CPU gather with GPU spread is broken due to different theta/dtheta layout. */
    pmeGpu->settings.performGPUGather = true;

    pme_gpu_set_testing(pmeGpu, false);

    pmeGpu->deviceInfo = gpuInfo;

    pme_gpu_init_internal(pmeGpu, persistent);
    pme_gpu_alloc_energy_virial(pmeGpu);

    pme_gpu_copy_common_data_from(pme);

    GMX_ASSERT(pmeGpu->common->epsilon_r != 0.0f, "PME GPU: bad electrostatic coefficient");

    auto *kernelParamsPtr = pme_gpu_get_kernel_params_base_ptr(pmeGpu);
    kernelParamsPtr->constants.elFactor = ONE_4PI_EPS0 / pmeGpu->common->epsilon_r;
}
#if 0 //moved to pme-ocl.cpp
void pme_gpu_transform_spline_atom_data(const PmeGpu *pmeGpu, const pme_atomcomm_t *atc,
                                        PmeSplineDataType type, int dimIndex, PmeLayoutTransform transform)
{
    // The GPU atom spline data is laid out in a different way currently than the CPU one.
    // This function converts the data from GPU to CPU layout (in the host memory).
    // It is only intended for testing purposes so far.
    // Ideally we should use similar layouts on CPU and GPU if we care about mixed modes and their performance
    // (e.g. spreading on GPU, gathering on CPU).
    GMX_RELEASE_ASSERT(atc->nthread == 1, "Only the serial PME data layout is supported");
    const uintmax_t threadIndex  = 0;
    const auto      atomCount    = pme_gpu_get_kernel_params_base_ptr(pmeGpu)->atoms.nAtoms;
    const auto      atomsPerWarp = pme_gpu_get_atoms_per_warp(pmeGpu);
    const auto      pmeOrder     = pmeGpu->common->pme_order;

    real           *cpuSplineBuffer;
    float          *h_splineBuffer;
    switch (type)
    {
        case PmeSplineDataType::Values:
            cpuSplineBuffer = atc->spline[threadIndex].theta[dimIndex];
            h_splineBuffer  = pmeGpu->staging.h_theta;
            break;

        case PmeSplineDataType::Derivatives:
            cpuSplineBuffer = atc->spline[threadIndex].dtheta[dimIndex];
            h_splineBuffer  = pmeGpu->staging.h_dtheta;
            break;

        default:
            GMX_THROW(gmx::InternalError("Unknown spline data type"));
    }

    for (auto atomIndex = 0; atomIndex < atomCount; atomIndex++)
    {
        auto atomWarpIndex = atomIndex % atomsPerWarp;
        auto warpIndex     = atomIndex / atomsPerWarp;
        for (auto orderIndex = 0; orderIndex < pmeOrder; orderIndex++)
        {
            const auto gpuValueIndex = ((pmeOrder * warpIndex + orderIndex) * DIM + dimIndex) * atomsPerWarp + atomWarpIndex;
            const auto cpuValueIndex = atomIndex * pmeOrder + orderIndex;
            GMX_ASSERT(cpuValueIndex < atomCount * pmeOrder, "Atom spline data index out of bounds (while transforming GPU data layout for host)");
            switch (transform)
            {
                case PmeLayoutTransform::GpuToHost:
                    cpuSplineBuffer[cpuValueIndex] = h_splineBuffer[gpuValueIndex];
                    break;

                case PmeLayoutTransform::HostToGpu:
                    h_splineBuffer[gpuValueIndex] = cpuSplineBuffer[cpuValueIndex];
                    break;

                default:
                    GMX_THROW(gmx::InternalError("Unknown layout transform"));
            }
        }
    }
}

void pme_gpu_get_real_grid_sizes(const PmeGpu *pmeGpu, gmx::IVec *gridSize, gmx::IVec *paddedGridSize)
{
    GMX_ASSERT(gridSize != nullptr, "");
    GMX_ASSERT(paddedGridSize != nullptr, "");
    GMX_ASSERT(pmeGpu != nullptr, "");
    auto *kernelParamsPtr = pme_gpu_get_kernel_params_base_ptr(pmeGpu);
    for (int i = 0; i < DIM; i++)
    {
        (*gridSize)[i]       = kernelParamsPtr->grid.realGridSize[i];
        (*paddedGridSize)[i] = kernelParamsPtr->grid.realGridSizePadded[i];
    }
}

inline void pme_gpu_reinit(gmx_pme_t *pme, gmx_device_info_t *gpuInfo, PmePersistentDataHandle persistent)
{
    if (!pme_gpu_active(pme))
    {
        return;
    }

    if (!pme->gpu)
    {
        /* First-time initialization */
      pme_gpu_init(pme, gpuInfo, persistent);
    }
    else
    {
        /* After this call nothing in the GPU code should refer to the gmx_pme_t *pme itself - until the next pme_gpu_reinit */
        pme_gpu_copy_common_data_from(pme);
    }
    /* GPU FFT will only get used for a single rank.*/
    pme->gpu->settings.performGPUFFT   = (pme->gpu->common->runMode == PmeRunMode::GPU) && !pme_gpu_uses_dd(pme->gpu);
    pme->gpu->settings.performGPUSolve = (pme->gpu->common->runMode == PmeRunMode::GPU);

    /* Reinit active timers */
    pme_gpu_reinit_timings(pme->gpu);

    pme_gpu_reinit_grids(pme->gpu);
    pme_gpu_reinit_computation(pme->gpu);
    /* Clear the previous box - doesn't hurt, and forces the PME CPU recipbox
     * update for mixed mode on grid switch. TODO: use shared recipbox field.
     */
    std::memset(pme->gpu->common->previousBox, 0, sizeof(pme->gpu->common->previousBox));
}


inline void pme_gpu_destroy(PmeGpu *pmeGpu)
{
    /* Free lots of data */
    pme_gpu_free_energy_virial(pmeGpu);
    pme_gpu_free_bspline_values(pmeGpu);
    pme_gpu_free_forces(pmeGpu);
    pme_gpu_free_coordinates(pmeGpu);
    pme_gpu_free_coefficients(pmeGpu);
    pme_gpu_free_spline_data(pmeGpu);
    pme_gpu_free_grid_indices(pmeGpu);
    pme_gpu_free_fract_shifts(pmeGpu);
    pme_gpu_free_grids(pmeGpu);

    pme_gpu_destroy_3dfft(pmeGpu);

    /* Free the GPU-framework specific data last */
    pme_gpu_destroy_specific(pmeGpu);

    delete pmeGpu;
}


inline void pme_gpu_reinit_atoms(PmeGpu *pmeGpu, const int nAtoms, const real *charges)
{
    auto      *kernelParamsPtr = pme_gpu_get_kernel_params_base_ptr(pmeGpu);
    kernelParamsPtr->atoms.nAtoms = nAtoms;
    const int  alignment = pme_gpu_get_atom_data_alignment(pmeGpu);
    pmeGpu->nAtomsPadded = ((nAtoms + alignment - 1) / alignment) * alignment;
    const int  nAtomsAlloc   = c_usePadding ? pmeGpu->nAtomsPadded : nAtoms;
    const bool haveToRealloc = (pmeGpu->nAtomsAlloc < nAtomsAlloc); /* This check might be redundant, but is logical */
    pmeGpu->nAtomsAlloc = nAtomsAlloc;

#if GMX_DOUBLE
    GMX_RELEASE_ASSERT(false, "Only single precision supported");
    GMX_UNUSED_VALUE(charges);
#else
    pme_gpu_realloc_and_copy_input_coefficients(pmeGpu, reinterpret_cast<const float *>(charges));
    /* Could also be checked for haveToRealloc, but the copy always needs to be performed */
#endif

    if (haveToRealloc)
    {
        pme_gpu_realloc_coordinates(pmeGpu);
        pme_gpu_realloc_forces(pmeGpu);
        pme_gpu_realloc_spline_data(pmeGpu);
        pme_gpu_realloc_grid_indices(pmeGpu);
    }
}
#endif


#if 0
/*! \libinternal \brief
 * Returns the GPU gathering staging forces buffer.
 *
 * \param[in] pmeGpu             The PME GPU structure.
 * \returns                      The input/output forces.
 */
gmx::ArrayRef<gmx::RVec> pme_gpu_get_forces(PmeGpu *pmeGpu);

/*! \libinternal \brief
 * Returns the output virial and energy of the PME solving.
 *
 * \param[in] pmeGpu             The PME GPU structure.
 * \param[out] energy            The output energy.
 * \param[out] virial            The output virial matrix.
 */
void pme_gpu_get_energy_virial(const PmeGpu *pmeGpu, real *energy, matrix virial);

/*! \libinternal \brief
 * Updates the unit cell parameters. Does not check if update is necessary - that is done in pme_gpu_prepare_computation().
 *
 * \param[in] pmeGpu         The PME GPU structure.
 * \param[in] box            The unit cell box.
 */
void pme_gpu_update_input_box(PmeGpu *pmeGpu, const matrix box);

/*! \libinternal \brief
 * Finishes the PME GPU computation, waiting for the output forces and/or energy/virial to be copied to the host.
 * If forces were computed, they will have arrived at the external host buffer provided to gather.
 * If virial/energy were computed, they will have arrived into the internal staging buffer
 * (even though that should have already happened before even launching the gather).
 * Finally, cudaEvent_t based GPU timers get updated if enabled. They also need stream synchronization for correctness.
 * Additionally, device-side buffers are cleared asynchronously for the next computation.
 *
 * \param[in] pmeGpu         The PME GPU structure.
 */
void pme_gpu_finish_computation(const PmeGpu *pmeGpu);

/*! \libinternal \brief
 * Rearranges the atom spline data between the GPU and host layouts.
 * Only used for test purposes so far, likely to be horribly slow.
 *
 * \param[in]  pmeGpu     The PME GPU structure.
 * \param[out] atc        The PME CPU atom data structure (with a single-threaded layout).
 * \param[in]  type       The spline data type (values or derivatives).
 * \param[in]  dimIndex   Dimension index.
 * \param[in]  transform  Layout transform type
 */
void pme_gpu_transform_spline_atom_data(const PmeGpu *pmeGpu, const pme_atomcomm_t *atc,
                                        PmeSplineDataType type, int dimIndex, PmeLayoutTransform transform);

/*! \libinternal \brief
 * Get the normal/padded grid dimensions of the real-space PME grid on GPU. Only used in tests.
 *
 * \param[in] pmeGpu             The PME GPU structure.
 * \param[out] gridSize          Pointer to the grid dimensions to fill in.
 * \param[out] paddedGridSize    Pointer to the padded grid dimensions to fill in.
 */
void pme_gpu_get_real_grid_sizes(const PmeGpu *pmeGpu, gmx::IVec *gridSize, gmx::IVec *paddedGridSize);

/*! \libinternal \brief
 * (Re-)initializes the PME GPU data at the beginning of the run or on DLB.
 *
 * \param[in,out] pme       The PME structure.
 * \param[in,out] gpuInfo   The GPU information structure.
 * \throws gmx::NotImplementedError if this generally valid PME structure is not valid for GPU runs.
 */
void pme_gpu_reinit(gmx_pme_t *pme, gmx_device_info_t *gpuInfo);

/*! \libinternal \brief
 * Destroys the PME GPU data at the end of the run.
 *
 * \param[in] pmeGpu     The PME GPU structure.
 */
void pme_gpu_destroy(PmeGpu *pmeGpu);

/*! \libinternal \brief
 * Reallocates the local atoms data (charges, coordinates, etc.). Copies the charges to the GPU.
 *
 * \param[in] pmeGpu    The PME GPU structure.
 * \param[in] nAtoms    The number of particles.
 * \param[in] charges   The pointer to the host-side array of particle charges.
 *
 * This is a function that should only be called in the beginning of the run and on domain decomposition.
 * Should be called before the pme_gpu_set_io_ranges.
 */
void pme_gpu_reinit_atoms(PmeGpu           *pmeGpu,
                          const int         nAtoms,
                          const real       *charges);

/*! \brief \libinternal
 * The PME GPU reinitialization function that is called both at the end of any PME computation and on any load balancing.
 *
 * This clears the device-side working buffers in preparation for new computation.
 *
 * \param[in] pmeGpu            The PME GPU structure.
 */
void pme_gpu_reinit_computation(const PmeGpu *pmeGpu);
#endif

#endif
