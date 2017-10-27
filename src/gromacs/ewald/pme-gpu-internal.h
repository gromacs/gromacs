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

#include "gromacs/fft/fft.h"                   // for the gmx_fft_direction enum
#include "gromacs/gpu_utils/gpu_macros.h"      // for the CUDA_FUNC_ macros
#include "gromacs/utility/arrayref.h"

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

//! Type of spline data
enum class PmeSplineDataType
{
    Values,      // theta
    Derivatives, // dtheta
};               //TODO move this into new and shiny pme.h (pme-types.h?)

//! PME grid dimension ordering (from major to minor)
enum class GridOrdering
{
    YZX,
    XYZ
};

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
const int c_virialAndEnergyCount = 7;

/* A block of CUDA-only functions that live in pme.cu */

/*! \libinternal \brief
 * Returns the number of atoms per chunk in the atom charges/coordinates data layout.
 * Depends on CUDA-specific block sizes, needed for the atom data padding.
 *
 * \param[in] pmeGpu            The PME GPU structure.
 * \returns   Number of atoms in a single GPU atom data chunk.
 */
CUDA_FUNC_QUALIFIER int pme_gpu_get_atom_data_alignment(const PmeGpu *CUDA_FUNC_ARGUMENT(pmeGpu)) CUDA_FUNC_TERM_WITH_RETURN(1)

/*! \libinternal \brief
 * Returns the number of atoms per chunk in the atom spline theta/dtheta data layout.
 *
 * \param[in] pmeGpu            The PME GPU structure.
 * \returns   Number of atoms in a single GPU atom spline data chunk.
 */
CUDA_FUNC_QUALIFIER int pme_gpu_get_atoms_per_warp(const PmeGpu *CUDA_FUNC_ARGUMENT(pmeGpu)) CUDA_FUNC_TERM_WITH_RETURN(1)

/*! \libinternal \brief
 * Synchronizes the current computation, waiting for the GPU kernels/transfers to finish.
 *
 * \param[in] pmeGpu            The PME GPU structure.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_synchronize(const PmeGpu *CUDA_FUNC_ARGUMENT(pmeGpu)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * Allocates the fixed size energy and virial buffer both on GPU and CPU.
 *
 * \param[in] pmeGpu            The PME GPU structure.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_alloc_energy_virial(const PmeGpu *CUDA_FUNC_ARGUMENT(pmeGpu)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * Frees the energy and virial memory both on GPU and CPU.
 *
 * \param[in] pmeGpu            The PME GPU structure.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_free_energy_virial(PmeGpu *CUDA_FUNC_ARGUMENT(pmeGpu)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * Clears the energy and virial memory on GPU with 0.
 * Should be called at the end of PME computation which returned energy/virial.
 *
 * \param[in] pmeGpu            The PME GPU structure.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_clear_energy_virial(const PmeGpu *CUDA_FUNC_ARGUMENT(pmeGpu)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * Reallocates and copies the pre-computed B-spline values to the GPU.
 *
 * \param[in] pmeGpu             The PME GPU structure.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_realloc_and_copy_bspline_values(const PmeGpu *CUDA_FUNC_ARGUMENT(pmeGpu)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * Frees the pre-computed B-spline values on the GPU (and the transfer CPU buffers).
 *
 * \param[in] pmeGpu             The PME GPU structure.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_free_bspline_values(const PmeGpu *CUDA_FUNC_ARGUMENT(pmeGpu)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * Reallocates the GPU buffer for the PME forces.
 *
 * \param[in] pmeGpu             The PME GPU structure.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_realloc_forces(PmeGpu *CUDA_FUNC_ARGUMENT(pmeGpu)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * Frees the GPU buffer for the PME forces.
 *
 * \param[in] pmeGpu             The PME GPU structure.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_free_forces(const PmeGpu *CUDA_FUNC_ARGUMENT(pmeGpu)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * Copies the forces from the CPU buffer to the GPU (to reduce them with the PME GPU gathered forces).
 * To be called e.g. after the bonded calculations.
 *
 * \param[in] pmeGpu             The PME GPU structure.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_copy_input_forces(PmeGpu *CUDA_FUNC_ARGUMENT(pmeGpu)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * Copies the forces from the GPU to the CPU buffer. To be called after the gathering stage.
 *
 * \param[in] pmeGpu             The PME GPU structure.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_copy_output_forces(PmeGpu *CUDA_FUNC_ARGUMENT(pmeGpu)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * Checks whether work in the PME GPU stream has completed.
 *
 * \param[in] pmeGpu            The PME GPU structure.
 *
 * \returns                     True if work in the PME stream has completed.
 */
CUDA_FUNC_QUALIFIER bool pme_gpu_stream_query(const PmeGpu *CUDA_FUNC_ARGUMENT(pmeGpu)) CUDA_FUNC_TERM_WITH_RETURN(0)

/*! \libinternal \brief
 * Reallocates the input coordinates buffer on the GPU (and clears the padded part if needed).
 *
 * \param[in] pmeGpu            The PME GPU structure.
 *
 * Needs to be called on every DD step/in the beginning.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_realloc_coordinates(const PmeGpu *CUDA_FUNC_ARGUMENT(pmeGpu)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * Copies the input coordinates from the CPU buffer onto the GPU.
 *
 * \param[in] pmeGpu            The PME GPU structure.
 * \param[in] h_coordinates     Input coordinates (XYZ rvec array).
 *
 * Needs to be called for every PME computation. The coordinates are then used in the spline calculation.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_copy_input_coordinates(const PmeGpu    *CUDA_FUNC_ARGUMENT(pmeGpu),
                                                        const rvec      *CUDA_FUNC_ARGUMENT(h_coordinates)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * Frees the coordinates on the GPU.
 *
 * \param[in] pmeGpu            The PME GPU structure.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_free_coordinates(const PmeGpu *CUDA_FUNC_ARGUMENT(pmeGpu)) CUDA_FUNC_TERM

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
CUDA_FUNC_QUALIFIER void pme_gpu_realloc_and_copy_input_coefficients(const PmeGpu    *CUDA_FUNC_ARGUMENT(pmeGpu),
                                                                     const float     *CUDA_FUNC_ARGUMENT(h_coefficients)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * Frees the charges/coefficients on the GPU.
 *
 * \param[in] pmeGpu             The PME GPU structure.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_free_coefficients(const PmeGpu *CUDA_FUNC_ARGUMENT(pmeGpu)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * Reallocates the buffers on the GPU and the host for the atoms spline data.
 *
 * \param[in] pmeGpu            The PME GPU structure.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_realloc_spline_data(const PmeGpu *CUDA_FUNC_ARGUMENT(pmeGpu)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * Frees the buffers on the GPU for the atoms spline data.
 *
 * \param[in] pmeGpu            The PME GPU structure.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_free_spline_data(const PmeGpu *CUDA_FUNC_ARGUMENT(pmeGpu)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * Reallocates the buffers on the GPU and the host for the particle gridline indices.
 *
 * \param[in] pmeGpu            The PME GPU structure.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_realloc_grid_indices(const PmeGpu *CUDA_FUNC_ARGUMENT(pmeGpu)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * Frees the buffer on the GPU for the particle gridline indices.
 *
 * \param[in] pmeGpu            The PME GPU structure.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_free_grid_indices(const PmeGpu *CUDA_FUNC_ARGUMENT(pmeGpu)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * Reallocates the real space grid and the complex reciprocal grid (if needed) on the GPU.
 *
 * \param[in] pmeGpu            The PME GPU structure.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_realloc_grids(PmeGpu *CUDA_FUNC_ARGUMENT(pmeGpu)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * Frees the real space grid and the complex reciprocal grid (if needed) on the GPU.
 *
 * \param[in] pmeGpu            The PME GPU structure.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_free_grids(const PmeGpu *CUDA_FUNC_ARGUMENT(pmeGpu)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * Clears the real space grid on the GPU.
 * Should be called at the end of each computation.
 *
 * \param[in] pmeGpu            The PME GPU structure.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_clear_grids(const PmeGpu *CUDA_FUNC_ARGUMENT(pmeGpu)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * Reallocates and copies the pre-computed fractional coordinates' shifts to the GPU.
 *
 * \param[in] pmeGpu            The PME GPU structure.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_realloc_and_copy_fract_shifts(PmeGpu *CUDA_FUNC_ARGUMENT(pmeGpu)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * Frees the pre-computed fractional coordinates' shifts on the GPU.
 *
 * \param[in] pmeGpu            The PME GPU structure.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_free_fract_shifts(const PmeGpu *CUDA_FUNC_ARGUMENT(pmeGpu)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * Copies the input real-space grid from the host to the GPU.
 *
 * \param[in] pmeGpu   The PME GPU structure.
 * \param[in] h_grid   The host-side grid buffer.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_copy_input_gather_grid(const PmeGpu    *CUDA_FUNC_ARGUMENT(pmeGpu),
                                                        float           *CUDA_FUNC_ARGUMENT(h_grid)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * Copies the output real-space grid from the GPU to the host.
 *
 * \param[in] pmeGpu   The PME GPU structure.
 * \param[out] h_grid  The host-side grid buffer.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_copy_output_spread_grid(const PmeGpu    *CUDA_FUNC_ARGUMENT(pmeGpu),
                                                         float           *CUDA_FUNC_ARGUMENT(h_grid)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * Copies the spread output spline data and gridline indices from the GPU to the host.
 *
 * \param[in] pmeGpu   The PME GPU structure.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_copy_output_spread_atom_data(const PmeGpu *CUDA_FUNC_ARGUMENT(pmeGpu)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * Copies the gather input spline data and gridline indices from the host to the GPU.
 *
 * \param[in] pmeGpu   The PME GPU structure.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_copy_input_gather_atom_data(const PmeGpu *CUDA_FUNC_ARGUMENT(pmeGpu)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * Waits for the grid copying to the host-side buffer after spreading to finish.
 *
 * \param[in] pmeGpu  The PME GPU structure.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_sync_spread_grid(const PmeGpu *CUDA_FUNC_ARGUMENT(pmeGpu)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * Does the one-time GPU-framework specific PME initialization.
 * For CUDA, the PME stream is created with the highest priority.
 *
 * \param[in] pmeGpu  The PME GPU structure.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_init_internal(PmeGpu *CUDA_FUNC_ARGUMENT(pmeGpu)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * Destroys the PME GPU-framework specific data.
 * Should be called last in the PME GPU destructor.
 *
 * \param[in] pmeGpu  The PME GPU structure.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_destroy_specific(const PmeGpu *CUDA_FUNC_ARGUMENT(pmeGpu)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * Initializes the PME GPU synchronization events.
 *
 * \param[in] pmeGpu  The PME GPU structure.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_init_sync_events(const PmeGpu *CUDA_FUNC_ARGUMENT(pmeGpu)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * Destroys the PME GPU synchronization events.
 *
 * \param[in] pmeGpu  The PME GPU structure.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_destroy_sync_events(const PmeGpu *CUDA_FUNC_ARGUMENT(pmeGpu)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * Initializes the CUDA FFT structures.
 *
 * \param[in] pmeGpu  The PME GPU structure.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_reinit_3dfft(const PmeGpu *CUDA_FUNC_ARGUMENT(pmeGpu)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * Destroys the CUDA FFT structures.
 *
 * \param[in] pmeGpu  The PME GPU structure.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_destroy_3dfft(const PmeGpu *CUDA_FUNC_ARGUMENT(pmeGpu)) CUDA_FUNC_TERM

/* Several CUDA event-based timing functions that live in pme-timings.cu */

/*! \libinternal \brief
 * Finalizes all the active PME GPU stage timings for the current computation. Should be called at the end of every computation.
 *
 * \param[in] pmeGpu         The PME GPU structure.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_update_timings(const PmeGpu *CUDA_FUNC_ARGUMENT(pmeGpu)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * Updates the internal list of active PME GPU stages (if timings are enabled).
 *
 * \param[in] pmeGpu         The PME GPU data structure.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_reinit_timings(const PmeGpu *CUDA_FUNC_ARGUMENT(pmeGpu)) CUDA_FUNC_TERM

/*! \brief
 * Resets the PME GPU timings. To be called at the reset MD step.
 *
 * \param[in] pmeGpu         The PME GPU structure.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_reset_timings(const PmeGpu *CUDA_FUNC_ARGUMENT(pmeGpu)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * Copies the PME GPU timings to the gmx_wallclock_gpu_t structure (for log output). To be called at the run end.
 *
 * \param[in] pmeGpu         The PME GPU structure.
 * \param[in] timings        The gmx_wallclock_gpu_pme_t structure.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_get_timings(const PmeGpu            *CUDA_FUNC_ARGUMENT(pmeGpu),
                                             gmx_wallclock_gpu_pme_t *CUDA_FUNC_ARGUMENT(timings)) CUDA_FUNC_TERM

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
CUDA_FUNC_QUALIFIER void pme_gpu_spread(const PmeGpu    *CUDA_FUNC_ARGUMENT(pmeGpu),
                                        int              CUDA_FUNC_ARGUMENT(gridIndex),
                                        real            *CUDA_FUNC_ARGUMENT(h_grid),
                                        bool             CUDA_FUNC_ARGUMENT(computeSplines),
                                        bool             CUDA_FUNC_ARGUMENT(spreadCharges)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * 3D FFT R2C/C2R routine.
 *
 * \param[in]  pmeGpu          The PME GPU structure.
 * \param[in]  direction       Transform direction (real-to-complex or complex-to-real)
 * \param[in]  gridIndex       Index of the PME grid - unused, assumed to be 0.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_3dfft(const PmeGpu          *CUDA_FUNC_ARGUMENT(pmeGpu),
                                       enum gmx_fft_direction CUDA_FUNC_ARGUMENT(direction),
                                       const int              CUDA_FUNC_ARGUMENT(gridIndex)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * A GPU Fourier space solving function.
 *
 * \param[in]     pmeGpu                  The PME GPU structure.
 * \param[in,out] h_grid                  The host-side input and output Fourier grid buffer (used only with testing or host-side FFT)
 * \param[in]     gridOrdering            Specifies the dimenion ordering of the complex grid. TODO: store this information?
 * \param[in]     computeEnergyAndVirial  Tells if the energy and virial computation should also be performed.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_solve(const PmeGpu    *CUDA_FUNC_ARGUMENT(pmeGpu),
                                       t_complex       *CUDA_FUNC_ARGUMENT(h_grid),
                                       GridOrdering     CUDA_FUNC_ARGUMENT(gridOrdering),
                                       bool             CUDA_FUNC_ARGUMENT(computeEnergyAndVirial)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * A GPU force gathering function.
 *
 * \param[in]     pmeGpu           The PME GPU structure.
 * \param[in]     forceTreatment   Tells how data in h_forces should be treated.
 *                                 TODO: determine efficiency/balance of host/device-side reductions.
 * \param[in]     h_grid           The host-side grid buffer (used only in testing mode)
 */
CUDA_FUNC_QUALIFIER void pme_gpu_gather(PmeGpu                *CUDA_FUNC_ARGUMENT(pmeGpu),
                                        PmeForceOutputHandling CUDA_FUNC_ARGUMENT(forceTreatment),
                                        const float           *CUDA_FUNC_ARGUMENT(h_grid)
                                        ) CUDA_FUNC_TERM


/* The inlined convenience PME GPU status getters */

/*! \libinternal \brief
 * Tells if PME runs on multiple GPUs with the decomposition.
 *
 * \param[in] pmeGpu         The PME GPU structure.
 * \returns                  True if PME runs on multiple GPUs, false otherwise.
 */
gmx_inline bool pme_gpu_uses_dd(const PmeGpu *pmeGpu)
{
    return !pmeGpu->settings.useDecomposition;
}

/*! \libinternal \brief
 * Tells if PME performs the gathering stage on GPU.
 *
 * \param[in] pmeGpu         The PME GPU structure.
 * \returns                  True if the gathering is performed on GPU, false otherwise.
 */
gmx_inline bool pme_gpu_performs_gather(const PmeGpu *pmeGpu)
{
    return pmeGpu->settings.performGPUGather;
}

/*! \libinternal \brief
 * Tells if PME performs the FFT stages on GPU.
 *
 * \param[in] pmeGpu         The PME GPU structure.
 * \returns                  True if FFT is performed on GPU, false otherwise.
 */
gmx_inline bool pme_gpu_performs_FFT(const PmeGpu *pmeGpu)
{
    return pmeGpu->settings.performGPUFFT;
}

/*! \libinternal \brief
 * Tells if PME performs the grid (un-)wrapping on GPU.
 *
 * \param[in] pmeGpu         The PME GPU structure.
 * \returns                  True if (un-)wrapping is performed on GPU, false otherwise.
 */
gmx_inline bool pme_gpu_performs_wrapping(const PmeGpu *pmeGpu)
{
    return pmeGpu->settings.useDecomposition;
}

/*! \libinternal \brief
 * Tells if PME performs the grid solving on GPU.
 *
 * \param[in] pmeGpu         The PME GPU structure.
 * \returns                  True if solving is performed on GPU, false otherwise.
 */
gmx_inline bool pme_gpu_performs_solve(const PmeGpu *pmeGpu)
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
gmx_inline void pme_gpu_set_testing(PmeGpu *pmeGpu, bool testing)
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
gmx_inline bool pme_gpu_is_testing(const PmeGpu *pmeGpu)
{
    return pmeGpu->settings.copyAllOutputs;
}

/* A block of C++ functions that live in pme-gpu-internal.cpp */

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

//! A binary enum for spline data layout transformation
enum class PmeLayoutTransform
{
    GpuToHost,
    HostToGpu
};

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
