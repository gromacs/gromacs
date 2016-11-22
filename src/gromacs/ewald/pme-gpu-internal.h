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
 *
 * \brief This file contains internal function definitions for performing the PME calculations on GPU.
 * These are not meant to be exposed outside of the PME GPU code.
 * As of now, their bodies are still in the common pme-gpu.cpp files.
 *
 * \author Aleksei Iupinov <a.yupinov@gmail.com>
 * \ingroup module_ewald
 */

#ifndef PMEGPUINTERNAL_H
#define PMEGPUINTERNAL_H

#include "gromacs/fft/fft.h"                   // for the enum gmx_fft_direction
#include "gromacs/gpu_utils/gpu_macros.h"      // for the CUDA_FUNC_ macros

#include "pme-gpu-types.h"                     // for the inline functions accessing pme_gpu_t members

struct gmx_hw_info_t;
struct gmx_gpu_opt_t;
struct gmx_pme_t;                              // only used in pme_gpu_reinit
struct gmx_wallclock_gpu_t;                    // only used in pme_gpu_get_timings

/* Some general defines for PME GPU behaviour follow.
 * Some of the might be possible to turn into booleans.
 */

/*! \brief \libinternal
 * 0: The atom data GPU buffers are sized precisely according to the number of atoms.
 *    The atom index checks in the spread/gather code potentially hinder the performance.
 * 1: The atom data GPU buffers are padded with zeroes so that the number of atoms
 *    potentially fitting is divisible by PME_SPREADGATHER_ATOMS_PER_BLOCK (currently always 8).
 *    The atom index checks are not performed. There should be a performance win, but how big is it, remains to be seen.
 *    Additional cudaMemsetAsync calls are done occasionally (only charges/coordinates; spline data is always recalculated now).
 * \todo Estimate performance differences, consider turning into a boolean.
 */
#define PME_GPU_USE_PADDING 1

/*! \brief \libinternal
 * 0: Atoms with zero charges are processed by PME. Could introduce some overhead.
 * 1: Atoms with zero charges are not processed by PME. Adds branching to the spread/gather.
 *    Could be good for performance in specific systems with lots of neutral atoms.
 * \todo Estimate performance differences, consider turning into a boolean.
 */
#define PME_GPU_SKIP_ZEROES 0

/*! \brief \libinternal
 * This is a number of output floats of PME solve.
 * 6 floats for symmetric virial matrix + 1 float for reciprocal energy.
 * Better to have a magic number like this defined in one place.
 * Works better as a define - for more concise CUDA kernel.
 */
#define PME_GPU_VIRIAL_AND_ENERGY_COUNT 7

/* A block of CUDA-only functions that live in pme.cu */

/*! \libinternal \brief
 * Returns the number of atoms per execution block in the atom data layout.
 * Depends on CUDA-specific block sizes, needed for the atom data padding.
 *
 * \param[in] pmeGPU            The PME GPU structure.
 * \returns   Number of atoms in a single GPU atom data chunk.
 */
CUDA_FUNC_QUALIFIER int pme_gpu_get_atom_data_alignment(const pme_gpu_t *CUDA_FUNC_ARGUMENT(pmeGPU)) CUDA_FUNC_TERM_WITH_RETURN(1)

/*! \libinternal \brief
 * Synchronizes the current step, waiting for the GPU kernels/transfers to finish.
 *
 * \param[in] pmeGPU            The PME GPU structure.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_synchronize(const pme_gpu_t *CUDA_FUNC_ARGUMENT(pmeGPU)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * Allocates the fixed size energy and virial buffer both on GPU and CPU.
 *
 * \param[in] pmeGPU            The PME GPU structure.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_alloc_energy_virial(const pme_gpu_t *CUDA_FUNC_ARGUMENT(pmeGPU)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * Frees the energy and virial memory both on GPU and CPU.
 *
 * \param[in] pmeGPU            The PME GPU structure.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_free_energy_virial(pme_gpu_t *CUDA_FUNC_ARGUMENT(pmeGPU)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * Clears the energy and virial memory on GPU with 0.
 * Should be called at the end of the energy/virial calculation step.
 *
 * \param[in] pmeGPU            The PME GPU structure.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_clear_energy_virial(const pme_gpu_t *CUDA_FUNC_ARGUMENT(pmeGPU)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * Reallocates and copies the pre-computed B-spline values to the GPU.
 * FIXME: currently uses just a global memory, could be using texture memory/ldg.
 *
 * \param[in] pmeGPU             The PME GPU structure.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_realloc_and_copy_bspline_values(const pme_gpu_t *CUDA_FUNC_ARGUMENT(pmeGPU)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * Frees the pre-computed B-spline values on the GPU (and the transfer CPU buffers).
 *
 * \param[in] pmeGPU             The PME GPU structure.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_free_bspline_values(const pme_gpu_t *CUDA_FUNC_ARGUMENT(pmeGPU)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * Reallocates the GPU buffer for the PME forces.
 *
 * \param[in] pmeGPU             The PME GPU structure.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_realloc_forces(const pme_gpu_t *CUDA_FUNC_ARGUMENT(pmeGPU)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * Frees the GPU buffer for the PME forces.
 *
 * \param[in] pmeGPU             The PME GPU structure.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_free_forces(const pme_gpu_t *CUDA_FUNC_ARGUMENT(pmeGPU)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * Copies the forces from the CPU buffer to the GPU (to reduce them with the PME GPU gathered forces).
 * To be called e.g. after the bonded calculations.
 *
 * \param[in] pmeGPU             The PME GPU structure.
 * \param[in] h_forces           The input forces rvec buffer.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_copy_input_forces(const pme_gpu_t *CUDA_FUNC_ARGUMENT(pmeGPU),
                                                   const float     *CUDA_FUNC_ARGUMENT(h_forces)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * Copies the forces from the GPU to the CPU buffer. To be called after the gathering stage.
 *
 * \param[in] pmeGPU             The PME GPU structure.
 * \param[out] h_forces          The output forces rvec buffer.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_copy_output_forces(const pme_gpu_t *CUDA_FUNC_ARGUMENT(pmeGPU),
                                                    float           *CUDA_FUNC_ARGUMENT(h_forces)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * Waits for the PME GPU output forces copying to the CPU buffer to finish.
 *
 * \param[in] pmeGPU            The PME GPU structure.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_sync_output_forces(const pme_gpu_t *CUDA_FUNC_ARGUMENT(pmeGPU)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * Reallocates the input coordinates buffer on the GPU (and clears the padded part if needed).
 *
 * \param[in] pmeGPU            The PME GPU structure.
 *
 * Needs to be called on every DD step/in the beginning.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_realloc_coordinates(const pme_gpu_t *CUDA_FUNC_ARGUMENT(pmeGPU)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * Copies the input coordinates from the CPU buffer onto the GPU.
 *
 * \param[in] pmeGPU            The PME GPU structure.
 * \param[in] h_coordinates     Input coordinates (XYZ rvec array).
 *
 * Needs to be called every MD step. The coordinates are then used in the spline calculation.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_copy_input_coordinates(const pme_gpu_t *CUDA_FUNC_ARGUMENT(pmeGPU),
                                                        const rvec      *CUDA_FUNC_ARGUMENT(h_coordinates)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * Frees the coordinates on the GPU.
 *
 * \param[in] pmeGPU            The PME GPU structure.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_free_coordinates(const pme_gpu_t *CUDA_FUNC_ARGUMENT(pmeGPU)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * Reallocates the buffer on the GPU and copies the charges/coefficients from the CPU buffer.
 * Clears the padded part if needed.
 *
 * \param[in] pmeGPU            The PME GPU structure.
 * \param[in] h_coefficients    The input atom charges/coefficients.
 *
 * Does not need to be done every MD step, only whenever the local charges change.
 * (So, in the beginning of the run, or on DD step).
 */
CUDA_FUNC_QUALIFIER void pme_gpu_realloc_and_copy_input_coefficients(const pme_gpu_t *CUDA_FUNC_ARGUMENT(pmeGPU),
                                                                     const float     *CUDA_FUNC_ARGUMENT(h_coefficients)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * Frees the charges/coefficients on the GPU.
 *
 * \param[in] pmeGPU             The PME GPU structure.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_free_coefficients(const pme_gpu_t *CUDA_FUNC_ARGUMENT(pmeGPU)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * Reallocates the buffers on the GPU for the atoms spline data.
 *
 * \param[in] pmeGPU            The PME GPU structure.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_realloc_spline_data(const pme_gpu_t *CUDA_FUNC_ARGUMENT(pmeGPU)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * Frees the buffers on the GPU for the atoms spline data.
 *
 * \param[in] pmeGPU            The PME GPU structure.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_free_spline_data(const pme_gpu_t *CUDA_FUNC_ARGUMENT(pmeGPU)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * Reallocates the buffer on the GPU for the particle gridline indices.
 *
 * \param[in] pmeGPU            The PME GPU structure.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_realloc_grid_indices(const pme_gpu_t *CUDA_FUNC_ARGUMENT(pmeGPU)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * Frees the buffer on the GPU for the particle gridline indices.
 *
 * \param[in] pmeGPU            The PME GPU structure.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_free_grid_indices(const pme_gpu_t *CUDA_FUNC_ARGUMENT(pmeGPU)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * Reallocates the real space grid and the complex reciprocal grid (if needed) on the GPU.
 *
 * \param[in] pmeGPU            The PME GPU structure.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_realloc_grids(pme_gpu_t *CUDA_FUNC_ARGUMENT(pmeGPU)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * Frees the real space grid and the complex reciprocal grid (if needed) on the GPU.
 *
 * \param[in] pmeGPU            The PME GPU structure.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_free_grids(const pme_gpu_t *CUDA_FUNC_ARGUMENT(pmeGPU)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * Clears the real space grid on the GPU.
 * Should be called at the end of each MD step.
 *
 * \param[in] pmeGPU            The PME GPU structure.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_clear_grids(const pme_gpu_t *CUDA_FUNC_ARGUMENT(pmeGPU)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * Reallocates and copies the pre-computed fractional coordinates' shifts to the GPU.
 *
 * \param[in] pmeGPU            The PME GPU structure.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_realloc_and_copy_fract_shifts(pme_gpu_t *CUDA_FUNC_ARGUMENT(pmeGPU)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * Frees the pre-computed fractional coordinates' shifts on the GPU.
 *
 * \param[in] pmeGPU            The PME GPU structure.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_free_fract_shifts(const pme_gpu_t *CUDA_FUNC_ARGUMENT(pmeGPU)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * Waits for the PME GPU output virial/energy copying to the intermediate CPU buffer to finish.
 *
 * \param[in] pmeGPU  The PME GPU structure.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_sync_output_energy_virial(const pme_gpu_t *CUDA_FUNC_ARGUMENT(pmeGPU)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * Waits for the PME GPU grid copying to the host-side buffer to finish.
 *
 * \param[in] pmeGPU  The PME GPU structure.
 * \param[in] dir     The FFT direction.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_sync_grid(const pme_gpu_t        *CUDA_FUNC_ARGUMENT(pmeGPU),
                                           const gmx_fft_direction CUDA_FUNC_ARGUMENT(dir)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * Does the one-time GPU-framework specific PME initialization.
 * Should be called early.
 *
 * \param[in] pmeGPU  The PME GPU structure.
 * \param[in] hwinfo  The hardware information structure.
 * \param[in] gpu_opt The GPU information structure.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_init_specific(pme_gpu_t           *CUDA_FUNC_ARGUMENT(pmeGPU),
                                               const gmx_hw_info_t *CUDA_FUNC_ARGUMENT(hwinfo),
                                               const gmx_gpu_opt_t *CUDA_FUNC_ARGUMENT(gpu_opt)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * Destroys the PME GPU-framework specific data.
 * Should be called last in the PME GPU destructor.
 *
 * \param[in] pmeGPU  The PME GPU structure.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_destroy_specific(const pme_gpu_t *CUDA_FUNC_ARGUMENT(pmeGPU)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * Initializes the PME GPU synchronization events.
 *
 * \param[in] pmeGPU  The PME GPU structure.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_init_sync_events(const pme_gpu_t *CUDA_FUNC_ARGUMENT(pmeGPU)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * Destroys the PME GPU synchronization events.
 *
 * \param[in] pmeGPU  The PME GPU structure.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_destroy_sync_events(const pme_gpu_t *CUDA_FUNC_ARGUMENT(pmeGPU)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * Initializes the CUDA FFT structures.
 *
 * \param[in] pmeGPU  The PME GPU structure.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_reinit_3dfft(const pme_gpu_t *CUDA_FUNC_ARGUMENT(pmeGPU)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * Destroys the CUDA FFT structures.
 *
 * \param[in] pmeGPU  The PME GPU structure.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_destroy_3dfft(const pme_gpu_t *CUDA_FUNC_ARGUMENT(pmeGPU)) CUDA_FUNC_TERM

/* Several CUDA event-based timing functions that live in pme-timings.cu */

/*! \libinternal \brief
 * Allocates and initializes the PME GPU timings.
 *
 * \param[in] pmeGPU         The PME GPU structure.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_init_timings(const pme_gpu_t *CUDA_FUNC_ARGUMENT(pmeGPU)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * Destroys the PME GPU timings.
 *
 * \param[in] pmeGPU         The PME GPU structure.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_destroy_timings(const pme_gpu_t *CUDA_FUNC_ARGUMENT(pmeGPU)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * Finalizes all the PME GPU stage timings for the current step. Should be called at the end of every step.
 *
 * \param[in] pmeGPU         The PME GPU structure.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_update_timings(const pme_gpu_t *CUDA_FUNC_ARGUMENT(pmeGPU)) CUDA_FUNC_TERM

/*! \brief
 * Resets the PME GPU timings. To be called at the reset step.
 *
 * \param[in] pmeGPU         The PME GPU structure.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_reset_timings(const pme_gpu_t *CUDA_FUNC_ARGUMENT(pmeGPU)) CUDA_FUNC_TERM

/*! \libinternal \brief
 * Copies the PME GPU timings to the gmx_wallclock_gpu_t structure (for log output). To be called at the run end.
 *
 * \param[in] pmeGPU         The PME GPU structure.
 * \param[in] timings        The gmx_wallclock_gpu_t structure.
 */
CUDA_FUNC_QUALIFIER void pme_gpu_get_timings(const pme_gpu_t      *CUDA_FUNC_ARGUMENT(pmeGPU),
                                             gmx_wallclock_gpu_t **CUDA_FUNC_ARGUMENT(timings)) CUDA_FUNC_TERM

/* Separate PME GPU stages, living in ther own *.cu files */

#include "pme-internal.h"
//#include "gromacs/math/gmxcomplex.h"
// A GPU counterpart to the spread_on_grid
CUDA_FUNC_QUALIFIER void pme_gpu_spread(const gmx_pme_t *CUDA_FUNC_ARGUMENT(pme),
                                        pme_atomcomm_t  *CUDA_FUNC_ARGUMENT(atc),
                                        const int        CUDA_FUNC_ARGUMENT(grid_index),
                                        pmegrid_t       *CUDA_FUNC_ARGUMENT(pmegrid),
                                        const gmx_bool   CUDA_FUNC_ARGUMENT(bCalcSplines),
                                        const gmx_bool   CUDA_FUNC_ARGUMENT(bSpread)) CUDA_FUNC_TERM

// A GPU counterpart to gmx_parallel_3dfft_execute
CUDA_FUNC_QUALIFIER void pme_gpu_3dfft(const pme_gpu_t       *CUDA_FUNC_ARGUMENT(pmeGPU),
                                       enum gmx_fft_direction CUDA_FUNC_ARGUMENT(dir),
                                       const int              CUDA_FUNC_ARGUMENT(grid_index)) CUDA_FUNC_TERM

// A GPU counterpart to the solve_pme_yzx
CUDA_FUNC_QUALIFIER void pme_gpu_solve(
        gmx_pme_t     *CUDA_FUNC_ARGUMENT(pme),
        t_complex     *CUDA_FUNC_ARGUMENT(grid),
        const gmx_bool CUDA_FUNC_ARGUMENT(bEnerVir)) CUDA_FUNC_TERM

// A GPU counterpart to the gather_f_bsplines
//* \param[in] forces         The pointer to the host-side array of particle forces.
//*                           It will be used for output, but can also be used for input,
//*                           if bClearForces is passed as false to the pme_gpu_launch_gather.

CUDA_FUNC_QUALIFIER void pme_gpu_gather(const gmx_pme_t *CUDA_FUNC_ARGUMENT(pme),
                                        const gmx_bool   CUDA_FUNC_ARGUMENT(bOverwriteForces),
                                        float           *CUDA_FUNC_ARGUMENT(h_forces)) CUDA_FUNC_TERM



/* The inlined convenience PME GPU status getters */

/*! \libinternal \brief
 * Tells if PME runs on multiple GPUs with the decomposition.
 *
 * \param[in] pmeGPU         The PME GPU structure.
 * \returns                  TRUE if PME runs on multiple GPUs, FALSE otherwise.
 */
gmx_inline bool pme_gpu_uses_dd(const pme_gpu_t *pmeGPU)
{
    return !pmeGPU->settings.useDecomposition;
}

/*! \libinternal \brief
 * Tells if PME performs the gathering stage on GPU.
 *
 * \param[in] pmeGPU         The PME GPU structure.
 * \returns                  TRUE if the gathering is performed on GPU, FALSE otherwise.
 */
gmx_inline bool pme_gpu_performs_gather(const pme_gpu_t *pmeGPU)
{
    return pmeGPU->settings.performGPUGather;
}

/*! \libinternal \brief
 * Tells if PME performs the FFT stages on GPU.
 *
 * \param[in] pmeGPU         The PME GPU structure.
 * \returns                  TRUE if FFT is performed on GPU, FALSE otherwise.
 */
gmx_inline bool pme_gpu_performs_FFT(const pme_gpu_t *pmeGPU)
{
    return pmeGPU->settings.performGPUFFT;
}

/*! \libinternal \brief
 * Tells if PME performs the grid (un-)wrapping on GPU.
 *
 * \param[in] pmeGPU         The PME GPU structure.
 * \returns                  TRUE if (un-)wrapping is performed on GPU, FALSE otherwise.
 */
gmx_inline bool pme_gpu_performs_wrapping(const pme_gpu_t *pmeGPU)
{
    return pmeGPU->settings.useDecomposition;
}

/*! \libinternal \brief
 * Tells if PME performs the grid solving on GPU.
 *
 * \param[in] pmeGPU         The PME GPU structure.
 * \returns                  TRUE if solving is performed on GPU, FALSE otherwise.
 */
gmx_inline bool pme_gpu_performs_solve(const pme_gpu_t *pmeGPU)
{
    return pmeGPU->settings.performGPUSolve;
}

/* A block of C++ functions that live in pme-gpu-internal.cpp */

/*! \libinternal \brief
 * Returns the output virial and energy of the PME solving.
 * Should be called after pme_gpu_finish_step.
 *
 * \param[in] pmeGPU             The PME GPU structure.
 * \param[out] energy            The output energy.
 * \param[out] virial            The output virial matrix.
 */
void pme_gpu_get_energy_virial(const pme_gpu_t *pmeGPU, real *energy, matrix virial);

/*! \libinternal \brief
 * Starts the PME GPU step (copies coordinates onto GPU, possibly sets the unit cell parameters).
 * Does nothing if PME on GPU is disabled.
 *
 * \param[in] pmeGPU         The PME GPU structure.
 * \param[in] box            The unit cell box which does not necessarily change every step (only with pressure coupling enabled).
 *                           Currently it is simply compared with the previous one to determine if it needs to be updated.
 * \param[in] h_coordinates  Input coordinates (XYZ rvec array).
 */
void pme_gpu_start_step(pme_gpu_t *pmeGPU, const matrix box, const rvec *h_coordinates);

/*! \libinternal \brief
 * Finishes the PME GPU step, waiting for the output forces and/or energy/virial to be copied to the host.
 *
 * \param[in] pmeGPU         The PME GPU structure.
 * \param[in] bCalcForces    The left-over flag from the CPU code which tells the function to copy the forces to the CPU side. Should be passed to the launch call instead.
 * \param[in] bCalcEnerVir   The left-over flag from the CPU code which tells the function to copy the energy/virial to the CPU side. Should be passed to the launch call instead.
 */
void pme_gpu_finish_step(const pme_gpu_t *pmeGPU,
                         const bool       bCalcForces,
                         const bool       bCalcEnerVir);

/*! \libinternal \brief
 * (Re-)initializes the PME GPU data at the beginning of the run or on DLB.
 *
 * \param[in] pme     The PME structure.
 * \param[in] hwinfo  The hardware information structure.
 * \param[in] gpu_opt The GPU information structure.
 */
void pme_gpu_reinit(gmx_pme_t           *pme,
                    const gmx_hw_info_t *hwinfo,
                    const gmx_gpu_opt_t *gpu_opt);

/*! \libinternal \brief
 * Destroys the PME GPU data at the end of the run.
 *
 * \param[in] pmeGPU     The PME GPU structure.
 */
void pme_gpu_destroy(pme_gpu_t *pmeGPU);

/*! \libinternal \brief
 * Reallocates the local atoms data (charges, coordinates, etc.). Copies the charges to the GPU.
 *
 * \param[in] pmeGPU         The PME GPU structure.
 * \param[in] nAtoms         The number of particles.
 * \param[in] coefficients   The pointer to the host-side array of particle charges.
 *
 * This is a function that should only be called in the beginning of the run and on domain decomposition.
 * Should be called before the pme_gpu_set_io_ranges.
 */
void pme_gpu_reinit_atoms(pme_gpu_t        *pmeGPU,
                          const int         nAtoms,
                          const real       *coefficients);

/*! \brief
 * Finds out if PME is set to run on GPU. Not to be called on per-step basis.
 * TODO: remove this entirely with introduction of task/hardware scheduler.
 *
 * \param[in] pme  The PME structure.
 * \returns        True if PME runs on GPU currently, false otherwise.
 */
inline bool pme_gpu_active(const gmx_pme_t *pme)
{
    /* Something to think about: should this function be called from all the CUDA_FUNC_QUALIFIER functions?
     * In other words, should we plan for dynamic toggling of the PME GPU?
     */
    return (pme != NULL) && pme->useGPU;
}

#endif // PMEGPUINTERNAL_H
