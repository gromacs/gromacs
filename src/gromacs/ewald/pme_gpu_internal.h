/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2017,2018,2019, by the GROMACS development team, led by
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
 * As of now, their bodies are still in the common pme_gpu.cpp files.
 *
 * \author Aleksei Iupinov <a.yupinov@gmail.com>
 * \ingroup module_ewald
 */

#ifndef GMX_EWALD_PME_GPU_INTERNAL_H
#define GMX_EWALD_PME_GPU_INTERNAL_H

#include "gromacs/fft/fft.h"              // for the gmx_fft_direction enum
#include "gromacs/gpu_utils/gpu_macros.h" // for the GPU_FUNC_ macros
#include "gromacs/utility/arrayref.h"

#include "pme_gpu_types_host.h" // for the inline functions accessing PmeGpu members

struct gmx_hw_info_t;
struct gmx_gpu_opt_t;
struct gmx_pme_t; // only used in pme_gpu_reinit
struct gmx_wallclock_gpu_pme_t;
class PmeAtomComm;
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
};               // TODO move this into new and shiny pme.h (pme-types.h?)

//! PME grid dimension ordering (from major to minor)
enum class GridOrdering
{
    YZX,
    XYZ
};

/*! \libinternal \brief
 * Returns the number of atoms per chunk in the atom charges/coordinates data layout.
 * Depends on CUDA-specific block sizes, needed for the atom data padding.
 *
 * \param[in] pmeGpu            The PME GPU structure.
 * \returns   Number of atoms in a single GPU atom data chunk.
 */
int pme_gpu_get_atom_data_alignment(const PmeGpu* pmeGpu);

/*! \libinternal \brief
 * Returns the number of atoms per chunk in the atom spline theta/dtheta data layout.
 *
 * \param[in] pmeGpu            The PME GPU structure.
 * \returns   Number of atoms in a single GPU atom spline data chunk.
 */
int pme_gpu_get_atoms_per_warp(const PmeGpu* pmeGpu);

/*! \libinternal \brief
 * Synchronizes the current computation, waiting for the GPU kernels/transfers to finish.
 *
 * \param[in] pmeGpu            The PME GPU structure.
 */
GPU_FUNC_QUALIFIER void pme_gpu_synchronize(const PmeGpu* GPU_FUNC_ARGUMENT(pmeGpu)) GPU_FUNC_TERM;

/*! \libinternal \brief
 * Allocates the fixed size energy and virial buffer both on GPU and CPU.
 *
 * \param[in,out] pmeGpu            The PME GPU structure.
 */
void pme_gpu_alloc_energy_virial(PmeGpu* pmeGpu);

/*! \libinternal \brief
 * Frees the energy and virial memory both on GPU and CPU.
 *
 * \param[in] pmeGpu            The PME GPU structure.
 */
void pme_gpu_free_energy_virial(PmeGpu* pmeGpu);

/*! \libinternal \brief
 * Clears the energy and virial memory on GPU with 0.
 * Should be called at the end of PME computation which returned energy/virial.
 *
 * \param[in] pmeGpu            The PME GPU structure.
 */
void pme_gpu_clear_energy_virial(const PmeGpu* pmeGpu);

/*! \libinternal \brief
 * Reallocates and copies the pre-computed B-spline values to the GPU.
 *
 * \param[in,out] pmeGpu             The PME GPU structure.
 */
void pme_gpu_realloc_and_copy_bspline_values(PmeGpu* pmeGpu);

/*! \libinternal \brief
 * Frees the pre-computed B-spline values on the GPU (and the transfer CPU buffers).
 *
 * \param[in] pmeGpu             The PME GPU structure.
 */
void pme_gpu_free_bspline_values(const PmeGpu* pmeGpu);

/*! \libinternal \brief
 * Reallocates the GPU buffer for the PME forces.
 *
 * \param[in] pmeGpu             The PME GPU structure.
 */
void pme_gpu_realloc_forces(PmeGpu* pmeGpu);

/*! \libinternal \brief
 * Frees the GPU buffer for the PME forces.
 *
 * \param[in] pmeGpu             The PME GPU structure.
 */
void pme_gpu_free_forces(const PmeGpu* pmeGpu);

/*! \libinternal \brief
 * Copies the forces from the CPU buffer to the GPU (to reduce them with the PME GPU gathered
 * forces). To be called e.g. after the bonded calculations.
 *
 * \param[in] pmeGpu             The PME GPU structure.
 */
void pme_gpu_copy_input_forces(PmeGpu* pmeGpu);

/*! \libinternal \brief
 * Copies the forces from the GPU to the CPU buffer. To be called after the gathering stage.
 *
 * \param[in] pmeGpu             The PME GPU structure.
 */
void pme_gpu_copy_output_forces(PmeGpu* pmeGpu);

/*! \libinternal \brief
 * Checks whether work in the PME GPU stream has completed.
 *
 * \param[in] pmeGpu            The PME GPU structure.
 *
 * \returns                     True if work in the PME stream has completed.
 */
bool pme_gpu_stream_query(const PmeGpu* pmeGpu);

/*! \libinternal \brief
 * Reallocates the input coordinates buffer on the GPU (and clears the padded part if needed).
 *
 * \param[in] pmeGpu            The PME GPU structure.
 *
 * Needs to be called on every DD step/in the beginning.
 */
void pme_gpu_realloc_coordinates(const PmeGpu* pmeGpu);

/*! \libinternal \brief
 * Frees the coordinates on the GPU.
 *
 * \param[in] pmeGpu            The PME GPU structure.
 */
void pme_gpu_free_coordinates(const PmeGpu* pmeGpu);

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
void pme_gpu_realloc_and_copy_input_coefficients(const PmeGpu* pmeGpu, const float* h_coefficients);

/*! \libinternal \brief
 * Frees the charges/coefficients on the GPU.
 *
 * \param[in] pmeGpu             The PME GPU structure.
 */
void pme_gpu_free_coefficients(const PmeGpu* pmeGpu);

/*! \libinternal \brief
 * Reallocates the buffers on the GPU and the host for the atoms spline data.
 *
 * \param[in,out] pmeGpu            The PME GPU structure.
 */
void pme_gpu_realloc_spline_data(PmeGpu* pmeGpu);

/*! \libinternal \brief
 * Frees the buffers on the GPU for the atoms spline data.
 *
 * \param[in] pmeGpu            The PME GPU structure.
 */
void pme_gpu_free_spline_data(const PmeGpu* pmeGpu);

/*! \libinternal \brief
 * Reallocates the buffers on the GPU and the host for the particle gridline indices.
 *
 * \param[in,out] pmeGpu            The PME GPU structure.
 */
void pme_gpu_realloc_grid_indices(PmeGpu* pmeGpu);

/*! \libinternal \brief
 * Frees the buffer on the GPU for the particle gridline indices.
 *
 * \param[in] pmeGpu            The PME GPU structure.
 */
void pme_gpu_free_grid_indices(const PmeGpu* pmeGpu);

/*! \libinternal \brief
 * Reallocates the real space grid and the complex reciprocal grid (if needed) on the GPU.
 *
 * \param[in] pmeGpu            The PME GPU structure.
 */
void pme_gpu_realloc_grids(PmeGpu* pmeGpu);

/*! \libinternal \brief
 * Frees the real space grid and the complex reciprocal grid (if needed) on the GPU.
 *
 * \param[in] pmeGpu            The PME GPU structure.
 */
void pme_gpu_free_grids(const PmeGpu* pmeGpu);

/*! \libinternal \brief
 * Clears the real space grid on the GPU.
 * Should be called at the end of each computation.
 *
 * \param[in] pmeGpu            The PME GPU structure.
 */
void pme_gpu_clear_grids(const PmeGpu* pmeGpu);

/*! \libinternal \brief
 * Reallocates and copies the pre-computed fractional coordinates' shifts to the GPU.
 *
 * \param[in] pmeGpu            The PME GPU structure.
 */
void pme_gpu_realloc_and_copy_fract_shifts(PmeGpu* pmeGpu);

/*! \libinternal \brief
 * Frees the pre-computed fractional coordinates' shifts on the GPU.
 *
 * \param[in] pmeGpu            The PME GPU structure.
 */
void pme_gpu_free_fract_shifts(const PmeGpu* pmeGpu);

/*! \libinternal \brief
 * Copies the input real-space grid from the host to the GPU.
 *
 * \param[in] pmeGpu   The PME GPU structure.
 * \param[in] h_grid   The host-side grid buffer.
 */
void pme_gpu_copy_input_gather_grid(const PmeGpu* pmeGpu, float* h_grid);

/*! \libinternal \brief
 * Copies the output real-space grid from the GPU to the host.
 *
 * \param[in] pmeGpu   The PME GPU structure.
 * \param[out] h_grid  The host-side grid buffer.
 */
void pme_gpu_copy_output_spread_grid(const PmeGpu* pmeGpu, float* h_grid);

/*! \libinternal \brief
 * Copies the spread output spline data and gridline indices from the GPU to the host.
 *
 * \param[in] pmeGpu   The PME GPU structure.
 */
void pme_gpu_copy_output_spread_atom_data(const PmeGpu* pmeGpu);

/*! \libinternal \brief
 * Copies the gather input spline data and gridline indices from the host to the GPU.
 *
 * \param[in] pmeGpu   The PME GPU structure.
 */
void pme_gpu_copy_input_gather_atom_data(const PmeGpu* pmeGpu);

/*! \libinternal \brief
 * Waits for the grid copying to the host-side buffer after spreading to finish.
 *
 * \param[in] pmeGpu  The PME GPU structure.
 */
void pme_gpu_sync_spread_grid(const PmeGpu* pmeGpu);

/*! \libinternal \brief
 * Does the one-time GPU-framework specific PME initialization.
 * For CUDA, the PME stream is created with the highest priority.
 *
 * \param[in] pmeGpu  The PME GPU structure.
 */
void pme_gpu_init_internal(PmeGpu* pmeGpu);

/*! \libinternal \brief
 * Destroys the PME GPU-framework specific data.
 * Should be called last in the PME GPU destructor.
 *
 * \param[in] pmeGpu  The PME GPU structure.
 */
void pme_gpu_destroy_specific(const PmeGpu* pmeGpu);

/*! \libinternal \brief
 * Initializes the CUDA FFT structures.
 *
 * \param[in] pmeGpu  The PME GPU structure.
 */
void pme_gpu_reinit_3dfft(const PmeGpu* pmeGpu);

/*! \libinternal \brief
 * Destroys the CUDA FFT structures.
 *
 * \param[in] pmeGpu  The PME GPU structure.
 */
void pme_gpu_destroy_3dfft(const PmeGpu* pmeGpu);

/* Several GPU event-based timing functions that live in pme_gpu_timings.cpp */

/*! \libinternal \brief
 * Finalizes all the active PME GPU stage timings for the current computation. Should be called at the end of every computation.
 *
 * \param[in] pmeGpu         The PME GPU structure.
 */
void pme_gpu_update_timings(const PmeGpu* pmeGpu);

/*! \libinternal \brief
 * Updates the internal list of active PME GPU stages (if timings are enabled).
 *
 * \param[in] pmeGpu         The PME GPU data structure.
 */
void pme_gpu_reinit_timings(const PmeGpu* pmeGpu);

/*! \brief
 * Resets the PME GPU timings. To be called at the reset MD step.
 *
 * \param[in] pmeGpu         The PME GPU structure.
 */
void pme_gpu_reset_timings(const PmeGpu* pmeGpu);

/*! \libinternal \brief
 * Copies the PME GPU timings to the gmx_wallclock_gpu_t structure (for log output). To be called at the run end.
 *
 * \param[in] pmeGpu         The PME GPU structure.
 * \param[in] timings        The gmx_wallclock_gpu_pme_t structure.
 */
void pme_gpu_get_timings(const PmeGpu* pmeGpu, gmx_wallclock_gpu_pme_t* timings);

/* The PME stages themselves */

/*! \libinternal \brief
 * A GPU spline computation and charge spreading function.
 *
 * \param[in]  pmeGpu          The PME GPU structure.
 * \param[in]  xReadyOnDevice  Event synchronizer indicating that the coordinates are ready in the device memory;
 *                             can be nullptr when invoked on a separate PME rank or from PME tests.
 * \param[in]  gridIndex       Index of the PME grid - unused, assumed to be 0.
 * \param[out] h_grid          The host-side grid buffer (used only if the result of the spread is expected on the host,
 *                             e.g. testing or host-side FFT)
 * \param[in]  computeSplines  Should the computation of spline parameters and gridline indices be performed.
 * \param[in]  spreadCharges   Should the charges/coefficients be spread on the grid.
 */
GPU_FUNC_QUALIFIER void pme_gpu_spread(const PmeGpu*         GPU_FUNC_ARGUMENT(pmeGpu),
                                       GpuEventSynchronizer* GPU_FUNC_ARGUMENT(xReadyOnDevice),
                                       int                   GPU_FUNC_ARGUMENT(gridIndex),
                                       real*                 GPU_FUNC_ARGUMENT(h_grid),
                                       bool                  GPU_FUNC_ARGUMENT(computeSplines),
                                       bool GPU_FUNC_ARGUMENT(spreadCharges)) GPU_FUNC_TERM;

/*! \libinternal \brief
 * 3D FFT R2C/C2R routine.
 *
 * \param[in]  pmeGpu          The PME GPU structure.
 * \param[in]  direction       Transform direction (real-to-complex or complex-to-real)
 * \param[in]  gridIndex       Index of the PME grid - unused, assumed to be 0.
 */
void pme_gpu_3dfft(const PmeGpu* pmeGpu, enum gmx_fft_direction direction, int gridIndex);

/*! \libinternal \brief
 * A GPU Fourier space solving function.
 *
 * \param[in]     pmeGpu                  The PME GPU structure.
 * \param[in,out] h_grid                  The host-side input and output Fourier grid buffer (used only with testing or host-side FFT)
 * \param[in]     gridOrdering            Specifies the dimenion ordering of the complex grid. TODO: store this information?
 * \param[in]     computeEnergyAndVirial  Tells if the energy and virial computation should also be performed.
 */
GPU_FUNC_QUALIFIER void pme_gpu_solve(const PmeGpu* GPU_FUNC_ARGUMENT(pmeGpu),
                                      t_complex*    GPU_FUNC_ARGUMENT(h_grid),
                                      GridOrdering  GPU_FUNC_ARGUMENT(gridOrdering),
                                      bool GPU_FUNC_ARGUMENT(computeEnergyAndVirial)) GPU_FUNC_TERM;

/*! \libinternal \brief
 * A GPU force gathering function.
 *
 * \param[in]     pmeGpu           The PME GPU structure.
 * \param[in]     forceTreatment   Tells how data in h_forces should be treated.
 *                                 TODO: determine efficiency/balance of host/device-side
 * reductions. \param[in]     h_grid           The host-side grid buffer (used only in testing mode)
 */
GPU_FUNC_QUALIFIER void pme_gpu_gather(PmeGpu*                GPU_FUNC_ARGUMENT(pmeGpu),
                                       PmeForceOutputHandling GPU_FUNC_ARGUMENT(forceTreatment),
                                       const float* GPU_FUNC_ARGUMENT(h_grid)) GPU_FUNC_TERM;

/*! \brief Return pointer to device copy of coordinate data.
 * \param[in] pmeGpu         The PME GPU structure.
 * \returns                  Pointer to coordinate data
 */
GPU_FUNC_QUALIFIER DeviceBuffer<float> pme_gpu_get_kernelparam_coordinates(const PmeGpu* GPU_FUNC_ARGUMENT(pmeGpu))
        GPU_FUNC_TERM_WITH_RETURN(DeviceBuffer<float>{});

/*! \brief Sets the device pointer to coordinate data
 * \param[in] pmeGpu         The PME GPU structure.
 * \param[in] d_x            Pointer to coordinate data
 */
GPU_FUNC_QUALIFIER void pme_gpu_set_kernelparam_coordinates(const PmeGpu* GPU_FUNC_ARGUMENT(pmeGpu),
                                                            DeviceBuffer<float> GPU_FUNC_ARGUMENT(d_x)) GPU_FUNC_TERM;

/*! \brief Return pointer to device copy of force data.
 * \param[in] pmeGpu         The PME GPU structure.
 * \returns                  Pointer to force data
 */
GPU_FUNC_QUALIFIER void* pme_gpu_get_kernelparam_forces(const PmeGpu* GPU_FUNC_ARGUMENT(pmeGpu))
        GPU_FUNC_TERM_WITH_RETURN(nullptr);

/*! \brief Return pointer to GPU stream.
 * \param[in] pmeGpu         The PME GPU structure.
 * \returns                  Pointer to stream object.
 */
GPU_FUNC_QUALIFIER void* pme_gpu_get_stream(const PmeGpu* GPU_FUNC_ARGUMENT(pmeGpu))
        GPU_FUNC_TERM_WITH_RETURN(nullptr);

/*! \brief Return pointer to GPU context (for OpenCL builds).
 * \param[in] pmeGpu         The PME GPU structure.
 * \returns                  Pointer to context object.
 */
GPU_FUNC_QUALIFIER void* pme_gpu_get_context(const PmeGpu* GPU_FUNC_ARGUMENT(pmeGpu))
        GPU_FUNC_TERM_WITH_RETURN(nullptr);

/*! \brief Return pointer to the sync object triggered after the PME force calculation completion
 * \param[in] pmeGpu         The PME GPU structure.
 * \returns                  Pointer to sync object
 */
GPU_FUNC_QUALIFIER GpuEventSynchronizer* pme_gpu_get_forces_ready_synchronizer(
        const PmeGpu* GPU_FUNC_ARGUMENT(pmeGpu)) GPU_FUNC_TERM_WITH_RETURN(nullptr);

/* The inlined convenience PME GPU status getters */

/*! \libinternal \brief
 * Tells if PME runs on multiple GPUs with the decomposition.
 *
 * \param[in] pmeGpu         The PME GPU structure.
 * \returns                  True if PME runs on multiple GPUs, false otherwise.
 */
inline bool pme_gpu_uses_dd(const PmeGpu* pmeGpu)
{
    return !pmeGpu->settings.useDecomposition;
}

/*! \libinternal \brief
 * Tells if PME performs the gathering stage on GPU.
 *
 * \param[in] pmeGpu         The PME GPU structure.
 * \returns                  True if the gathering is performed on GPU, false otherwise.
 */
inline bool pme_gpu_performs_gather(const PmeGpu* pmeGpu)
{
    return pmeGpu->settings.performGPUGather;
}

/*! \libinternal \brief
 * Tells if PME performs the FFT stages on GPU.
 *
 * \param[in] pmeGpu         The PME GPU structure.
 * \returns                  True if FFT is performed on GPU, false otherwise.
 */
inline bool pme_gpu_performs_FFT(const PmeGpu* pmeGpu)
{
    return pmeGpu->settings.performGPUFFT;
}

/*! \libinternal \brief
 * Tells if PME performs the grid (un-)wrapping on GPU.
 *
 * \param[in] pmeGpu         The PME GPU structure.
 * \returns                  True if (un-)wrapping is performed on GPU, false otherwise.
 */
inline bool pme_gpu_performs_wrapping(const PmeGpu* pmeGpu)
{
    return pmeGpu->settings.useDecomposition;
}

/*! \libinternal \brief
 * Tells if PME performs the grid solving on GPU.
 *
 * \param[in] pmeGpu         The PME GPU structure.
 * \returns                  True if solving is performed on GPU, false otherwise.
 */
inline bool pme_gpu_performs_solve(const PmeGpu* pmeGpu)
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
inline void pme_gpu_set_testing(PmeGpu* pmeGpu, bool testing)
{
    if (pmeGpu)
    {
        pmeGpu->settings.copyAllOutputs = testing;
        pmeGpu->settings.transferKind = testing ? GpuApiCallBehavior::Sync : GpuApiCallBehavior::Async;
    }
}

/*! \libinternal \brief
 * Tells if PME is in the testing mode.
 *
 * \param[in] pmeGpu             The PME GPU structure.
 * \returns                      true if testing mode is enabled, false otherwise.
 */
inline bool pme_gpu_is_testing(const PmeGpu* pmeGpu)
{
    return pmeGpu->settings.copyAllOutputs;
}

/* A block of C++ functions that live in pme_gpu_internal.cpp */

/*! \libinternal \brief
 * Returns the energy and virial GPU outputs, useful for testing.
 *
 * It is the caller's responsibility to be aware of whether the GPU
 * handled the solve stage.
 *
 * \param[in] pme                The PME structure.
 * \param[out] output            Pointer to output where energy and virial should be stored.
 */
GPU_FUNC_QUALIFIER void pme_gpu_getEnergyAndVirial(const gmx_pme_t& GPU_FUNC_ARGUMENT(pme),
                                                   PmeOutput* GPU_FUNC_ARGUMENT(output)) GPU_FUNC_TERM;

/*! \libinternal \brief
 * Returns the GPU outputs (forces, energy and virial)
 *
 * \param[in] pme                The PME structure.
 * \param[in] flags              The combination of flags that affected this PME computation.
 *                               The flags are the GMX_PME_ flags from pme.h.
 * \returns                      The output object.
 */
GPU_FUNC_QUALIFIER PmeOutput pme_gpu_getOutput(const gmx_pme_t& GPU_FUNC_ARGUMENT(pme),
                                               int              GPU_FUNC_ARGUMENT(flags))
        GPU_FUNC_TERM_WITH_RETURN(PmeOutput{});

/*! \libinternal \brief
 * Updates the unit cell parameters. Does not check if update is necessary - that is done in pme_gpu_prepare_computation().
 *
 * \param[in] pmeGpu         The PME GPU structure.
 * \param[in] box            The unit cell box.
 */
GPU_FUNC_QUALIFIER void pme_gpu_update_input_box(PmeGpu*      GPU_FUNC_ARGUMENT(pmeGpu),
                                                 const matrix GPU_FUNC_ARGUMENT(box)) GPU_FUNC_TERM;

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
void pme_gpu_finish_computation(const PmeGpu* pmeGpu);

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
GPU_FUNC_QUALIFIER void pme_gpu_transform_spline_atom_data(const PmeGpu* GPU_FUNC_ARGUMENT(pmeGpu),
                                                           const PmeAtomComm* GPU_FUNC_ARGUMENT(atc),
                                                           PmeSplineDataType GPU_FUNC_ARGUMENT(type),
                                                           int GPU_FUNC_ARGUMENT(dimIndex),
                                                           PmeLayoutTransform GPU_FUNC_ARGUMENT(transform)) GPU_FUNC_TERM;

/*! \libinternal \brief
 * Gets a unique index to an element in a spline parameter buffer (theta/dtheta),
 * which is laid out for GPU spread/gather kernels. The index is wrt the execution block,
 * in range(0, atomsPerBlock * order * DIM).
 * This is a wrapper, only used in unit tests.
 * \param[in] order            PME order
 * \param[in] splineIndex      Spline contribution index (from 0 to \p order - 1)
 * \param[in] dimIndex         Dimension index (from 0 to 2)
 * \param[in] atomIndex        Atom index wrt the block.
 * \param[in] atomsPerWarp     Number of atoms processed by a warp.
 *
 * \returns Index into theta or dtheta array using GPU layout.
 */
int getSplineParamFullIndex(int order, int splineIndex, int dimIndex, int atomIndex, int atomsPerWarp);

/*! \libinternal \brief
 * Get the normal/padded grid dimensions of the real-space PME grid on GPU. Only used in tests.
 *
 * \param[in] pmeGpu             The PME GPU structure.
 * \param[out] gridSize          Pointer to the grid dimensions to fill in.
 * \param[out] paddedGridSize    Pointer to the padded grid dimensions to fill in.
 */
GPU_FUNC_QUALIFIER void pme_gpu_get_real_grid_sizes(const PmeGpu* GPU_FUNC_ARGUMENT(pmeGpu),
                                                    gmx::IVec*    GPU_FUNC_ARGUMENT(gridSize),
                                                    gmx::IVec* GPU_FUNC_ARGUMENT(paddedGridSize)) GPU_FUNC_TERM;

/*! \libinternal \brief
 * (Re-)initializes the PME GPU data at the beginning of the run or on DLB.
 *
 * \param[in,out] pme             The PME structure.
 * \param[in]     gpuInfo         The GPU information structure.
 * \param[in]     pmeGpuProgram   The PME GPU program data
 * \throws gmx::NotImplementedError if this generally valid PME structure is not valid for GPU runs.
 */
GPU_FUNC_QUALIFIER void pme_gpu_reinit(gmx_pme_t*               GPU_FUNC_ARGUMENT(pme),
                                       const gmx_device_info_t* GPU_FUNC_ARGUMENT(gpuInfo),
                                       PmeGpuProgramHandle GPU_FUNC_ARGUMENT(pmeGpuProgram)) GPU_FUNC_TERM;

/*! \libinternal \brief
 * Destroys the PME GPU data at the end of the run.
 *
 * \param[in] pmeGpu     The PME GPU structure.
 */
GPU_FUNC_QUALIFIER void pme_gpu_destroy(PmeGpu* GPU_FUNC_ARGUMENT(pmeGpu)) GPU_FUNC_TERM;

/*! \libinternal \brief
 * Reallocates the local atoms data (charges, coordinates, etc.). Copies the charges to the GPU.
 *
 * \param[in] pmeGpu    The PME GPU structure.
 * \param[in] nAtoms    The number of particles.
 * \param[in] charges   The pointer to the host-side array of particle charges.
 *
 * This is a function that should only be called in the beginning of the run and on domain
 * decomposition. Should be called before the pme_gpu_set_io_ranges.
 */
GPU_FUNC_QUALIFIER void pme_gpu_reinit_atoms(PmeGpu*     GPU_FUNC_ARGUMENT(pmeGpu),
                                             int         GPU_FUNC_ARGUMENT(nAtoms),
                                             const real* GPU_FUNC_ARGUMENT(charges)) GPU_FUNC_TERM;

/*! \brief \libinternal
 * The PME GPU reinitialization function that is called both at the end of any PME computation and on any load balancing.
 *
 * This clears the device-side working buffers in preparation for new computation.
 *
 * \param[in] pmeGpu            The PME GPU structure.
 */
void pme_gpu_reinit_computation(const PmeGpu* pmeGpu);

/*! \brief
 * Blocks until PME GPU tasks are completed, and gets the output forces and virial/energy
 * (if they were to be computed).
 *
 * \param[in]  pme            The PME data structure.
 * \param[in]  flags          The combination of flags to affect this PME computation.
 *                            The flags are the GMX_PME_ flags from pme.h.
 * \param[out] wcycle         The wallclock counter.
 * \return     The output forces, energy and virial
 */
GPU_FUNC_QUALIFIER PmeOutput pme_gpu_wait_finish_task(gmx_pme_t*     GPU_FUNC_ARGUMENT(pme),
                                                      int            GPU_FUNC_ARGUMENT(flags),
                                                      gmx_wallcycle* GPU_FUNC_ARGUMENT(wcycle))
        GPU_FUNC_TERM_WITH_RETURN(PmeOutput{});

#endif
