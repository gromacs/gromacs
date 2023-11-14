/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2016- The GROMACS Authors
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

#include "gromacs/fft/fft.h" // for the gmx_fft_direction enum
#include "gromacs/gpu_utils/devicebuffer_datatype.h"
#include "gromacs/gpu_utils/gpu_macros.h" // for the GPU_FUNC_ macros

#include "pme_gpu_types_host.h"
#include "pme_output.h"

class DeviceContext;
struct DeviceInformation;
class DeviceStream;
class GpuEventSynchronizer;
struct gmx_hw_info_t;
struct gmx_pme_t; // only used in pme_gpu_reinit
struct gmx_wallcycle;
class PmeAtomComm;
enum class PmeForceOutputHandling;
struct PmeGpu;
class PmeGpuProgram;
struct PmeGpuStaging;
struct PmeGpuSettings;
struct t_complex;
typedef struct gmx_parallel_3dfft* gmx_parallel_3dfft_t;

#ifndef FEP_STATE_A
//! Grid index of FEP state A (or unperturbed system)
#    define FEP_STATE_A 0
#endif
#ifndef FEP_STATE_B
//! Grid index of FEP state B
#    define FEP_STATE_B 1
#endif

namespace gmx
{
template<typename>
class ArrayRef;
class MDLogger;
} // namespace gmx

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
    XYZ,
    Count
};

/*! \libinternal \brief
 * Returns the size of the block size requirement
 *
 * The GPU version of PME requires that the coordinates array have a
 * size divisible by the returned number.
 *
 * \returns Number of atoms in a single GPU atom data chunk, which
 * determines a minimum divisor of the size of the memory allocated.
 */
int pme_gpu_get_atom_data_block_size();

/*!\brief Return the number of atoms per warp */
GPU_FUNC_QUALIFIER int pme_gpu_get_atoms_per_warp(const PmeGpu* GPU_FUNC_ARGUMENT(pmeGpu))
        GPU_FUNC_TERM_WITH_RETURN(0);

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
 * \param[in] pmeGpu                          The PME GPU structure.
 * \param[in] gpuGraphWithSeparatePmeRank     Whether MD GPU Graph with separate PME rank is in use.
 */
void pme_gpu_clear_energy_virial(const PmeGpu* pmeGpu, bool gpuGraphWithSeparatePmeRank);

/*! \libinternal \brief
 * Reallocates and copies the pre-computed B-spline values to the GPU.
 *
 * \param[in,out] pmeGpu             The PME GPU structure.
 * \param[in]     gridIndex          The index of the grid to use. 0 is Coulomb in the normal
 *                                   state or FEP state A and 1 is Coulomb in FEP state B.
 */
void pme_gpu_realloc_and_copy_bspline_values(PmeGpu* pmeGpu, int gridIndex = 0);

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
 * Reallocates the buffer on the GPU and copies the charges/coefficients from the CPU buffer.
 * Clears the padded part if needed.
 *
 * \param[in] pmeGpu            The PME GPU structure.
 * \param[in] h_coefficients    The input atom charges/coefficients.
 * \param[in] gridIndex         The index of the grid to use. 0 is Coulomb in the normal
 *                              state or FEP state A and 1 is Coulomb in FEP state B.
 *
 * Does not need to be done for every PME computation, only whenever the local charges change.
 * (So, in the beginning of the run, or on DD step).
 */
void pme_gpu_realloc_and_copy_input_coefficients(const PmeGpu* pmeGpu,
                                                 const float*  h_coefficients,
                                                 int           gridIndex = 0);

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
 * Reinitialize PME halo exchange parameters and staging device buffers for MPI communication.
 *
 * \param[in] pmeGpu            The PME GPU structure.
 */
void pme_gpu_reinit_haloexchange(PmeGpu* pmeGpu);

/*! \libinternal \brief
 * Frees device staging buffers used for PME halo exchange.
 *
 * \param[in] pmeGpu            The PME GPU structure.
 */
void pme_gpu_free_haloexchange(const PmeGpu* pmeGpu);

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
 * \param[in] pmeGpu    The PME GPU structure.
 * \param[in] h_grid    The host-side grid buffer.
 * \param[in] gridIndex The index of the grid to use. 0 is Coulomb in the normal
 *                      state or FEP state A and 1 is Coulomb in FEP state B.
 */
void pme_gpu_copy_input_gather_grid(const PmeGpu* pmeGpu, const float* h_grid, int gridIndex = 0);

/*! \libinternal \brief
 * Copies the output real-space grid from the GPU to the host.
 *
 * \param[in] pmeGpu    The PME GPU structure.
 * \param[out] h_grid   The host-side grid buffer.
 * \param[in] gridIndex The index of the grid to use. 0 is Coulomb in the normal
 *                      state or FEP state A and 1 is Coulomb in FEP state B.
 */
void pme_gpu_copy_output_spread_grid(const PmeGpu* pmeGpu, float* h_grid, int gridIndex = 0);

/*! \libinternal \brief
 * Copies the spread output spline data and gridline indices from the GPU to the host.
 *
 * \param[in] pmeGpu    The PME GPU structure.
 */
void pme_gpu_copy_output_spread_atom_data(const PmeGpu* pmeGpu);

/*! \libinternal \brief
 * Copies the gather input spline data and gridline indices from the host to the GPU.
 *
 * \param[in] pmeGpu    The PME GPU structure.
 */
void pme_gpu_copy_input_gather_atom_data(const PmeGpu* pmeGpu);

/*! \libinternal \brief
 * Waits for the grid copying to the host-side buffer after spreading to finish.
 *
 * \param[in] pmeGpu  The PME GPU structure.
 */
void pme_gpu_sync_spread_grid(const PmeGpu* pmeGpu);

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

/* The PME stages themselves */

/*! \libinternal \brief
 * A GPU spline computation and charge spreading function.
 *
 * \param[in]  pmeGpu                    The PME GPU structure.
 * \param[in]  xReadyOnDevice            Event synchronizer indicating that the coordinates are
 *                                       ready in the device memory; can be nullptr when invoked
 *                                       on a separate PME rank or from PME tests.
 * \param[out] h_grids                   The host-side grid buffers (used only if the result
 *                                       of the spread is expected on the host, e.g. testing
 *                                       or host-side FFT)
 * \param[in]  fftSetup                  Host-side FFT setup structure used in Mixed mode
 * \param[in]  computeSplines            Should the computation of spline parameters and gridline
 *                                       indices be performed.
 * \param[in]  spreadCharges             Should the charges/coefficients be spread on the grid.
 * \param[in]  lambda                    The lambda value of the current system state.
 * \param[in]  useGpuDirectComm          Whether direct GPU PME-PP communication is active
 * \param[in]  pmeCoordinateReceiverGpu  Coordinate receiver object, which must be valid when
 *                                       direct GPU PME-PP communication is active
 * \param[in]  useMdGpuGraph             Whether MD GPU Graph is in use.
 * \param[in]  wcycle                    The wallclock counter.
 */
GPU_FUNC_QUALIFIER void pme_gpu_spread(const PmeGpu*         GPU_FUNC_ARGUMENT(pmeGpu),
                                       GpuEventSynchronizer* GPU_FUNC_ARGUMENT(xReadyOnDevice),
                                       float**               GPU_FUNC_ARGUMENT(h_grids),
                                       gmx_parallel_3dfft_t* GPU_FUNC_ARGUMENT(fftSetup),
                                       bool                  GPU_FUNC_ARGUMENT(computeSplines),
                                       bool                  GPU_FUNC_ARGUMENT(spreadCharges),
                                       real                  GPU_FUNC_ARGUMENT(lambda),
                                       bool                  GPU_FUNC_ARGUMENT(useGpuDirectComm),
                                       gmx::PmeCoordinateReceiverGpu* GPU_FUNC_ARGUMENT(pmeCoordinateReceiverGpu),
                                       bool           GPU_FUNC_ARGUMENT(useMdGpuGraph),
                                       gmx_wallcycle* GPU_FUNC_ARGUMENT(wcycle)) GPU_FUNC_TERM;

/*! \libinternal \brief
 * 3D FFT R2C/C2R routine.
 *
 * \param[in]  pmeGpu          The PME GPU structure.
 * \param[in]  direction       Transform direction (real-to-complex or complex-to-real)
 * \param[in]  gridIndex       The index of the grid to use. 0 is Coulomb in the normal
 *                             state or FEP state A and 1 is Coulomb in FEP state B.
 */
void pme_gpu_3dfft(const PmeGpu* pmeGpu, enum gmx_fft_direction direction, int gridIndex = 0);

/*! \libinternal \brief
 * A GPU Fourier space solving function.
 *
 * \param[in]     pmeGpu                  The PME GPU structure.
 * \param[in]     gridIndex               The index of the grid to use. 0 is Coulomb in the normal
 *                                        state or FEP state A and 1 is Coulomb in FEP state B.
 * \param[in,out] h_grid                  The host-side input and output Fourier grid buffer (used only with testing or host-side FFT)
 * \param[in]     gridOrdering            Specifies the dimenion ordering of the complex grid. TODO: store this information?
 * \param[in]     computeEnergyAndVirial  Tells if the energy and virial computation should be performed.
 */
GPU_FUNC_QUALIFIER void pme_gpu_solve(const PmeGpu* GPU_FUNC_ARGUMENT(pmeGpu),
                                      int           GPU_FUNC_ARGUMENT(gridIndex),
                                      t_complex*    GPU_FUNC_ARGUMENT(h_grid),
                                      GridOrdering  GPU_FUNC_ARGUMENT(gridOrdering),
                                      bool GPU_FUNC_ARGUMENT(computeEnergyAndVirial)) GPU_FUNC_TERM;

/*! \libinternal \brief
 * A GPU force gathering function.
 *
 * \param[in]     pmeGpu                   The PME GPU structure.
 * \param[in]     h_grids                  The host-side grid buffer (used only in testing mode).
 * \param[in]     fftSetup                 Host-side FFT setup structure used in Mixed mode
 * \param[in]     lambda                   The lambda value to use.
 * \param[in]     wcycle                   The wallclock counter.
 * \param[in]     computeVirial            Whether this is a virial step.
 */
GPU_FUNC_QUALIFIER void pme_gpu_gather(PmeGpu*               GPU_FUNC_ARGUMENT(pmeGpu),
                                       float**               GPU_FUNC_ARGUMENT(h_grids),
                                       gmx_parallel_3dfft_t* GPU_FUNC_ARGUMENT(fftSetup),
                                       float                 GPU_FUNC_ARGUMENT(lambda),
                                       gmx_wallcycle*        GPU_FUNC_ARGUMENT(wcycle),
                                       bool GPU_FUNC_ARGUMENT(computeVirial)) GPU_FUNC_TERM;


/*! \brief Sets the device pointer to coordinate data
 * \param[in] pmeGpu         The PME GPU structure.
 * \param[in] d_x            Pointer to coordinate data
 */
GPU_FUNC_QUALIFIER void pme_gpu_set_kernelparam_coordinates(const PmeGpu* GPU_FUNC_ARGUMENT(pmeGpu),
                                                            DeviceBuffer<gmx::RVec> GPU_FUNC_ARGUMENT(d_x)) GPU_FUNC_TERM;

/*! \brief Return pointer to device copy of force data.
 * \param[in] pmeGpu         The PME GPU structure.
 * \returns                  Pointer to force data
 */
GPU_FUNC_QUALIFIER DeviceBuffer<gmx::RVec> pme_gpu_get_kernelparam_forces(const PmeGpu* GPU_FUNC_ARGUMENT(pmeGpu))
        GPU_FUNC_TERM_WITH_RETURN(DeviceBuffer<gmx::RVec>{});

GPU_FUNC_QUALIFIER void pme_gpu_set_kernelparam_useNvshmem(const PmeGpu* GPU_FUNC_ARGUMENT(pmeGpu),
                                                           bool GPU_FUNC_ARGUMENT(useNvshmem)) GPU_FUNC_TERM;

/*! \brief Return pointer to the sync object triggered after the PME force calculation completion
 * \param[in] pmeGpu         The PME GPU structure.
 * \returns                  Pointer to sync object
 */
GPU_FUNC_QUALIFIER GpuEventSynchronizer* pme_gpu_get_forces_ready_synchronizer(
        const PmeGpu* GPU_FUNC_ARGUMENT(pmeGpu)) GPU_FUNC_TERM_WITH_RETURN(nullptr);

/*! \libinternal \brief
 * Returns the PME GPU settings
 *
 * \param[in] pmeGpu         The PME GPU structure.
 * \returns                  The settings for PME on GPU
 */
inline const PmeGpuSettings& pme_gpu_settings(const PmeGpu* pmeGpu)
{
    return pmeGpu->settings;
}

/*! \libinternal \brief
 * Returns the PME GPU staging object
 *
 * \param[in] pmeGpu         The PME GPU structure.
 * \returns                  The staging object for PME on GPU
 */
inline const PmeGpuStaging& pme_gpu_staging(const PmeGpu* pmeGpu)
{
    return pmeGpu->staging;
}

/*! \libinternal \brief
 * Sets whether the PME module is running in testing mode
 *
 * \param[in] pmeGpu         The PME GPU structure.
 * \param[in] testing        Whether testing mode is on.
 */
inline void pme_gpu_set_testing(PmeGpu* pmeGpu, bool testing)
{
    if (pmeGpu)
    {
        pmeGpu->settings.copyAllOutputs = testing;
        pmeGpu->settings.transferKind = testing ? GpuApiCallBehavior::Sync : GpuApiCallBehavior::Async;
    }
}

/* A block of C++ functions that live in pme_gpu_internal.cpp */

/*! \libinternal \brief
 * Returns the energy and virial GPU outputs, useful for testing.
 *
 * It is the caller's responsibility to be aware of whether the GPU
 * handled the solve stage.
 *
 * \param[in] pme                The PME structure.
 * \param[in] lambda             The lambda value to use when calculating the results.
 * \param[out] output            Pointer to output where energy and virial should be stored.
 */
GPU_FUNC_QUALIFIER void pme_gpu_getEnergyAndVirial(const gmx_pme_t& GPU_FUNC_ARGUMENT(pme),
                                                   float            GPU_FUNC_ARGUMENT(lambda),
                                                   PmeOutput* GPU_FUNC_ARGUMENT(output)) GPU_FUNC_TERM;

/*! \libinternal \brief
 * Returns the GPU outputs (forces, energy and virial)
 *
 * \param[in] pme                     The PME structure.
 * \param[in] computeEnergyAndVirial  Whether the energy and virial are being computed
 * \param[in] lambdaQ            The Coulomb lambda to use when finalizing the output.
 * \returns                           The output object.
 */
GPU_FUNC_QUALIFIER PmeOutput pme_gpu_getOutput(const gmx_pme_t& GPU_FUNC_ARGUMENT(pme),
                                               bool GPU_FUNC_ARGUMENT(computeEnergyAndVirial),
                                               real GPU_FUNC_ARGUMENT(lambdaQ))
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
 * \param[in,out] pme               The PME structure.
 * \param[in]     deviceContext     The GPU context.
 * \param[in]     deviceStream      The GPU stream.
 * \param[in,out] pmeGpuProgram     The handle to the program/kernel data created outside (e.g. in unit tests/runner)
 * \param[in]     useMdGpuGraph     Whether MD GPU Graph is in use
 * \throws gmx::NotImplementedError if this generally valid PME structure is not valid for GPU runs.
 */
GPU_FUNC_QUALIFIER void pme_gpu_reinit(gmx_pme_t*           GPU_FUNC_ARGUMENT(pme),
                                       const DeviceContext* GPU_FUNC_ARGUMENT(deviceContext),
                                       const DeviceStream*  GPU_FUNC_ARGUMENT(deviceStream),
                                       const PmeGpuProgram* GPU_FUNC_ARGUMENT(pmeGpuProgram),
                                       bool GPU_FUNC_ARGUMENT(useMdGpuGraph)) GPU_FUNC_TERM;

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
 * \param[in] chargesA  The pointer to the host-side array of particle charges in the unperturbed state or FEP state A.
 * \param[in] chargesB  The pointer to the host-side array of particle charges in FEP state B.
 *
 * This is a function that should only be called in the beginning of the run and on domain
 * decomposition. Should be called before the pme_gpu_set_io_ranges.
 */
GPU_FUNC_QUALIFIER void pme_gpu_reinit_atoms(PmeGpu*     GPU_FUNC_ARGUMENT(pmeGpu),
                                             int         GPU_FUNC_ARGUMENT(nAtoms),
                                             const real* GPU_FUNC_ARGUMENT(chargesA),
                                             const real* GPU_FUNC_ARGUMENT(chargesB) = nullptr) GPU_FUNC_TERM;

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
 * \param[in]  pme                     The PME data structure.
 * \param[in]  computeEnergyAndVirial  Tells if the energy and virial computation should be performed.
 * \param[in]  lambdaQ                 The Coulomb lambda to use when calculating the results.
 * \param[out] wcycle                  The wallclock counter.
 * \return                             The output forces, energy and virial
 */
GPU_FUNC_QUALIFIER PmeOutput pme_gpu_wait_finish_task(gmx_pme_t* GPU_FUNC_ARGUMENT(pme),
                                                      bool GPU_FUNC_ARGUMENT(computeEnergyAndVirial),
                                                      real           GPU_FUNC_ARGUMENT(lambdaQ),
                                                      gmx_wallcycle* GPU_FUNC_ARGUMENT(wcycle))
        GPU_FUNC_TERM_WITH_RETURN(PmeOutput{});

#endif
