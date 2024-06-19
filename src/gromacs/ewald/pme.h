/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
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
/*! \libinternal \file
 *
 * \brief This file contains function declarations necessary for
 * computing energies and forces for the PME long-ranged part (Coulomb
 * and LJ).
 *
 * \author Berk Hess <hess@kth.se>
 * \inlibraryapi
 * \ingroup module_ewald
 */

#ifndef GMX_EWALD_PME_H
#define GMX_EWALD_PME_H

#include <memory>
#include <string>
#include <vector>

#include "gromacs/gpu_utils/devicebuffer_datatype.h"
#include "gromacs/gpu_utils/gpu_macros.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

struct gmx_hw_info_t;
struct t_commrec;
struct t_inputrec;
struct t_nrnb;
struct PmeGpu;
struct gmx_wallclock_gpu_pme_t;
struct gmx_enerdata_t;
struct gmx_mtop_t;
struct gmx_pme_t;
struct gmx_wallcycle;
struct NumPmeDomains;
struct PmeGridsStorage;

class DeviceContext;
class DeviceStream;
enum class GpuTaskCompletion;
class PmeGpuProgram;
class GpuEventSynchronizer;

namespace gmx
{
template<typename>
class ArrayRef;
class ForceWithVirial;
class MDLogger;
enum class PinningPolicy : int;
class StepWorkload;

/*! \libinternal \brief Class for managing usage of separate PME-only ranks
 *
 * Used for checking if some parts of the code could not use PME-only ranks
 *
 */
class SeparatePmeRanksPermitted
{
public:
    //! Disables PME ranks permitted flag with a reason
    void disablePmeRanks(const std::string& reason);
    //! Return status of PME ranks usage
    bool permitSeparatePmeRanks() const;
    //! Returns all reasons, for not using PME ranks
    std::string reasonsWhyDisabled() const;

private:
    //! Flag that informs whether simualtion could use dedicated PME ranks
    bool permitSeparatePmeRanks_ = true;
    //! Storage for all reasons, why PME ranks could not be used
    std::vector<std::string> reasons_;
};

class PmeCoordinateReceiverGpu;
} // namespace gmx

enum
{
    GMX_SUM_GRID_FORWARD,
    GMX_SUM_GRID_BACKWARD
};

/*! \brief Possible PME codepaths on a rank.
 * \todo: make this enum class with gmx_pme_t C++ refactoring
 */
enum class PmeRunMode
{
    None,  //!< No PME task is done
    CPU,   //!< Whole PME computation is done on CPU
    GPU,   //!< Whole PME computation is done on GPU
    Mixed, //!< Mixed mode: only spread and gather run on GPU; FFT and solving are done on CPU.
};

/*! \brief Return the smallest allowed PME grid size for \p pmeOrder */
int minimalPmeGridSize(int pmeOrder);

/*! \brief Calculates the number of grid lines used as halo region for PME decomposition.
 *
 * This function is only used for PME GPU decomposition case
 *
 * \param[in] pmeOrder                      PME interpolation order
 * \param[in] haloExtentForAtomDisplacement extent of halo region in nm to account for atom
 *                                          displacement
 * \param[in] gridSpacing                   spacing between grid lines
 *
 * \returns  extended halo region used in GPU PME decomposition
 */
int numGridLinesForExtendedHaloRegion(int pmeOrder, real haloExtentForAtomDisplacement, real gridSpacing);

/*! \brief
 * Gets grid spacing from simulation box and grid dimension
 *
 * \param[in]     box       Simulation box
 * \param[in]     gridDim   Fourier grid dimension
 *
 * \returns  Grid spacing value
 */
real getGridSpacingFromBox(const matrix box, const ivec gridDim);

//! Return whether the grid of \c pme is identical to \c grid_size.
bool gmx_pme_grid_matches(const gmx_pme_t& pme, const ivec grid_size);

/*! \brief Check restrictions on pme_order and the PME grid nkx,nky,nkz.
 *
 * With errorsAreFatal=true, an exception or fatal error is generated
 * on violation of restrictions.
 * With errorsAreFatal=false, false is returned on violation of restrictions.
 * When all restrictions are obeyed, true is returned.
 * Argument useThreads tells if any MPI rank doing PME uses more than 1 threads.
 * If at calling useThreads is unknown, pass true for conservative checking.
 *
 * The PME GPU restrictions are checked separately during pme_gpu_init().
 */
bool gmx_pme_check_restrictions(int  pme_order,
                                int  nkx,
                                int  nky,
                                int  nkz,
                                int  numPmeDomainsAlongX,
                                int  numPmeDomainsAlongY,
                                int  extendedHaloRegion,
                                bool useGpuPme,
                                bool useThreads,
                                bool errorsAreFatal);

/*! \brief Construct PME data
 *
 * \throws   gmx::InconsistentInputError if input grid sizes/PME order are inconsistent.
 * \returns  Pointer to newly allocated and initialized PME data.
 *
 * \p pmeGridsStorage can be nullptr, in which case new PmeGridsStorage is allocated.
 * When \p pmeGridsStorage is not a nullptr, grid storage is taken from there.
 *
 * \todo We should evolve something like a \c GpuManager that holds \c
 * DeviceInformation* and \c PmeGpuProgram* and perhaps other
 * related things whose lifetime can/should exceed that of a task (or
 * perhaps task manager). See Issue #2522.
 */
gmx_pme_t* gmx_pme_init(const t_commrec*                 cr,
                        const NumPmeDomains&             numPmeDomains,
                        const t_inputrec*                ir,
                        const matrix                     box,
                        real                             haloExtentForAtomDisplacement,
                        gmx_bool                         bFreeEnergy_q,
                        gmx_bool                         bFreeEnergy_lj,
                        gmx_bool                         bReproducible,
                        real                             ewaldcoeff_q,
                        real                             ewaldcoeff_lj,
                        int                              nthread,
                        PmeRunMode                       runMode,
                        PmeGpu*                          pmeGpu,
                        const DeviceContext*             deviceContext,
                        const DeviceStream*              deviceStream,
                        const PmeGpuProgram*             pmeGpuProgram,
                        const gmx::MDLogger&             mdlog,
                        std::shared_ptr<PmeGridsStorage> pmeGridsStoragePtr);

/*! \brief As gmx_pme_init, but takes most settings, except the grid/Ewald coefficients,
 * and the shared grid storage from pme_src.
 */
void gmx_pme_reinit(gmx_pme_t**       pmedata,
                    const t_commrec*  cr,
                    gmx_pme_t*        pme_src,
                    const t_inputrec* ir,
                    const ivec        grid_size,
                    real              ewaldcoeff_q,
                    real              ewaldcoeff_lj);

/*! \brief Destroys the PME data structure (including GPU data). */
void gmx_pme_destroy(gmx_pme_t* pme);

/*! \brief Destroys the PME data structure.
 *
 * \param pme             The data structure to destroy.
 * \param destroyGpuData  Set to \c false if \p pme is a copy created by \ref gmx_pme_reinit,
 * and only it should be destroyed, while shared GPU data should be preserved. If \c true,
 * the shared data will be destroyed, like in \c gmx_pme_destroy(gmx_pme_t* pme).
 */
void gmx_pme_destroy(gmx_pme_t* pme, bool destroyGpuData);

/*! \brief Do a PME calculation on a CPU for the long range electrostatics and/or LJ.
 *
 * Computes the PME forces and the energy and viral, when requested,
 * for all atoms in \p coordinates. Forces, when requested, are added
 * to the buffer \p forces, which is allowed to contain more elements
 * than the number of elements in \p coordinates.
 * The meaning of \p flags is defined above, and determines which
 * parts of the calculation are performed.
 *
 * \return 0 indicates all well, non zero is an error code.
 */
int gmx_pme_do(struct gmx_pme_t*              pme,
               gmx::ArrayRef<const gmx::RVec> coordinates,
               gmx::ArrayRef<gmx::RVec>       forces,
               gmx::ArrayRef<const real>      chargeA,
               gmx::ArrayRef<const real>      chargeB,
               gmx::ArrayRef<const real>      c6A,
               gmx::ArrayRef<const real>      c6B,
               gmx::ArrayRef<const real>      sigmaA,
               gmx::ArrayRef<const real>      sigmaB,
               const matrix                   box,
               const t_commrec*               cr,
               int                            maxshift_x,
               int                            maxshift_y,
               t_nrnb*                        nrnb,
               gmx_wallcycle*                 wcycle,
               matrix                         vir_q,
               matrix                         vir_lj,
               real*                          energy_q,
               real*                          energy_lj,
               real                           lambda_q,
               real                           lambda_lj,
               real*                          dvdlambda_q,
               real*                          dvdlambda_lj,
               const gmx::StepWorkload&       stepWork);

/*! \brief Calculate the PME grid energy V for n charges.
 *
 * The potential (found in \p pme) must have been found already with a
 * call to gmx_pme_do(). Note that the charges are not spread on the grid in the
 * pme struct. Currently does not work in parallel or with free
 * energy.
 */
real gmx_pme_calc_energy(gmx_pme_t* pme, gmx::ArrayRef<const gmx::RVec> x, gmx::ArrayRef<const real> q);

/*! \brief
 * This function updates the local atom data on GPU after DD (charges, coordinates, etc.).
 * TODO: it should update the PME CPU atom data as well.
 * (currently PME CPU call gmx_pme_do() gets passed the input pointers for each computation).
 *
 * \param[in,out] pme        The PME structure.
 * \param[in]     numAtoms   The number of particles.
 * \param[in]     chargesA   The pointer to the array of particle charges in the normal state or FEP
 * state A. Can be nullptr if PME is not performed on the GPU.
 * \param[in]     chargesB   The pointer to the array of particle charges in state B. Only used if
 * charges are perturbed and can otherwise be nullptr.
 */
void gmx_pme_reinit_atoms(gmx_pme_t*                pme,
                          int                       numAtoms,
                          gmx::ArrayRef<const real> chargesA,
                          gmx::ArrayRef<const real> chargesB);

/* A block of PME GPU functions */

/*! \brief Checks whether the GROMACS build allows to run PME on GPU.
 * TODO: this partly duplicates an internal PME assert function
 * pme_gpu_check_restrictions(), except that works with a
 * formed gmx_pme_t structure. Should that one go away/work with inputrec?
 *
 * \param[out] error   If non-null, the error message when PME is not supported on GPU.
 *
 * \returns true if PME can run on GPU on this build, false otherwise.
 */
bool pme_gpu_supports_build(std::string* error);

/*! \brief Checks whether the input system allows to run PME on GPU.
 * TODO: this partly duplicates an internal PME assert function
 * pme_gpu_check_restrictions(), except that works with a
 * formed gmx_pme_t structure. Should that one go away/work with inputrec?
 *
 * \param[in]  ir     Input system.
 * \param[out] error  If non-null, the error message if the input is not supported on GPU.
 *
 * \returns true if PME can run on GPU with this input, false otherwise.
 */
bool pme_gpu_supports_input(const t_inputrec& ir, std::string* error);

/*! \brief Checks whether the input system allows to run PME on GPU in Mixed mode.
 * Assumes that the input system is compatible with GPU PME otherwise, that is,
 * before calling this function one should check that \ref pme_gpu_supports_input returns \c true.
 *
 * \param[in]  ir     Input system.
 * \param[out] error  If non-null, the error message if the input is not supported.
 *
 * \returns true if PME can run on GPU in Mixed mode with this input, false otherwise.
 */
bool pme_gpu_mixed_mode_supports_input(const t_inputrec& ir, std::string* error);

/*! \brief
 * Returns the active PME codepath (CPU, GPU, mixed).
 * \todo This is a rather static data that should be managed by the higher level task scheduler.
 *
 * \param[in]  pme            The PME data structure.
 * \returns active PME codepath.
 */
PmeRunMode pme_run_mode(const gmx_pme_t* pme);

/*! \libinternal \brief
 * Return the pinning policy appropriate for this build configuration
 * for relevant buffers used for PME task on this rank (e.g. running
 * on a GPU). */
gmx::PinningPolicy pme_get_pinning_policy();

/*! \brief
 * Tells if PME is enabled to run on GPU (not necessarily active at the moment).
 * \todo This is a rather static data that should be managed by the hardware assignment manager.
 * For now, it is synonymous with the active PME codepath (in the absence of dynamic switching).
 *
 * \param[in]  pme            The PME data structure.
 * \returns true if PME can run on GPU, false otherwise.
 */
inline bool pme_gpu_task_enabled(const gmx_pme_t* pme)
{
    return (pme != nullptr) && (pme_run_mode(pme) != PmeRunMode::CPU);
}

/*! \libinternal \brief
 * Sets the nvshmem usage status, and allocates required structs
 * if NVSHMEM should be used.
 *
 * \param[in] pmeGpu             The PME GPU structure.
 * \param[in] useNvshmem         should use NVSHMEM.
 */
GPU_FUNC_QUALIFIER void pme_gpu_use_nvshmem(PmeGpu* GPU_FUNC_ARGUMENT(pmeGpu),
                                            bool    GPU_FUNC_ARGUMENT(useNvshmem)) GPU_FUNC_TERM;

/*! \brief Returns the block size requirement
 *
 * The GPU version of PME requires that the coordinates array have a
 * size divisible by the returned number.
 *
 * \param[in]  pme  The PME data structure.
 */
GPU_FUNC_QUALIFIER int pme_gpu_get_block_size(const gmx_pme_t* GPU_FUNC_ARGUMENT(pme))
        GPU_FUNC_TERM_WITH_RETURN(0);

// The following functions are all the PME GPU entry points,
// currently inlining to nothing on non-CUDA builds.

/*! \brief
 * Resets the PME GPU timings. To be called at the reset step.
 *
 * \param[in] pme            The PME structure.
 */
GPU_FUNC_QUALIFIER void pme_gpu_reset_timings(const gmx_pme_t* GPU_FUNC_ARGUMENT(pme)) GPU_FUNC_TERM;

/*! \brief
 * Copies the PME GPU timings to the gmx_wallclock_gpu_pme_t structure (for log output). To be called at the run end.
 *
 * \param[in] pme               The PME structure.
 * \param[in] timings           The gmx_wallclock_gpu_pme_t structure.
 */
GPU_FUNC_QUALIFIER void pme_gpu_get_timings(const gmx_pme_t* GPU_FUNC_ARGUMENT(pme),
                                            gmx_wallclock_gpu_pme_t* GPU_FUNC_ARGUMENT(timings)) GPU_FUNC_TERM;

/* The main PME GPU functions */

/*! \brief
 * Prepares PME on GPU computation (updating the box if needed)
 * \param[in] pme               The PME data structure.
 * \param[in] box               The unit cell box.
 * \param[in] wcycle            The wallclock counter.
 * \param[in] stepWork          The required work for this simulation step
 */
GPU_FUNC_QUALIFIER void pme_gpu_prepare_computation(gmx_pme_t*     GPU_FUNC_ARGUMENT(pme),
                                                    const matrix   GPU_FUNC_ARGUMENT(box),
                                                    gmx_wallcycle* GPU_FUNC_ARGUMENT(wcycle),
                                                    const gmx::StepWorkload& GPU_FUNC_ARGUMENT(stepWork)) GPU_FUNC_TERM;

/*! \brief
 * Launches first stage of PME on GPU - spreading kernel.
 *
 * \param[in] pme                            The PME data structure.
 * \param[in] xReadyOnDevice                 Event synchronizer indicating that the coordinates
 *                                           are ready in the device memory; nullptr allowed only
 *                                           on separate PME ranks.
 * \param[in] wcycle                         The wallclock counter.
 * \param[in] lambdaQ                        The Coulomb lambda of the current state of the
 *                                           system. Only used if FEP of Coulomb is active.
 * \param[in] useGpuDirectComm               Whether direct GPU PME-PP communication is active
 * \param[in]  pmeCoordinateReceiverGpu      Coordinate receiver object, which must be valid when
 *                                           direct GPU PME-PP communication is active
 * \param[in] useMdGpuGraph                  Whether MD GPU Graph is in use.
 */
GPU_FUNC_QUALIFIER void pme_gpu_launch_spread(gmx_pme_t*            GPU_FUNC_ARGUMENT(pme),
                                              GpuEventSynchronizer* GPU_FUNC_ARGUMENT(xReadyOnDevice),
                                              gmx_wallcycle*        GPU_FUNC_ARGUMENT(wcycle),
                                              real                  GPU_FUNC_ARGUMENT(lambdaQ),
                                              bool GPU_FUNC_ARGUMENT(useGpuDirectComm),
                                              gmx::PmeCoordinateReceiverGpu* GPU_FUNC_ARGUMENT(pmeCoordinateReceiverGpu),
                                              bool GPU_FUNC_ARGUMENT(useMdGpuGraph)) GPU_FUNC_TERM;

/*! \brief
 * Launches middle stages of PME (FFT R2C, solving, FFT C2R) either on GPU or on CPU, depending on the run mode.
 *
 * \param[in] pme               The PME data structure.
 * \param[in] wcycle            The wallclock counter.
 * \param[in] stepWork          The required work for this simulation step
 */
GPU_FUNC_QUALIFIER void
pme_gpu_launch_complex_transforms(gmx_pme_t*               GPU_FUNC_ARGUMENT(pme),
                                  gmx_wallcycle*           GPU_FUNC_ARGUMENT(wcycle),
                                  const gmx::StepWorkload& GPU_FUNC_ARGUMENT(stepWork)) GPU_FUNC_TERM;

/*! \brief
 * Launches last stage of PME on GPU - force gathering and D2H force transfer.
 *
 * \param[in] pme               The PME data structure.
 * \param[in] wcycle            The wallclock counter.
 * \param[in] lambdaQ           The Coulomb lambda to use when calculating the results.
 * \param[in] computeVirial     Whether this is a virial step.
 */
GPU_FUNC_QUALIFIER void pme_gpu_launch_gather(gmx_pme_t*     GPU_FUNC_ARGUMENT(pme),
                                              gmx_wallcycle* GPU_FUNC_ARGUMENT(wcycle),
                                              real           GPU_FUNC_ARGUMENT(lambdaQ),
                                              bool GPU_FUNC_ARGUMENT(computeVirial)) GPU_FUNC_TERM;

/*! \brief
 * Attempts to complete PME GPU tasks.
 *
 * The \p completionKind argument controls whether the function blocks until all
 * PME GPU tasks enqueued completed (as pme_gpu_wait_finish_task() does) or only
 * checks and returns immediately if they did not.
 * When blocking or the tasks have completed it also gets the output forces
 * by assigning the ArrayRef to the \p forces pointer passed in.
 * Virial/energy are also outputs if they were to be computed.
 *
 * \param[in]  pme             The PME data structure.
 * \param[in]  stepWork        The required work for this simulation step
 * \param[in]  wcycle          The wallclock counter.
 * \param[out] forceWithVirial The output force and virial
 * \param[out] enerd           The output energies
 * \param[in]  lambdaQ         The Coulomb lambda to use when calculating the results.
 * \param[in]  completionKind  Indicates whether PME task completion should only be checked rather
 *                             than waited for
 * \returns                    True if the PME GPU tasks have completed
 */
GPU_FUNC_QUALIFIER bool pme_gpu_try_finish_task(gmx_pme_t*               GPU_FUNC_ARGUMENT(pme),
                                                const gmx::StepWorkload& GPU_FUNC_ARGUMENT(stepWork),
                                                gmx_wallcycle*           GPU_FUNC_ARGUMENT(wcycle),
                                                gmx::ForceWithVirial* GPU_FUNC_ARGUMENT(forceWithVirial),
                                                gmx_enerdata_t*       GPU_FUNC_ARGUMENT(enerd),
                                                real                  GPU_FUNC_ARGUMENT(lambdaQ),
                                                GpuTaskCompletion GPU_FUNC_ARGUMENT(completionKind))
        GPU_FUNC_TERM_WITH_RETURN(false);

/*! \brief
 * Blocks until PME GPU tasks are completed, and gets the output forces and virial/energy
 * (if they were to be computed).
 *
 * \param[in]  pme             The PME data structure.
 * \param[in]  stepWork        The required work for this simulation step
 * \param[in]  wcycle          The wallclock counter.
 * \param[out] forceWithVirial The output force and virial
 * \param[out] enerd           The output energies
 * \param[in]  lambdaQ         The Coulomb lambda to use when calculating the results.
 */
GPU_FUNC_QUALIFIER void pme_gpu_wait_and_reduce(gmx_pme_t*               GPU_FUNC_ARGUMENT(pme),
                                                const gmx::StepWorkload& GPU_FUNC_ARGUMENT(stepWork),
                                                gmx_wallcycle*           GPU_FUNC_ARGUMENT(wcycle),
                                                gmx::ForceWithVirial* GPU_FUNC_ARGUMENT(forceWithVirial),
                                                gmx_enerdata_t*       GPU_FUNC_ARGUMENT(enerd),
                                                real GPU_FUNC_ARGUMENT(lambdaQ)) GPU_FUNC_TERM;

/*! \brief
 * The PME GPU reinitialization function that is called both at the end of any PME computation and on any load balancing.
 *
 * Clears the internal grid and energy/virial buffers; it is not safe to start
 * the PME computation without calling this.
 * Note that unlike in the nbnxn module, the force buffer does not need clearing.
 *
 * \todo Rename this function to *clear* -- it clearly only does output resetting
 * and we should be clear about what the function does..
 *
 * \param[in] pme                            The PME data structure.
 * \param[in] gpuGraphWithSeparatePmeRank    Whether MD GPU Graph with separate PME rank is in use.
 * \param[in] wcycle                         The wallclock counter.
 */
GPU_FUNC_QUALIFIER void pme_gpu_reinit_computation(const gmx_pme_t* GPU_FUNC_ARGUMENT(pme),
                                                   bool GPU_FUNC_ARGUMENT(gpuGraphWithSeparatePmeRank),
                                                   gmx_wallcycle* GPU_FUNC_ARGUMENT(wcycle)) GPU_FUNC_TERM;

/*! \brief Set pointer to device copy of coordinate data.
 * \param[in] pme            The PME data structure.
 * \param[in] d_x            The pointer to the positions buffer to be set
 */
GPU_FUNC_QUALIFIER void pme_gpu_set_device_x(const gmx_pme_t*        GPU_FUNC_ARGUMENT(pme),
                                             DeviceBuffer<gmx::RVec> GPU_FUNC_ARGUMENT(d_x)) GPU_FUNC_TERM;

/*! \brief Get pointer to device copy of force data.
 * \param[in] pme            The PME data structure.
 * \returns                  Pointer to force data
 */
GPU_FUNC_QUALIFIER DeviceBuffer<gmx::RVec> pme_gpu_get_device_f(const gmx_pme_t* GPU_FUNC_ARGUMENT(pme))
        GPU_FUNC_TERM_WITH_RETURN(DeviceBuffer<gmx::RVec>{});

/*! \brief Get pointer to the device synchronizer object that allows syncing on PME force calculation completion
 * \param[in] pme            The PME data structure.
 * \returns                  Pointer to synchronizer
 */
GPU_FUNC_QUALIFIER GpuEventSynchronizer* pme_gpu_get_f_ready_synchronizer(const gmx_pme_t* GPU_FUNC_ARGUMENT(pme))
        GPU_FUNC_TERM_WITH_RETURN(nullptr);

#endif
