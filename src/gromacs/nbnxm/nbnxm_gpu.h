/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2012- The GROMACS Authors
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
 *  \brief Declare interface for GPU execution for NBNXN module
 *
 *  \author Szilard Pall <pall.szilard@gmail.com>
 *  \author Mark Abraham <mark.j.abraham@gmail.com>
 *  \ingroup module_nbnxm
 */

#ifndef GMX_NBNXM_NBNXM_GPU_H
#define GMX_NBNXM_NBNXM_GPU_H

#include "gromacs/gpu_utils/gpu_macros.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/locality.h"
#include "gromacs/nbnxm/atomdata.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

#include "nbnxm.h"

struct interaction_const_t;
struct nbnxn_atomdata_t;
struct gmx_wallcycle;
enum class GpuTaskCompletion;

namespace gmx
{
class ListedForcesGpu;
class StepWorkload;
} // namespace gmx

namespace Nbnxm
{

class Grid;

/*! \brief Returns true if LJ combination rules are used in the non-bonded kernels.
 *
 *  \param[in] vdwType  The VdW interaction/implementation type as defined by VdwType
 *                      enumeration.
 *
 * \returns Whether combination rules are used by the run.
 */
static inline bool useLjCombRule(const enum VdwType vdwType)
{
    return (vdwType == VdwType::CutCombGeom || vdwType == VdwType::CutCombLB);
}

/*! \brief
 * Launch asynchronously the xq buffer host to device copy.
 *
 * The nonlocal copy is skipped if there is no dependent work to do,
 * neither non-local nonbonded interactions nor bonded GPU work.
 *
 * \param [in]    nb        GPU nonbonded data.
 * \param [in]    nbdata    Host-side atom data structure.
 * \param [in]    aloc      Atom locality flag.
 */
GPU_FUNC_QUALIFIER
void gpu_copy_xq_to_gpu(NbnxmGpu gmx_unused*          nb,
                        const struct nbnxn_atomdata_t gmx_unused* nbdata,
                        gmx::AtomLocality gmx_unused              aloc) GPU_FUNC_TERM;

/*! \brief
 * Launch asynchronously the nonbonded force calculations.
 *
 *  Also launches the initial pruning of a fresh list after search.
 *
 *  The local and non-local interaction calculations are launched in two
 *  separate streams. If there is no work (i.e. empty pair list), the
 *  force kernel launch is omitted.
 *
 */
GPU_FUNC_QUALIFIER
void gpu_launch_kernel(NbnxmGpu gmx_unused*    nb,
                       const gmx::StepWorkload gmx_unused& stepWork,
                       gmx::InteractionLocality gmx_unused iloc) GPU_FUNC_TERM;

/*! \brief
 * Launch asynchronously the nonbonded prune-only kernel.
 *
 *  The local and non-local list pruning are launched in their separate streams.
 *
 *  Notes for future scheduling tuning:
 *  Currently we schedule the dynamic pruning between two MD steps *after* both local and
 *  nonlocal force D2H transfers completed. We could launch already after the cpyback
 *  is launched, but we want to avoid prune kernels (especially in the non-local
 *  high prio-stream) competing with nonbonded work.
 *
 *  However, this is not ideal as this schedule does not expose the available
 *  concurrency. The dynamic pruning kernel:
 *    - should be allowed to overlap with any task other than force compute, including
 *      transfers (F D2H and the next step's x H2D as well as force clearing).
 *    - we'd prefer to avoid competition with non-bonded force kernels belonging
 *      to the same rank and ideally other ranks too.
 *
 *  In the most general case, the former would require scheduling pruning in a separate
 *  stream and adding additional event sync points to ensure that force kernels read
 *  consistent pair list data. This would lead to some overhead (due to extra
 *  cudaStreamWaitEvent calls, 3-5 us/call) which we might be able to live with.
 *  The gains from additional overlap might not be significant as long as
 *  update+constraints anyway takes longer than pruning, but there will still
 *  be use-cases where more overlap may help (e.g. multiple ranks per GPU,
 *  no/hbonds only constraints).
 *  The above second point is harder to address given that multiple ranks will often
 *  share a GPU. Ranks that complete their nonbondeds sooner can schedule pruning earlier
 *  and without a third priority level it is difficult to avoid some interference of
 *  prune kernels with force tasks (in particular preemption of low-prio local force task).
 *
 * \param [inout] nb        GPU nonbonded data.
 * \param [in]    iloc      Interaction locality flag.
 * \param [in]    numParts  Number of parts the pair list is split into in the rolling kernel.
 */
GPU_FUNC_QUALIFIER
void gpu_launch_kernel_pruneonly(NbnxmGpu gmx_unused*                nb,
                                 gmx::InteractionLocality gmx_unused iloc,
                                 int gmx_unused                      numParts) GPU_FUNC_TERM;

/*! \brief
 * Launch asynchronously the download of short-range forces from the GPU
 * (and energies/shift forces if required).
 */
GPU_FUNC_QUALIFIER
void gpu_launch_cpyback(NbnxmGpu gmx_unused* nb,
                        nbnxn_atomdata_t gmx_unused* nbatom,
                        const gmx::StepWorkload gmx_unused& stepWork,
                        gmx::AtomLocality gmx_unused        aloc) GPU_FUNC_TERM;

/*! \brief Attempts to complete nonbonded GPU task.
 *
 *  This function attempts to complete the nonbonded task (both GPU and CPU auxiliary work).
 *  Success, i.e. that the tasks completed and results are ready to be consumed, is signaled
 *  by the return value (always true if blocking wait mode requested).
 *
 *  The \p completionKind parameter controls whether the behavior is non-blocking
 *  (achieved by passing GpuTaskCompletion::Check) or blocking wait until the results
 *  are ready (when GpuTaskCompletion::Wait is passed).
 *  As the "Check" mode the function will return immediately if the GPU stream
 *  still contain tasks that have not completed, it allows more flexible overlapping
 *  of work on the CPU with GPU execution.
 *
 *  Note that it is only safe to use the results, and to continue to the next MD
 *  step when this function has returned true which indicates successful completion of
 *  - All nonbonded GPU tasks: both compute and device transfer(s)
 *  - auxiliary tasks: updating the internal module state (timing accumulation, list pruning states) and
 *  - internal staging reduction of (\p fshift, \p e_el, \p e_lj).
 *
 * In GpuTaskCompletion::Check mode this function does the timing and keeps correct count
 * for the nonbonded task (incrementing only once per task), in the GpuTaskCompletion::Wait mode
 * timing is expected to be done in the caller.
 *
 *  TODO: improve the handling of outputs e.g. by ensuring that this function explcitly returns the
 *  force buffer (instead of that being passed only to nbnxn_gpu_launch_cpyback()) and by returning
 *  the energy and Fshift contributions for some external/centralized reduction.
 *
 * \param[in]  nb             The nonbonded data GPU structure
 * \param[in]  stepWork       Step schedule flags
 * \param[in]  aloc           Atom locality identifier
 * \param[out] e_lj           Pointer to the LJ energy output to accumulate into
 * \param[out] e_el           Pointer to the electrostatics energy output to accumulate into
 * \param[out] shiftForces    Shift forces buffer to accumulate into
 * \param[in]  completionKind Indicates whether nnbonded task completion should only be checked rather than waited for
 * \returns                   True if the nonbonded tasks associated with \p aloc locality have completed
 */
GPU_FUNC_QUALIFIER
bool gpu_try_finish_task(NbnxmGpu gmx_unused*    nb,
                         const gmx::StepWorkload gmx_unused& stepWork,
                         gmx::AtomLocality gmx_unused        aloc,
                         real gmx_unused* e_lj,
                         real gmx_unused*                    e_el,
                         gmx::ArrayRef<gmx::RVec> gmx_unused shiftForces,
                         GpuTaskCompletion gmx_unused completionKind) GPU_FUNC_TERM_WITH_RETURN(false);

/*! \brief  Completes the nonbonded GPU task blocking until GPU tasks and data
 * transfers to finish.
 *
 * Also does timing accounting and reduction of the internal staging buffers.
 * As this is called at the end of the step, it also resets the pair list and
 * pruning flags.
 *
 * \param[in] nb The nonbonded data GPU structure
 * \param[in]  stepWork        Step schedule flags
 * \param[in] aloc Atom locality identifier
 * \param[out] e_lj Pointer to the LJ energy output to accumulate into
 * \param[out] e_el Pointer to the electrostatics energy output to accumulate into
 * \param[out] shiftForces Shift forces buffer to accumulate into
 * \param[out] wcycle         Pointer to wallcycle data structure               */
GPU_FUNC_QUALIFIER
float gpu_wait_finish_task(NbnxmGpu gmx_unused*    nb,
                           const gmx::StepWorkload gmx_unused& stepWork,
                           gmx::AtomLocality gmx_unused        aloc,
                           real gmx_unused* e_lj,
                           real gmx_unused*                    e_el,
                           gmx::ArrayRef<gmx::RVec> gmx_unused shiftForces,
                           gmx_wallcycle gmx_unused* wcycle) GPU_FUNC_TERM_WITH_RETURN(0.0);

/*! \brief Initialization for X buffer operations on GPU.
 * Called on the NS step and performs (re-)allocations and memory copies. !*/
GPU_FUNC_QUALIFIER
void nbnxn_gpu_init_x_to_nbat_x(const Nbnxm::GridSet gmx_unused& gridSet,
                                NbnxmGpu gmx_unused* gpu_nbv) GPU_FUNC_TERM;

/*! \brief X buffer operations on GPU: performs conversion from rvec to nb format.
 *
 * \param[in]     grid             Grid to be converted.
 * \param[in,out] gpu_nbv          The nonbonded data GPU structure.
 * \param[in]     d_x              Device-side coordinates in plain rvec format.
 * \param[in]     xReadyOnDevice   Event synchronizer indicating that the coordinates are ready in
 * the device memory.
 * \param[in]     locality         Copy coordinates for local or non-local atoms.
 * \param[in]     gridId           Index of the grid being converted.
 * \param[in]     numColumnsMax    Maximum number of columns in the grid.
 * \param[in]     mustInsertNonLocalDependency Whether synchronization between local and non-local
 * streams should be added. Typically, true if and only if that is the last grid in gridset.
 */
GPU_FUNC_QUALIFIER
void nbnxn_gpu_x_to_nbat_x(const Nbnxm::Grid gmx_unused& grid,
                           NbnxmGpu gmx_unused*               gpu_nbv,
                           DeviceBuffer<gmx::RVec> gmx_unused d_x,
                           GpuEventSynchronizer gmx_unused* xReadyOnDevice,
                           gmx::AtomLocality gmx_unused     locality,
                           int gmx_unused                   gridId,
                           int gmx_unused                   numColumnsMax,
                           bool gmx_unused mustInsertNonLocalDependency) GPU_FUNC_TERM;

/*! \brief Sync the nonlocal stream with dependent tasks in the local queue.
 *
 *  As the point where the local stream tasks can be considered complete happens
 *  at the same call point where the nonlocal stream should be synced with the
 *  the local, this function records the event if called with the local stream as
 *  argument and inserts in the GPU stream a wait on the event on the nonlocal.
 *
 * \param[in] nb                   The nonbonded data GPU structure
 * \param[in] interactionLocality  Local or NonLocal sync point
 */
GPU_FUNC_QUALIFIER
void nbnxnInsertNonlocalGpuDependency(NbnxmGpu gmx_unused*                nb,
                                      gmx::InteractionLocality gmx_unused interactionLocality) GPU_FUNC_TERM;

/*! \brief Set up internal flags that indicate what type of short-range work there is.
 *
 * As nonbondeds and bondeds share input/output buffers and GPU queues,
 * both are considered when checking for work in the current domain.
 *
 * This function is expected to be called every time the work-distribution
 * can change (i.e. at search/domain decomposition steps).
 *
 * \param[inout]  nb               Pointer to the nonbonded GPU data structure
 * \param[in]     listedForcesGpu  Pointer to the GPU bonded data structure
 * \param[in]     iLocality        Interaction locality identifier
 */
GPU_FUNC_QUALIFIER
void setupGpuShortRangeWork(NbnxmGpu gmx_unused*       nb,
                            const gmx::ListedForcesGpu gmx_unused* listedForcesGpu,
                            gmx::InteractionLocality gmx_unused    iLocality) GPU_FUNC_TERM;

/*! \brief Returns true if there is GPU short-range work for the given interaction locality.
 *
 * Note that as, unlike nonbonded tasks, bonded tasks are not split into local/nonlocal,
 * and therefore if there are GPU offloaded bonded interactions, this function will return
 * true for both local and nonlocal atom range.
 *
 * \param[inout]  nb                   Pointer to the nonbonded GPU data structure
 * \param[in]     interactionLocality  Interaction locality identifier
 *
 * \return Whether there is short range work for a given locality.
 */
GPU_FUNC_QUALIFIER
bool haveGpuShortRangeWork(const NbnxmGpu gmx_unused* nb, gmx::InteractionLocality gmx_unused interactionLocality)
        GPU_FUNC_TERM_WITH_RETURN(false);

} // namespace Nbnxm
#endif
