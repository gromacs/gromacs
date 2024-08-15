/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2014- The GROMACS Authors
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
 *  \brief Declare interface for GPU data transfer for NBNXN module
 *
 *  \author Szilard Pall <pall.szilard@gmail.com>
 *  \author Mark Abraham <mark.j.abraham@gmail.com>
 *  \ingroup module_nbnxm
 *  \inlibraryapi
 */

#ifndef GMX_NBNXN_GPU_DATA_MGMT_H
#define GMX_NBNXN_GPU_DATA_MGMT_H

#include <memory>

#include "gromacs/gpu_utils/gpu_macros.h"
#include "gromacs/mdtypes/locality.h"

#include "nbnxm.h"

struct gmx_wallclock_gpu_nbnxn_t;
struct interaction_const_t;
class DeviceStream;

namespace gmx
{
struct NbnxmGpu;
struct NBAtomDataGpu;
struct nbnxn_atomdata_t;
struct NbnxnPairlistGpu;
struct PairlistParams;
class DeviceStreamManager;

class GpuPairlist;

/** Initializes the data structures related to GPU nonbonded calculations. */
GPU_FUNC_QUALIFIER
NbnxmGpu* gpu_init(const DeviceStreamManager gmx_unused& deviceStreamManager,
                   const interaction_const_t gmx_unused* ic,
                   const PairlistParams gmx_unused& listParams,
                   const nbnxn_atomdata_t gmx_unused* nbat,
                   /* true if both local and non-local are done on GPU */
                   bool gmx_unused bLocalAndNonlocal) GPU_FUNC_TERM_WITH_RETURN(nullptr);

/** Initializes pair-list data for GPU, called at every pair search step. */
GPU_FUNC_QUALIFIER
void gpu_init_pairlist(NbnxmGpu gmx_unused*          nb,
                       const struct NbnxnPairlistGpu gmx_unused* h_nblist,
                       InteractionLocality gmx_unused            iloc) GPU_FUNC_TERM;

/** Initializes atom-data on the GPU, called at every pair search step. */
GPU_FUNC_QUALIFIER
void gpu_init_atomdata(NbnxmGpu gmx_unused* nb, const nbnxn_atomdata_t gmx_unused* nbat) GPU_FUNC_TERM;

/*! \brief Re-generate the GPU Ewald force table, resets rlist, and update the
 *  electrostatic type switching to twin cut-off (or back) if needed.
 */
GPU_FUNC_QUALIFIER
void gpu_pme_loadbal_update_param(struct nonbonded_verlet_t gmx_unused* nbv,
                                  const interaction_const_t gmx_unused& ic) GPU_FUNC_TERM;

/** Uploads shift vector to the GPU if the box is dynamic (otherwise just returns). */
GPU_FUNC_QUALIFIER
void gpu_upload_shiftvec(NbnxmGpu gmx_unused* nb, const nbnxn_atomdata_t gmx_unused* nbatom) GPU_FUNC_TERM;

/** Clears GPU outputs: nonbonded force, shift force and energy. */
GPU_FUNC_QUALIFIER
void gpu_clear_outputs(NbnxmGpu gmx_unused* nb, bool gmx_unused computeVirial) GPU_FUNC_TERM;

/** Frees all GPU resources used for the nonbonded calculations. */
GPU_FUNC_QUALIFIER
void gpu_free(NbnxmGpu gmx_unused* nb) GPU_FUNC_TERM;

/** Returns the GPU timings structure or NULL if GPU is not used or timing is off. */
GPU_FUNC_QUALIFIER
struct gmx_wallclock_gpu_nbnxn_t* gpu_get_timings(NbnxmGpu gmx_unused* nb)
        GPU_FUNC_TERM_WITH_RETURN(nullptr);

/** Resets nonbonded GPU timings. */
GPU_FUNC_QUALIFIER
void gpu_reset_timings(struct nonbonded_verlet_t gmx_unused* nbv) GPU_FUNC_TERM;

/** Calculates the minimum size of proximity lists to improve SM load balance
 *  with GPU non-bonded kernels. */
GPU_FUNC_QUALIFIER
int gpu_min_ci_balanced(NbnxmGpu gmx_unused* nb) GPU_FUNC_TERM_WITH_RETURN(-1);

/** Returns if analytical Ewald GPU kernels are used. */
GPU_FUNC_QUALIFIER
bool gpu_is_kernel_ewald_analytical(const NbnxmGpu gmx_unused* nb) GPU_FUNC_TERM_WITH_RETURN(FALSE);

/** Returns an opaque pointer to the GPU NBNXM atom data.
 */
GPU_FUNC_QUALIFIER
NBAtomDataGpu* gpuGetNBAtomData(NbnxmGpu gmx_unused* nb) GPU_FUNC_TERM_WITH_RETURN(nullptr);

/** Returns forces device buffer.
 */
GPU_FUNC_QUALIFIER
DeviceBuffer<RVec> gpu_get_f(NbnxmGpu gmx_unused* nb) GPU_FUNC_TERM_WITH_RETURN(DeviceBuffer<RVec>{});

/*! \brief Calculates working memory required for exclusive sum, used in neighbour list sorting on GPU.
 *
 * This is only used for CUDA/HIP, where the actual size is calculate based on the list.
 * For SYCL, the default value of 0 is important for the code to work correctly, this is why we have it set here.
 * */
CUDA_FUNC_QUALIFIER
size_t getExclusiveScanWorkingArraySize(GpuPairlist*        CUDA_FUNC_ARGUMENT(plist),
                                        const DeviceStream& CUDA_FUNC_ARGUMENT(deviceStream))
        CUDA_FUNC_TERM_WITH_RETURN(0);

/*! \brief Perform exclusive scan to obtain input for sci sorting. */
CUDA_FUNC_QUALIFIER
void performExclusiveScan(size_t              CUDA_FUNC_ARGUMENT(temporaryBufferSize),
                          char*               CUDA_FUNC_ARGUMENT(temporaryBuffer),
                          GpuPairlist*        CUDA_FUNC_ARGUMENT(plist),
                          const DeviceStream& CUDA_FUNC_ARGUMENT(deviceStream)) CUDA_FUNC_TERM;

} // namespace gmx

#endif
