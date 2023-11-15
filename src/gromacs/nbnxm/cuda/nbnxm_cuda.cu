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
/*! \file
 *  \brief Define CUDA implementation of nbnxn_gpu.h
 *
 *  \author Szilard Pall <pall.szilard@gmail.com>
 */
#include "gmxpre.h"

#include "config.h"

#include <cassert>
#include <cstdlib>

#include <cub/device/device_scan.cuh>

#include "gromacs/nbnxm/gpu_types_common.h"
#include "gromacs/nbnxm/nbnxm_gpu.h"

#if defined(_MSVC)
#    include <limits>
#endif


#include "gromacs/gpu_utils/devicebuffer.h"
#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/gpu_utils/gpueventsynchronizer.h"
#include "gromacs/gpu_utils/typecasts.cuh"
#include "gromacs/gpu_utils/vectype_ops.cuh"
#include "gromacs/hardware/device_information.h"
#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/nbnxm/atomdata.h"
#include "gromacs/nbnxm/gpu_common.h"
#include "gromacs/nbnxm/gpu_common_utils.h"
#include "gromacs/nbnxm/gpu_data_mgmt.h"
#include "gromacs/nbnxm/grid.h"
#include "gromacs/nbnxm/nbnxm.h"
#include "gromacs/nbnxm/pairlist.h"
#include "gromacs/timing/gpu_timing.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/gmxassert.h"

#include "nbnxm_cuda.h"
#include "nbnxm_cuda_types.h"

/***** The kernel declarations/definitions come here *****/


/* Top-level kernel declaration generation: will generate through multiple
 * inclusion the following flavors for all kernel declarations:
 * - force-only output;
 * - force and energy output;
 * - force-only with pair list pruning;
 * - force and energy output with pair list pruning.
 */
#define FUNCTION_DECLARATION_ONLY
/** Force only **/
#include "nbnxm_cuda_kernels.cuh"
/** Force & energy **/
#define CALC_ENERGIES
#include "nbnxm_cuda_kernels.cuh"
#undef CALC_ENERGIES

/*** Pair-list pruning kernels ***/
/** Force only **/
#define PRUNE_NBL
#include "nbnxm_cuda_kernels.cuh"
/** Force & energy **/
#define CALC_ENERGIES
#include "nbnxm_cuda_kernels.cuh"
#undef CALC_ENERGIES
#undef PRUNE_NBL

/* Prune-only kernels */
#include "nbnxm_cuda_kernel_pruneonly.cuh"
#undef FUNCTION_DECLARATION_ONLY

/* Now generate the function definitions if we are using a single compilation unit. */
#if GMX_CUDA_NB_SINGLE_COMPILATION_UNIT
#    include "nbnxm_cuda_kernel_F_noprune.cu"
#    include "nbnxm_cuda_kernel_F_prune.cu"
#    include "nbnxm_cuda_kernel_VF_noprune.cu"
#    include "nbnxm_cuda_kernel_VF_prune.cu"
#    include "nbnxm_cuda_kernel_pruneonly.cu"
#endif /* GMX_CUDA_NB_SINGLE_COMPILATION_UNIT */

#include "nbnxm_cuda_kernel_sci_sort.cuh"

namespace Nbnxm
{

/*! Nonbonded kernel function pointer type */
typedef void (*nbnxn_cu_kfunc_ptr_t)(const NBAtomDataGpu, const NBParamGpu, const gpu_plist, bool);

/*********************************/

/*! Returns the number of blocks to be used for the nonbonded GPU kernel. */
static inline int calc_nb_kernel_nblock(int nwork_units, const DeviceInformation* deviceInfo)
{
    int max_grid_x_size;

    assert(deviceInfo);
    /* CUDA does not accept grid dimension of 0 (which can happen e.g. with an
       empty domain) and that case should be handled before this point. */
    assert(nwork_units > 0);

    max_grid_x_size = deviceInfo->prop.maxGridSize[0];

    /* do we exceed the grid x dimension limit? */
    if (nwork_units > max_grid_x_size)
    {
        gmx_fatal(FARGS,
                  "Watch out, the input system is too large to simulate!\n"
                  "The number of nonbonded work units (=number of super-clusters) exceeds the"
                  "maximum grid size in x dimension (%d > %d)!",
                  nwork_units,
                  max_grid_x_size);
    }

    return nwork_units;
}


/* Constant arrays listing all kernel function pointers and enabling selection
   of a kernel in an elegant manner. */

/*! Pointers to the non-bonded kernels organized in 2-dim arrays by:
 *  electrostatics and VDW type.
 *
 *  Note that the row- and column-order of function pointers has to match the
 *  order of corresponding enumerated electrostatics and vdw types, resp.,
 *  defined in nbnxn_cuda_types.h.
 */

/*! Force-only kernel function pointers. */
static const nbnxn_cu_kfunc_ptr_t nb_kfunc_noener_noprune_ptr[c_numElecTypes][c_numVdwTypes] = {
    { nbnxn_kernel_ElecCut_VdwLJ_F_cuda,
      nbnxn_kernel_ElecCut_VdwLJCombGeom_F_cuda,
      nbnxn_kernel_ElecCut_VdwLJCombLB_F_cuda,
      nbnxn_kernel_ElecCut_VdwLJFsw_F_cuda,
      nbnxn_kernel_ElecCut_VdwLJPsw_F_cuda,
      nbnxn_kernel_ElecCut_VdwLJEwCombGeom_F_cuda,
      nbnxn_kernel_ElecCut_VdwLJEwCombLB_F_cuda },
    { nbnxn_kernel_ElecRF_VdwLJ_F_cuda,
      nbnxn_kernel_ElecRF_VdwLJCombGeom_F_cuda,
      nbnxn_kernel_ElecRF_VdwLJCombLB_F_cuda,
      nbnxn_kernel_ElecRF_VdwLJFsw_F_cuda,
      nbnxn_kernel_ElecRF_VdwLJPsw_F_cuda,
      nbnxn_kernel_ElecRF_VdwLJEwCombGeom_F_cuda,
      nbnxn_kernel_ElecRF_VdwLJEwCombLB_F_cuda },
    { nbnxn_kernel_ElecEwQSTab_VdwLJ_F_cuda,
      nbnxn_kernel_ElecEwQSTab_VdwLJCombGeom_F_cuda,
      nbnxn_kernel_ElecEwQSTab_VdwLJCombLB_F_cuda,
      nbnxn_kernel_ElecEwQSTab_VdwLJFsw_F_cuda,
      nbnxn_kernel_ElecEwQSTab_VdwLJPsw_F_cuda,
      nbnxn_kernel_ElecEwQSTab_VdwLJEwCombGeom_F_cuda,
      nbnxn_kernel_ElecEwQSTab_VdwLJEwCombLB_F_cuda },
    { nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJ_F_cuda,
      nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJCombGeom_F_cuda,
      nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJCombLB_F_cuda,
      nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJFsw_F_cuda,
      nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJPsw_F_cuda,
      nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJEwCombGeom_F_cuda,
      nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJEwCombLB_F_cuda },
    { nbnxn_kernel_ElecEw_VdwLJ_F_cuda,
      nbnxn_kernel_ElecEw_VdwLJCombGeom_F_cuda,
      nbnxn_kernel_ElecEw_VdwLJCombLB_F_cuda,
      nbnxn_kernel_ElecEw_VdwLJFsw_F_cuda,
      nbnxn_kernel_ElecEw_VdwLJPsw_F_cuda,
      nbnxn_kernel_ElecEw_VdwLJEwCombGeom_F_cuda,
      nbnxn_kernel_ElecEw_VdwLJEwCombLB_F_cuda },
    { nbnxn_kernel_ElecEwTwinCut_VdwLJ_F_cuda,
      nbnxn_kernel_ElecEwTwinCut_VdwLJCombGeom_F_cuda,
      nbnxn_kernel_ElecEwTwinCut_VdwLJCombLB_F_cuda,
      nbnxn_kernel_ElecEwTwinCut_VdwLJFsw_F_cuda,
      nbnxn_kernel_ElecEwTwinCut_VdwLJPsw_F_cuda,
      nbnxn_kernel_ElecEwTwinCut_VdwLJEwCombGeom_F_cuda,
      nbnxn_kernel_ElecEwTwinCut_VdwLJEwCombLB_F_cuda }
};

/*! Force + energy kernel function pointers. */
static const nbnxn_cu_kfunc_ptr_t nb_kfunc_ener_noprune_ptr[c_numElecTypes][c_numVdwTypes] = {
    { nbnxn_kernel_ElecCut_VdwLJ_VF_cuda,
      nbnxn_kernel_ElecCut_VdwLJCombGeom_VF_cuda,
      nbnxn_kernel_ElecCut_VdwLJCombLB_VF_cuda,
      nbnxn_kernel_ElecCut_VdwLJFsw_VF_cuda,
      nbnxn_kernel_ElecCut_VdwLJPsw_VF_cuda,
      nbnxn_kernel_ElecCut_VdwLJEwCombGeom_VF_cuda,
      nbnxn_kernel_ElecCut_VdwLJEwCombLB_VF_cuda },
    { nbnxn_kernel_ElecRF_VdwLJ_VF_cuda,
      nbnxn_kernel_ElecRF_VdwLJCombGeom_VF_cuda,
      nbnxn_kernel_ElecRF_VdwLJCombLB_VF_cuda,
      nbnxn_kernel_ElecRF_VdwLJFsw_VF_cuda,
      nbnxn_kernel_ElecRF_VdwLJPsw_VF_cuda,
      nbnxn_kernel_ElecRF_VdwLJEwCombGeom_VF_cuda,
      nbnxn_kernel_ElecRF_VdwLJEwCombLB_VF_cuda },
    { nbnxn_kernel_ElecEwQSTab_VdwLJ_VF_cuda,
      nbnxn_kernel_ElecEwQSTab_VdwLJCombGeom_VF_cuda,
      nbnxn_kernel_ElecEwQSTab_VdwLJCombLB_VF_cuda,
      nbnxn_kernel_ElecEwQSTab_VdwLJFsw_VF_cuda,
      nbnxn_kernel_ElecEwQSTab_VdwLJPsw_VF_cuda,
      nbnxn_kernel_ElecEwQSTab_VdwLJEwCombGeom_VF_cuda,
      nbnxn_kernel_ElecEwQSTab_VdwLJEwCombLB_VF_cuda },
    { nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJ_VF_cuda,
      nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJCombGeom_VF_cuda,
      nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJCombLB_VF_cuda,
      nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJFsw_VF_cuda,
      nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJPsw_VF_cuda,
      nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJEwCombGeom_VF_cuda,
      nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJEwCombLB_VF_cuda },
    { nbnxn_kernel_ElecEw_VdwLJ_VF_cuda,
      nbnxn_kernel_ElecEw_VdwLJCombGeom_VF_cuda,
      nbnxn_kernel_ElecEw_VdwLJCombLB_VF_cuda,
      nbnxn_kernel_ElecEw_VdwLJFsw_VF_cuda,
      nbnxn_kernel_ElecEw_VdwLJPsw_VF_cuda,
      nbnxn_kernel_ElecEw_VdwLJEwCombGeom_VF_cuda,
      nbnxn_kernel_ElecEw_VdwLJEwCombLB_VF_cuda },
    { nbnxn_kernel_ElecEwTwinCut_VdwLJ_VF_cuda,
      nbnxn_kernel_ElecEwTwinCut_VdwLJCombGeom_VF_cuda,
      nbnxn_kernel_ElecEwTwinCut_VdwLJCombLB_VF_cuda,
      nbnxn_kernel_ElecEwTwinCut_VdwLJFsw_VF_cuda,
      nbnxn_kernel_ElecEwTwinCut_VdwLJPsw_VF_cuda,
      nbnxn_kernel_ElecEwTwinCut_VdwLJEwCombGeom_VF_cuda,
      nbnxn_kernel_ElecEwTwinCut_VdwLJEwCombLB_VF_cuda }
};

/*! Force + pruning kernel function pointers. */
static const nbnxn_cu_kfunc_ptr_t nb_kfunc_noener_prune_ptr[c_numElecTypes][c_numVdwTypes] = {
    { nbnxn_kernel_ElecCut_VdwLJ_F_prune_cuda,
      nbnxn_kernel_ElecCut_VdwLJCombGeom_F_prune_cuda,
      nbnxn_kernel_ElecCut_VdwLJCombLB_F_prune_cuda,
      nbnxn_kernel_ElecCut_VdwLJFsw_F_prune_cuda,
      nbnxn_kernel_ElecCut_VdwLJPsw_F_prune_cuda,
      nbnxn_kernel_ElecCut_VdwLJEwCombGeom_F_prune_cuda,
      nbnxn_kernel_ElecCut_VdwLJEwCombLB_F_prune_cuda },
    { nbnxn_kernel_ElecRF_VdwLJ_F_prune_cuda,
      nbnxn_kernel_ElecRF_VdwLJCombGeom_F_prune_cuda,
      nbnxn_kernel_ElecRF_VdwLJCombLB_F_prune_cuda,
      nbnxn_kernel_ElecRF_VdwLJFsw_F_prune_cuda,
      nbnxn_kernel_ElecRF_VdwLJPsw_F_prune_cuda,
      nbnxn_kernel_ElecRF_VdwLJEwCombGeom_F_prune_cuda,
      nbnxn_kernel_ElecRF_VdwLJEwCombLB_F_prune_cuda },
    { nbnxn_kernel_ElecEwQSTab_VdwLJ_F_prune_cuda,
      nbnxn_kernel_ElecEwQSTab_VdwLJCombGeom_F_prune_cuda,
      nbnxn_kernel_ElecEwQSTab_VdwLJCombLB_F_prune_cuda,
      nbnxn_kernel_ElecEwQSTab_VdwLJFsw_F_prune_cuda,
      nbnxn_kernel_ElecEwQSTab_VdwLJPsw_F_prune_cuda,
      nbnxn_kernel_ElecEwQSTab_VdwLJEwCombGeom_F_prune_cuda,
      nbnxn_kernel_ElecEwQSTab_VdwLJEwCombLB_F_prune_cuda },
    { nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJ_F_prune_cuda,
      nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJCombGeom_F_prune_cuda,
      nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJCombLB_F_prune_cuda,
      nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJFsw_F_prune_cuda,
      nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJPsw_F_prune_cuda,
      nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJEwCombGeom_F_prune_cuda,
      nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJEwCombLB_F_prune_cuda },
    { nbnxn_kernel_ElecEw_VdwLJ_F_prune_cuda,
      nbnxn_kernel_ElecEw_VdwLJCombGeom_F_prune_cuda,
      nbnxn_kernel_ElecEw_VdwLJCombLB_F_prune_cuda,
      nbnxn_kernel_ElecEw_VdwLJFsw_F_prune_cuda,
      nbnxn_kernel_ElecEw_VdwLJPsw_F_prune_cuda,
      nbnxn_kernel_ElecEw_VdwLJEwCombGeom_F_prune_cuda,
      nbnxn_kernel_ElecEw_VdwLJEwCombLB_F_prune_cuda },
    { nbnxn_kernel_ElecEwTwinCut_VdwLJ_F_prune_cuda,
      nbnxn_kernel_ElecEwTwinCut_VdwLJCombGeom_F_prune_cuda,
      nbnxn_kernel_ElecEwTwinCut_VdwLJCombLB_F_prune_cuda,
      nbnxn_kernel_ElecEwTwinCut_VdwLJFsw_F_prune_cuda,
      nbnxn_kernel_ElecEwTwinCut_VdwLJPsw_F_prune_cuda,
      nbnxn_kernel_ElecEwTwinCut_VdwLJEwCombGeom_F_prune_cuda,
      nbnxn_kernel_ElecEwTwinCut_VdwLJEwCombLB_F_prune_cuda }
};

/*! Force + energy + pruning kernel function pointers. */
static const nbnxn_cu_kfunc_ptr_t nb_kfunc_ener_prune_ptr[c_numElecTypes][c_numVdwTypes] = {
    { nbnxn_kernel_ElecCut_VdwLJ_VF_prune_cuda,
      nbnxn_kernel_ElecCut_VdwLJCombGeom_VF_prune_cuda,
      nbnxn_kernel_ElecCut_VdwLJCombLB_VF_prune_cuda,
      nbnxn_kernel_ElecCut_VdwLJFsw_VF_prune_cuda,
      nbnxn_kernel_ElecCut_VdwLJPsw_VF_prune_cuda,
      nbnxn_kernel_ElecCut_VdwLJEwCombGeom_VF_prune_cuda,
      nbnxn_kernel_ElecCut_VdwLJEwCombLB_VF_prune_cuda },
    { nbnxn_kernel_ElecRF_VdwLJ_VF_prune_cuda,
      nbnxn_kernel_ElecRF_VdwLJCombGeom_VF_prune_cuda,
      nbnxn_kernel_ElecRF_VdwLJCombLB_VF_prune_cuda,
      nbnxn_kernel_ElecRF_VdwLJFsw_VF_prune_cuda,
      nbnxn_kernel_ElecRF_VdwLJPsw_VF_prune_cuda,
      nbnxn_kernel_ElecRF_VdwLJEwCombGeom_VF_prune_cuda,
      nbnxn_kernel_ElecRF_VdwLJEwCombLB_VF_prune_cuda },
    { nbnxn_kernel_ElecEwQSTab_VdwLJ_VF_prune_cuda,
      nbnxn_kernel_ElecEwQSTab_VdwLJCombGeom_VF_prune_cuda,
      nbnxn_kernel_ElecEwQSTab_VdwLJCombLB_VF_prune_cuda,
      nbnxn_kernel_ElecEwQSTab_VdwLJFsw_VF_prune_cuda,
      nbnxn_kernel_ElecEwQSTab_VdwLJPsw_VF_prune_cuda,
      nbnxn_kernel_ElecEwQSTab_VdwLJEwCombGeom_VF_prune_cuda,
      nbnxn_kernel_ElecEwQSTab_VdwLJEwCombLB_VF_prune_cuda },
    { nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJ_VF_prune_cuda,
      nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJCombGeom_VF_prune_cuda,
      nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJCombLB_VF_prune_cuda,
      nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJFsw_VF_prune_cuda,
      nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJPsw_VF_prune_cuda,
      nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJEwCombGeom_VF_prune_cuda,
      nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJEwCombLB_VF_prune_cuda },
    { nbnxn_kernel_ElecEw_VdwLJ_VF_prune_cuda,
      nbnxn_kernel_ElecEw_VdwLJCombGeom_VF_prune_cuda,
      nbnxn_kernel_ElecEw_VdwLJCombLB_VF_prune_cuda,
      nbnxn_kernel_ElecEw_VdwLJFsw_VF_prune_cuda,
      nbnxn_kernel_ElecEw_VdwLJPsw_VF_prune_cuda,
      nbnxn_kernel_ElecEw_VdwLJEwCombGeom_VF_prune_cuda,
      nbnxn_kernel_ElecEw_VdwLJEwCombLB_VF_prune_cuda },
    { nbnxn_kernel_ElecEwTwinCut_VdwLJ_VF_prune_cuda,
      nbnxn_kernel_ElecEwTwinCut_VdwLJCombGeom_VF_prune_cuda,
      nbnxn_kernel_ElecEwTwinCut_VdwLJCombLB_VF_prune_cuda,
      nbnxn_kernel_ElecEwTwinCut_VdwLJFsw_VF_prune_cuda,
      nbnxn_kernel_ElecEwTwinCut_VdwLJPsw_VF_prune_cuda,
      nbnxn_kernel_ElecEwTwinCut_VdwLJEwCombGeom_VF_prune_cuda,
      nbnxn_kernel_ElecEwTwinCut_VdwLJEwCombLB_VF_prune_cuda }
};

/*! Return a pointer to the kernel version to be executed at the current step. */
static inline nbnxn_cu_kfunc_ptr_t select_nbnxn_kernel(enum ElecType           elecType,
                                                       enum VdwType            vdwType,
                                                       bool                    bDoEne,
                                                       bool                    bDoPrune,
                                                       const DeviceInformation gmx_unused* deviceInfo)
{
    const int elecTypeIdx = static_cast<int>(elecType);
    const int vdwTypeIdx  = static_cast<int>(vdwType);

    GMX_ASSERT(elecTypeIdx < c_numElecTypes,
               "The electrostatics type requested is not implemented in the CUDA kernels.");
    GMX_ASSERT(vdwTypeIdx < c_numVdwTypes,
               "The VdW type requested is not implemented in the CUDA kernels.");

    /* assert assumptions made by the kernels */
    GMX_ASSERT(c_nbnxnGpuClusterSize * c_nbnxnGpuClusterSize / c_nbnxnGpuClusterpairSplit
                       == deviceInfo->prop.warpSize,
               "The CUDA kernels require the "
               "cluster_size_i*cluster_size_j/nbnxn_gpu_clusterpair_split to match the warp size "
               "of the architecture targeted.");

    if (bDoEne)
    {
        if (bDoPrune)
        {
            return nb_kfunc_ener_prune_ptr[elecTypeIdx][vdwTypeIdx];
        }
        else
        {
            return nb_kfunc_ener_noprune_ptr[elecTypeIdx][vdwTypeIdx];
        }
    }
    else
    {
        if (bDoPrune)
        {
            return nb_kfunc_noener_prune_ptr[elecTypeIdx][vdwTypeIdx];
        }
        else
        {
            return nb_kfunc_noener_noprune_ptr[elecTypeIdx][vdwTypeIdx];
        }
    }
}

/*! \brief Calculates the amount of shared memory required by the nonbonded kernel in use. */
static inline int calc_shmem_required_nonbonded(const int               num_threads_z,
                                                const DeviceInformation gmx_unused* deviceInfo,
                                                const NBParamGpu*                   nbp)
{
    int shmem;

    assert(deviceInfo);

    /* size of shmem (force-buffers/xq/atom type preloading) */
    /* NOTE: with the default kernel on sm3.0 we need shmem only for pre-loading */
    /* i-atom x+q in shared memory */
    shmem = c_nbnxnGpuNumClusterPerSupercluster * c_clSize * sizeof(float4);
    /* cj in shared memory, for each warp separately */
    shmem += num_threads_z * c_nbnxnGpuClusterpairSplit * c_nbnxnGpuJgroupSize * sizeof(int);

    if (nbp->vdwType == VdwType::CutCombGeom || nbp->vdwType == VdwType::CutCombLB)
    {
        /* i-atom LJ combination parameters in shared memory */
        shmem += c_nbnxnGpuNumClusterPerSupercluster * c_clSize * sizeof(float2);
    }
    else
    {
        /* i-atom types in shared memory */
        shmem += c_nbnxnGpuNumClusterPerSupercluster * c_clSize * sizeof(int);
    }
    /* for reducing prunedPairListCount over all warps in the block, to be used in plist sorting */
    shmem += 1 * sizeof(int);

    return shmem;
}


/*! \brief Calculates the amount of shared memory required by the nonbonded kernel in use.
 *
 * Take counts prepared in combined prune and interaction kernel and use them to sort plist.
 * Note that this sorted list is not available in the combined prune and interaction kernel
 * itself, which causes a performance degredation of 1-10% for that initial call */
static inline void gpuLaunchKernelSciSort(gpu_plist* plist, const DeviceStream& deviceStream)
{
    size_t scanTemporarySize = static_cast<size_t>(plist->sorting.nscanTemporary);

    cub::DeviceScan::ExclusiveSum(plist->sorting.scanTemporary,
                                  scanTemporarySize,
                                  plist->sorting.sciHistogram,
                                  plist->sorting.sciOffset,
                                  c_sciHistogramSize,
                                  deviceStream.stream());

    KernelLaunchConfig configSortSci;
    configSortSci.blockSize[0] = c_sciSortingThreadsPerBlock;
    configSortSci.blockSize[1] = 1;
    configSortSci.blockSize[2] = 1;
    configSortSci.gridSize[0] =
            (plist->nsci + c_sciSortingThreadsPerBlock - 1) / c_sciSortingThreadsPerBlock;
    configSortSci.sharedMemorySize = 0;

    const auto kernelSciSort = nbnxnKernelBucketSciSort;

    const auto kernelSciSortArgs = prepareGpuKernelArguments(kernelSciSort, configSortSci, plist);

    launchGpuKernel(kernelSciSort, configSortSci, deviceStream, nullptr, "nbnxn_kernel_sci_sort", kernelSciSortArgs);
}


/*! As we execute nonbonded workload in separate streams, before launching
   the kernel we need to make sure that he following operations have completed:
   - atomdata allocation and related H2D transfers (every nstlist step);
   - pair list H2D transfer (every nstlist step);
   - shift vector H2D transfer (every nstlist step);
   - force (+shift force and energy) output clearing (every step).

   These operations are issued in the local stream at the beginning of the step
   and therefore always complete before the local kernel launch. The non-local
   kernel is launched after the local on the same device/context hence it is
   inherently scheduled after the operations in the local stream (including the
   above "misc_ops") on pre-GK110 devices with single hardware queue, but on later
   devices with multiple hardware queues the dependency needs to be enforced.
   We use the misc_ops_and_local_H2D_done event to record the point where
   the local x+q H2D (and all preceding) tasks are complete and synchronize
   with this event in the non-local stream before launching the non-bonded kernel.
 */
void gpu_launch_kernel(NbnxmGpu* nb, const gmx::StepWorkload& stepWork, const InteractionLocality iloc)
{
    NBAtomDataGpu*      adat         = nb->atdat;
    NBParamGpu*         nbp          = nb->nbparam;
    gpu_plist*          plist        = nb->plist[iloc];
    Nbnxm::GpuTimers*   timers       = nb->timers;
    const DeviceStream& deviceStream = *nb->deviceStreams[iloc];

    bool bDoTime = nb->bDoTime;

    /* Don't launch the non-local kernel if there is no work to do.
       Doing the same for the local kernel is more complicated, since the
       local part of the force array also depends on the non-local kernel.
       So to avoid complicating the code and to reduce the risk of bugs,
       we always call the local kernel, and later (not in
       this function) the stream wait, local f copyback and the f buffer
       clearing. All these operations, except for the local interaction kernel,
       are needed for the non-local interactions. The skip of the local kernel
       call is taken care of later in this function. */
    if (canSkipNonbondedWork(*nb, iloc))
    {
        plist->haveFreshList = false;

        return;
    }

    if (nbp->useDynamicPruning && plist->haveFreshList)
    {
        /* Prunes for rlistOuter and rlistInner, sets plist->haveFreshList=false
           (TODO: ATM that's the way the timing accounting can distinguish between
           separate prune kernel and combined force+prune, maybe we need a better way?).
         */
        gpu_launch_kernel_pruneonly(nb, iloc, 1);
    }

    if (plist->nsci == 0)
    {
        /* Don't launch an empty local kernel (not allowed with CUDA) */
        return;
    }

    /* beginning of timed nonbonded calculation section */
    if (bDoTime)
    {
        timers->interaction[iloc].nb_k.openTimingRegion(deviceStream);
    }

    /* Kernel launch config:
     * - The thread block dimensions match the size of i-clusters, j-clusters,
     *   and j-cluster concurrency, in x, y, and z, respectively.
     * - The 1D block-grid contains as many blocks as super-clusters.
     */
    int num_threads_z = 1;
    if (nb->deviceContext_->deviceInfo().prop.major == 3 && nb->deviceContext_->deviceInfo().prop.minor == 7)
    {
        num_threads_z = 2;
    }
    int nblock = calc_nb_kernel_nblock(plist->nsci, &nb->deviceContext_->deviceInfo());


    KernelLaunchConfig config;
    config.blockSize[0] = c_clSize;
    config.blockSize[1] = c_clSize;
    config.blockSize[2] = num_threads_z;
    config.gridSize[0]  = nblock;
    config.sharedMemorySize =
            calc_shmem_required_nonbonded(num_threads_z, &nb->deviceContext_->deviceInfo(), nbp);

    if (debug)
    {
        fprintf(debug,
                "Non-bonded GPU launch configuration:\n\tThread block: %zux%zux%zu\n\t"
                "\tGrid: %zux%zu\n\t#Super-clusters/clusters: %d/%d (%d)\n"
                "\tShMem: %zu\n",
                config.blockSize[0],
                config.blockSize[1],
                config.blockSize[2],
                config.gridSize[0],
                config.gridSize[1],
                plist->nsci * c_nbnxnGpuNumClusterPerSupercluster,
                c_nbnxnGpuNumClusterPerSupercluster,
                plist->na_c,
                config.sharedMemorySize);
    }

    auto* timingEvent = bDoTime ? timers->interaction[iloc].nb_k.fetchNextEvent() : nullptr;

    /* Whether we need to call a combined prune and interaction kernel or just an interaction
     * kernel. bDoPrune being true implies we are not using dynamic pruning and are in the first
     * call to the interaction kernel after a neighbour list step */
    bool       bDoPrune = (plist->haveFreshList && !nb->timers->interaction[iloc].didPrune);
    const auto kernel   = select_nbnxn_kernel(
            nbp->elecType, nbp->vdwType, stepWork.computeEnergy, bDoPrune, &nb->deviceContext_->deviceInfo());
    const auto kernelArgs =
            prepareGpuKernelArguments(kernel, config, adat, nbp, plist, &stepWork.computeVirial);
    launchGpuKernel(kernel, config, deviceStream, timingEvent, "k_calc_nb", kernelArgs);

    if (bDoPrune)
    {
        gpuLaunchKernelSciSort(plist, deviceStream);
    }


    if (bDoTime)
    {
        timers->interaction[iloc].nb_k.closeTimingRegion(deviceStream);
    }

    if (GMX_NATIVE_WINDOWS)
    {
        /* Windows: force flushing WDDM queue */
        cudaStreamQuery(deviceStream.stream());
    }
}

/*! Calculates the amount of shared memory required by the CUDA kernel in use. */
static inline int calc_shmem_required_prune(const int num_threads_z)
{
    int shmem;

    /* i-atom x in shared memory */
    shmem = c_nbnxnGpuNumClusterPerSupercluster * c_clSize * sizeof(float4);
    /* cj in shared memory, for each warp separately */
    shmem += num_threads_z * c_nbnxnGpuClusterpairSplit * c_nbnxnGpuJgroupSize * sizeof(int);

    return shmem;
}

void gpu_launch_kernel_pruneonly(NbnxmGpu* nb, const InteractionLocality iloc, const int numParts)
{
    NBAtomDataGpu*      adat         = nb->atdat;
    NBParamGpu*         nbp          = nb->nbparam;
    gpu_plist*          plist        = nb->plist[iloc];
    Nbnxm::GpuTimers*   timers       = nb->timers;
    const DeviceStream& deviceStream = *nb->deviceStreams[iloc];

    bool bDoTime = nb->bDoTime;

    if (plist->haveFreshList)
    {
        GMX_ASSERT(numParts == 1, "With first pruning we expect 1 part");

        /* Set rollingPruningNumParts to signal that it is not set */
        plist->rollingPruningNumParts = 0;
    }
    else
    {
        if (plist->rollingPruningNumParts == 0)
        {
            plist->rollingPruningNumParts = numParts;
        }
        else
        {
            GMX_ASSERT(numParts == plist->rollingPruningNumParts,
                       "It is not allowed to change numParts in between list generation steps");
        }
    }

    /* Compute the max number of list entries to prune across all passes
     * Note that the actual number for a specific pass will be computed inside the kernel.
     * Also note that this CUDA implementation (parts tracking on device) differs from the
     * other backends (parts tracking on host, passed as kernel argument).
     */
    int numSciInPartMax = (plist->nsci) / numParts;

    /* Don't launch the kernel if there is no work to do (not allowed with CUDA) */
    if (numSciInPartMax <= 0)
    {
        plist->haveFreshList = false;

        return;
    }

    GpuRegionTimer* timer = nullptr;
    if (bDoTime)
    {
        timer = &(plist->haveFreshList ? timers->interaction[iloc].prune_k
                                       : timers->interaction[iloc].rollingPrune_k);
    }

    /* beginning of timed prune calculation section */
    if (bDoTime)
    {
        timer->openTimingRegion(deviceStream);
    }

    /* Kernel launch config:
     * - The thread block dimensions match the size of i-clusters, j-clusters,
     *   and j-cluster concurrency, in x, y, and z, respectively.
     * - The 1D block-grid contains as many blocks as super-clusters.
     */
    int num_threads_z = c_pruneKernelJPackedConcurrency;
    int nblock        = calc_nb_kernel_nblock(numSciInPartMax, &nb->deviceContext_->deviceInfo());

    KernelLaunchConfig config;
    config.blockSize[0]     = c_clSize;
    config.blockSize[1]     = c_clSize;
    config.blockSize[2]     = num_threads_z;
    config.gridSize[0]      = nblock;
    config.sharedMemorySize = calc_shmem_required_prune(num_threads_z);

    if (debug)
    {
        fprintf(debug,
                "Pruning GPU kernel launch configuration:\n\tThread block: %zux%zux%zu\n\t"
                "\tGrid: %zux%zu\n\t#Super-clusters/clusters: %d/%d (%d)\n"
                "\tShMem: %zu\n",
                config.blockSize[0],
                config.blockSize[1],
                config.blockSize[2],
                config.gridSize[0],
                config.gridSize[1],
                numSciInPartMax * c_nbnxnGpuNumClusterPerSupercluster,
                c_nbnxnGpuNumClusterPerSupercluster,
                plist->na_c,
                config.sharedMemorySize);
    }

    auto*          timingEvent  = bDoTime ? timer->fetchNextEvent() : nullptr;
    constexpr char kernelName[] = "k_pruneonly";
    const auto     kernel =
            plist->haveFreshList ? nbnxn_kernel_prune_cuda<true> : nbnxn_kernel_prune_cuda<false>;
    const auto kernelArgs = prepareGpuKernelArguments(kernel, config, adat, nbp, plist, &numParts);
    launchGpuKernel(kernel, config, deviceStream, timingEvent, kernelName, kernelArgs);

    if (plist->haveFreshList)
    {
        gpuLaunchKernelSciSort(plist, deviceStream);
    }

    /* TODO: consider a more elegant way to track which kernel has been called
       (combined or separate 1st pass prune, rolling prune). */
    if (plist->haveFreshList)
    {
        plist->haveFreshList = false;
        /* Mark that pruning has been done */
        nb->timers->interaction[iloc].didPrune = true;
    }
    else
    {
        /* Mark that rolling pruning has been done */
        nb->timers->interaction[iloc].didRollingPrune = true;
    }

    if (bDoTime)
    {
        timer->closeTimingRegion(deviceStream);
    }

    if (GMX_NATIVE_WINDOWS)
    {
        /* Windows: force flushing WDDM queue */
        cudaStreamQuery(deviceStream.stream());
    }
}

void cuda_set_cacheconfig()
{
    cudaError_t stat;

    for (int i = 0; i < c_numElecTypes; i++)
    {
        for (int j = 0; j < c_numVdwTypes; j++)
        {
            /* Default kernel 32/32 kB Shared/L1 */
            cudaFuncSetCacheConfig(nb_kfunc_ener_prune_ptr[i][j], cudaFuncCachePreferEqual);
            cudaFuncSetCacheConfig(nb_kfunc_ener_noprune_ptr[i][j], cudaFuncCachePreferEqual);
            cudaFuncSetCacheConfig(nb_kfunc_noener_prune_ptr[i][j], cudaFuncCachePreferEqual);
            stat = cudaFuncSetCacheConfig(nb_kfunc_noener_noprune_ptr[i][j], cudaFuncCachePreferEqual);
            CU_RET_ERR(stat, "cudaFuncSetCacheConfig failed");
        }
    }
}

} // namespace Nbnxm
