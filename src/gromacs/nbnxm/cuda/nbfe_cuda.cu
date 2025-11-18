/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2025- The GROMACS Authors
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
 *  \brief Define CUDA implementation of nonbonded free energy calculations
 *
 *  \author Yiqi Chen <yiqi.chen@metax-tech.com>
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
#include "gromacs/gpu_utils/typecasts_cuda_hip.h"
#include "gromacs/gpu_utils/vectype_ops_cuda.h"
#include "gromacs/hardware/device_information.h"
#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/gmxassert.h"

#include "nbnxm_cuda.h"
#include "nbnxm_cuda_kernel_utils.cuh"
#include "nbnxm_cuda_types.h"

/***** The kernel declarations/definitions come here *****/


/* Top-level kernel declaration generation: will generate through multiple
 * inclusion the following flavors for all kernel declarations:
 * - force-only output;
 * - force and energy output;
 */
#define FUNCTION_DECLARATION_ONLY
/** Force only **/
#include "nbfe_cuda_kernels.cuh"
/** Force & energy **/
#define CALC_ENERGIES
#include "nbfe_cuda_kernels.cuh"
#undef CALC_ENERGIES

#include "nbfe_foreign_cuda_kernels.cuh"
#undef FUNCTION_DECLARATION_ONLY

/* Now generate the function definitions if we are using a single compilation unit. */
#if GMX_CUDA_NB_SINGLE_COMPILATION_UNIT
#    include "nbfe_cuda_kernel_F.cu"
#    include "nbfe_cuda_kernel_VF.cu"
#    include "nbfe_foreign_cuda_kernel_V.cu"
#endif /* GMX_CUDA_NB_SINGLE_COMPILATION_UNIT */

namespace gmx
{
/*! Nonbonded FEP kernel function pointer type */
typedef void (*nbfe_cu_kfunc_ptr_t)(const NBAtomDataGpu, const NBParamGpu, const GpuFeplist, bool);
/*! Nonbonded foreign FEP kernel function pointer type */
typedef void (*nbfe_foreign_cu_kfunc_ptr_t)(const NBAtomDataGpu, const NBParamGpu, const GpuFeplist, int);
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

/*! Force only fep kernel function pointers. */
static const nbfe_cu_kfunc_ptr_t nb_fep_kfunc_noener_ptr[c_numElecTypes][c_numVdwTypes] = {
    { nbfe_kernel_ElecCut_VdwLJ_F_cuda,
      nbfe_kernel_ElecCut_VdwLJCombGeom_F_cuda,
      nbfe_kernel_ElecCut_VdwLJCombLB_F_cuda,
      nbfe_kernel_ElecCut_VdwLJFsw_F_cuda,
      nbfe_kernel_ElecCut_VdwLJPsw_F_cuda,
      nbfe_kernel_ElecCut_VdwLJEwCombGeom_F_cuda,
      nbfe_kernel_ElecCut_VdwLJEwCombLB_F_cuda },
    { nbfe_kernel_ElecRF_VdwLJ_F_cuda,
      nbfe_kernel_ElecRF_VdwLJCombGeom_F_cuda,
      nbfe_kernel_ElecRF_VdwLJCombLB_F_cuda,
      nbfe_kernel_ElecRF_VdwLJFsw_F_cuda,
      nbfe_kernel_ElecRF_VdwLJPsw_F_cuda,
      nbfe_kernel_ElecRF_VdwLJEwCombGeom_F_cuda,
      nbfe_kernel_ElecRF_VdwLJEwCombLB_F_cuda },
    { nbfe_kernel_ElecEwQSTab_VdwLJ_F_cuda,
      nbfe_kernel_ElecEwQSTab_VdwLJCombGeom_F_cuda,
      nbfe_kernel_ElecEwQSTab_VdwLJCombLB_F_cuda,
      nbfe_kernel_ElecEwQSTab_VdwLJFsw_F_cuda,
      nbfe_kernel_ElecEwQSTab_VdwLJPsw_F_cuda,
      nbfe_kernel_ElecEwQSTab_VdwLJEwCombGeom_F_cuda,
      nbfe_kernel_ElecEwQSTab_VdwLJEwCombLB_F_cuda },
    { nbfe_kernel_ElecEwQSTabTwinCut_VdwLJ_F_cuda,
      nbfe_kernel_ElecEwQSTabTwinCut_VdwLJCombGeom_F_cuda,
      nbfe_kernel_ElecEwQSTabTwinCut_VdwLJCombLB_F_cuda,
      nbfe_kernel_ElecEwQSTabTwinCut_VdwLJFsw_F_cuda,
      nbfe_kernel_ElecEwQSTabTwinCut_VdwLJPsw_F_cuda,
      nbfe_kernel_ElecEwQSTabTwinCut_VdwLJEwCombGeom_F_cuda,
      nbfe_kernel_ElecEwQSTabTwinCut_VdwLJEwCombLB_F_cuda },
    { nbfe_kernel_ElecEw_VdwLJ_F_cuda,
      nbfe_kernel_ElecEw_VdwLJCombGeom_F_cuda,
      nbfe_kernel_ElecEw_VdwLJCombLB_F_cuda,
      nbfe_kernel_ElecEw_VdwLJFsw_F_cuda,
      nbfe_kernel_ElecEw_VdwLJPsw_F_cuda,
      nbfe_kernel_ElecEw_VdwLJEwCombGeom_F_cuda,
      nbfe_kernel_ElecEw_VdwLJEwCombLB_F_cuda },
    { nbfe_kernel_ElecEwTwinCut_VdwLJ_F_cuda,
      nbfe_kernel_ElecEwTwinCut_VdwLJCombGeom_F_cuda,
      nbfe_kernel_ElecEwTwinCut_VdwLJCombLB_F_cuda,
      nbfe_kernel_ElecEwTwinCut_VdwLJFsw_F_cuda,
      nbfe_kernel_ElecEwTwinCut_VdwLJPsw_F_cuda,
      nbfe_kernel_ElecEwTwinCut_VdwLJEwCombGeom_F_cuda,
      nbfe_kernel_ElecEwTwinCut_VdwLJEwCombLB_F_cuda }
};

/*! Force + energy fep kernel function pointers. */
static const nbfe_cu_kfunc_ptr_t nb_fep_kfunc_ener_ptr[c_numElecTypes][c_numVdwTypes] = {
    { nbfe_kernel_ElecCut_VdwLJ_VF_cuda,
      nbfe_kernel_ElecCut_VdwLJCombGeom_VF_cuda,
      nbfe_kernel_ElecCut_VdwLJCombLB_VF_cuda,
      nbfe_kernel_ElecCut_VdwLJFsw_VF_cuda,
      nbfe_kernel_ElecCut_VdwLJPsw_VF_cuda,
      nbfe_kernel_ElecCut_VdwLJEwCombGeom_VF_cuda,
      nbfe_kernel_ElecCut_VdwLJEwCombLB_VF_cuda },
    { nbfe_kernel_ElecRF_VdwLJ_VF_cuda,
      nbfe_kernel_ElecRF_VdwLJCombGeom_VF_cuda,
      nbfe_kernel_ElecRF_VdwLJCombLB_VF_cuda,
      nbfe_kernel_ElecRF_VdwLJFsw_VF_cuda,
      nbfe_kernel_ElecRF_VdwLJPsw_VF_cuda,
      nbfe_kernel_ElecRF_VdwLJEwCombGeom_VF_cuda,
      nbfe_kernel_ElecRF_VdwLJEwCombLB_VF_cuda },
    { nbfe_kernel_ElecEwQSTab_VdwLJ_VF_cuda,
      nbfe_kernel_ElecEwQSTab_VdwLJCombGeom_VF_cuda,
      nbfe_kernel_ElecEwQSTab_VdwLJCombLB_VF_cuda,
      nbfe_kernel_ElecEwQSTab_VdwLJFsw_VF_cuda,
      nbfe_kernel_ElecEwQSTab_VdwLJPsw_VF_cuda,
      nbfe_kernel_ElecEwQSTab_VdwLJEwCombGeom_VF_cuda,
      nbfe_kernel_ElecEwQSTab_VdwLJEwCombLB_VF_cuda },
    { nbfe_kernel_ElecEwQSTabTwinCut_VdwLJ_VF_cuda,
      nbfe_kernel_ElecEwQSTabTwinCut_VdwLJCombGeom_VF_cuda,
      nbfe_kernel_ElecEwQSTabTwinCut_VdwLJCombLB_VF_cuda,
      nbfe_kernel_ElecEwQSTabTwinCut_VdwLJFsw_VF_cuda,
      nbfe_kernel_ElecEwQSTabTwinCut_VdwLJPsw_VF_cuda,
      nbfe_kernel_ElecEwQSTabTwinCut_VdwLJEwCombGeom_VF_cuda,
      nbfe_kernel_ElecEwQSTabTwinCut_VdwLJEwCombLB_VF_cuda },
    { nbfe_kernel_ElecEw_VdwLJ_VF_cuda,
      nbfe_kernel_ElecEw_VdwLJCombGeom_VF_cuda,
      nbfe_kernel_ElecEw_VdwLJCombLB_VF_cuda,
      nbfe_kernel_ElecEw_VdwLJFsw_VF_cuda,
      nbfe_kernel_ElecEw_VdwLJPsw_VF_cuda,
      nbfe_kernel_ElecEw_VdwLJEwCombGeom_VF_cuda,
      nbfe_kernel_ElecEw_VdwLJEwCombLB_VF_cuda },
    { nbfe_kernel_ElecEwTwinCut_VdwLJ_VF_cuda,
      nbfe_kernel_ElecEwTwinCut_VdwLJCombGeom_VF_cuda,
      nbfe_kernel_ElecEwTwinCut_VdwLJCombLB_VF_cuda,
      nbfe_kernel_ElecEwTwinCut_VdwLJFsw_VF_cuda,
      nbfe_kernel_ElecEwTwinCut_VdwLJPsw_VF_cuda,
      nbfe_kernel_ElecEwTwinCut_VdwLJEwCombGeom_VF_cuda,
      nbfe_kernel_ElecEwTwinCut_VdwLJEwCombLB_VF_cuda }
};

/*! Return a pointer to the nb fep kernel version to be executed at the current step. */
static inline nbfe_cu_kfunc_ptr_t select_nbfe_kernel(enum ElecType elecType, enum VdwType vdwType, bool bDoEne)
{
    const int elecTypeIdx = static_cast<int>(elecType);
    const int vdwTypeIdx  = static_cast<int>(vdwType);

    GMX_ASSERT(elecTypeIdx < c_numElecTypes,
               "The electrostatics type requested is not implemented in the CUDA kernels.");
    GMX_ASSERT(vdwTypeIdx < c_numVdwTypes,
               "The VdW type requested is not implemented in the CUDA kernels.");

    if (bDoEne)
    {
        return nb_fep_kfunc_ener_ptr[elecTypeIdx][vdwTypeIdx];
    }
    else
    {
        return nb_fep_kfunc_noener_ptr[elecTypeIdx][vdwTypeIdx];
    }
}

/*! foreign lambda kernel function pointers. */
static const nbfe_foreign_cu_kfunc_ptr_t nb_foreign_fep_kfunc_ptr[c_numElecTypes][c_numVdwTypes] = {
    { nbfe_foreign_kernel_ElecCut_VdwLJ_V_cuda,
      nbfe_foreign_kernel_ElecCut_VdwLJCombGeom_V_cuda,
      nbfe_foreign_kernel_ElecCut_VdwLJCombLB_V_cuda,
      nbfe_foreign_kernel_ElecCut_VdwLJFsw_V_cuda,
      nbfe_foreign_kernel_ElecCut_VdwLJPsw_V_cuda,
      nbfe_foreign_kernel_ElecCut_VdwLJEwCombGeom_V_cuda,
      nbfe_foreign_kernel_ElecCut_VdwLJEwCombLB_V_cuda },
    { nbfe_foreign_kernel_ElecRF_VdwLJ_V_cuda,
      nbfe_foreign_kernel_ElecRF_VdwLJCombGeom_V_cuda,
      nbfe_foreign_kernel_ElecRF_VdwLJCombLB_V_cuda,
      nbfe_foreign_kernel_ElecRF_VdwLJFsw_V_cuda,
      nbfe_foreign_kernel_ElecRF_VdwLJPsw_V_cuda,
      nbfe_foreign_kernel_ElecRF_VdwLJEwCombGeom_V_cuda,
      nbfe_foreign_kernel_ElecRF_VdwLJEwCombLB_V_cuda },
    { nbfe_foreign_kernel_ElecEwQSTab_VdwLJ_V_cuda,
      nbfe_foreign_kernel_ElecEwQSTab_VdwLJCombGeom_V_cuda,
      nbfe_foreign_kernel_ElecEwQSTab_VdwLJCombLB_V_cuda,
      nbfe_foreign_kernel_ElecEwQSTab_VdwLJFsw_V_cuda,
      nbfe_foreign_kernel_ElecEwQSTab_VdwLJPsw_V_cuda,
      nbfe_foreign_kernel_ElecEwQSTab_VdwLJEwCombGeom_V_cuda,
      nbfe_foreign_kernel_ElecEwQSTab_VdwLJEwCombLB_V_cuda },
    { nbfe_foreign_kernel_ElecEwQSTabTwinCut_VdwLJ_V_cuda,
      nbfe_foreign_kernel_ElecEwQSTabTwinCut_VdwLJCombGeom_V_cuda,
      nbfe_foreign_kernel_ElecEwQSTabTwinCut_VdwLJCombLB_V_cuda,
      nbfe_foreign_kernel_ElecEwQSTabTwinCut_VdwLJFsw_V_cuda,
      nbfe_foreign_kernel_ElecEwQSTabTwinCut_VdwLJPsw_V_cuda,
      nbfe_foreign_kernel_ElecEwQSTabTwinCut_VdwLJEwCombGeom_V_cuda,
      nbfe_foreign_kernel_ElecEwQSTabTwinCut_VdwLJEwCombLB_V_cuda },
    { nbfe_foreign_kernel_ElecEw_VdwLJ_V_cuda,
      nbfe_foreign_kernel_ElecEw_VdwLJCombGeom_V_cuda,
      nbfe_foreign_kernel_ElecEw_VdwLJCombLB_V_cuda,
      nbfe_foreign_kernel_ElecEw_VdwLJFsw_V_cuda,
      nbfe_foreign_kernel_ElecEw_VdwLJPsw_V_cuda,
      nbfe_foreign_kernel_ElecEw_VdwLJEwCombGeom_V_cuda,
      nbfe_foreign_kernel_ElecEw_VdwLJEwCombLB_V_cuda },
    { nbfe_foreign_kernel_ElecEwTwinCut_VdwLJ_V_cuda,
      nbfe_foreign_kernel_ElecEwTwinCut_VdwLJCombGeom_V_cuda,
      nbfe_foreign_kernel_ElecEwTwinCut_VdwLJCombLB_V_cuda,
      nbfe_foreign_kernel_ElecEwTwinCut_VdwLJFsw_V_cuda,
      nbfe_foreign_kernel_ElecEwTwinCut_VdwLJPsw_V_cuda,
      nbfe_foreign_kernel_ElecEwTwinCut_VdwLJEwCombGeom_V_cuda,
      nbfe_foreign_kernel_ElecEwTwinCut_VdwLJEwCombLB_V_cuda }
};

/*! Return a pointer to the foreign lambda kernel version to be executed at the current step. */
static inline nbfe_foreign_cu_kfunc_ptr_t select_nbfe_foreign_kernel(enum ElecType elecType,
                                                                     enum VdwType  vdwType,
                                                                     const DeviceInformation gmx_unused* deviceInfo)
{
    const int elecTypeIdx = static_cast<int>(elecType);
    const int vdwTypeIdx  = static_cast<int>(vdwType);

    GMX_ASSERT(elecTypeIdx < c_numElecTypes,
               "The electrostatics type requested is not implemented in the CUDA kernels.");
    GMX_ASSERT(vdwTypeIdx < c_numVdwTypes,
               "The VdW type requested is not implemented in the CUDA kernels.");

    return nb_foreign_fep_kfunc_ptr[elecTypeIdx][vdwTypeIdx];
}

/*! Launch the Nonbonded free energy GPU kernels. */
void gpu_launch_free_energy_kernel(NbnxmGpu*                      nb,
                                   const gmx::SimulationWorkload& simulationWork,
                                   const gmx::StepWorkload&       stepWork,
                                   const InteractionLocality      iloc)
{
    NBAtomDataGpu*      adat         = nb->atdat;
    NBParamGpu*         nbp          = nb->nbparam;
    auto*               feplist      = nb->feplist[iloc].get();
    const DeviceStream& deviceStream = *nb->deviceStreams[iloc];

    KernelLaunchConfig fepConfig;
    fepConfig.blockSize[0] = 64;
    fepConfig.blockSize[1] = 1;
    fepConfig.blockSize[2] = 1;

    const int bsize = fepConfig.blockSize[0] * fepConfig.blockSize[1] * fepConfig.blockSize[2];

    // one warp per nri
    const int nriPerBlock      = bsize / warp_size;
    int       nblock           = (feplist->numiAtoms + nriPerBlock - 1) / nriPerBlock;
    nblock                     = calc_nb_kernel_nblock(nblock, &nb->deviceContext_->deviceInfo());
    fepConfig.gridSize[0]      = nblock;
    fepConfig.gridSize[1]      = 1;
    fepConfig.gridSize[2]      = 1;
    fepConfig.sharedMemorySize = 0;

    if (debug)
    {
        fprintf(debug,
                "Non-bonded FEP GPU launch configuration:\n\tThread block: %zux%zux%zu\n\t"
                "\tGrid: %zux%zu\n\t#FEP nbl numiAtoms: %d \n",
                fepConfig.blockSize[0],
                fepConfig.blockSize[1],
                fepConfig.blockSize[2],
                fepConfig.gridSize[0],
                fepConfig.gridSize[1],
                feplist->numiAtoms);
    }
    const auto fepKernel = select_nbfe_kernel(nbp->elecType, nbp->vdwType, stepWork.computeEnergy);
    const auto fepKernelArgs = prepareGpuKernelArguments(
            fepKernel, fepConfig, adat, nbp, feplist, &stepWork.computeVirial);

    // Launch the FEP kernel
    launchGpuKernel(fepKernel, fepConfig, deviceStream, nullptr, "k_calc_nb_fep", fepKernelArgs);

    // Launch the foreign lambdas kernel
    if (simulationWork.useGpuForeignNonbondedFE && stepWork.computeDhdl)
    {
        const int          nLambda = nb->fephostdata->allLambdaCoul.size();
        KernelLaunchConfig fepForeignConfig;
        fepForeignConfig.blockSize[0] = 64;
        fepForeignConfig.blockSize[1] = fepConfig.blockSize[1];
        fepForeignConfig.blockSize[2] = fepConfig.blockSize[2];
        const int blockSize = fepForeignConfig.blockSize[0] * fepForeignConfig.blockSize[1]
                              * fepForeignConfig.blockSize[2];

        // one warp per nri
        const int nriPerBlock = blockSize / warp_size;
        nblock                = (feplist->numiAtoms + nriPerBlock - 1) / nriPerBlock;
        nblock                = calc_nb_kernel_nblock(nblock, &nb->deviceContext_->deviceInfo());
        fepForeignConfig.gridSize[0]      = nblock;
        fepForeignConfig.gridSize[1]      = 1;
        fepForeignConfig.gridSize[2]      = 1;
        fepForeignConfig.sharedMemorySize = (nLambda + 1) * 2 * sizeof(float);

        if (debug)
        {
            fprintf(debug,
                    "Foreign FEP GPU launch configuration:\n\tThread block: %zux%zux%zu\n\t"
                    "\tGrid: %zux%zu\n\t#Number of lambdas : %d \n"
                    "\tShMem: %zu\n",
                    fepForeignConfig.blockSize[0],
                    fepForeignConfig.blockSize[1],
                    fepForeignConfig.blockSize[2],
                    fepForeignConfig.gridSize[0],
                    fepForeignConfig.gridSize[1],
                    nLambda,
                    fepForeignConfig.sharedMemorySize);
        }

        const auto fepForeignKernel = select_nbfe_foreign_kernel(
                nbp->elecType, nbp->vdwType, &nb->deviceContext_->deviceInfo());
        const auto fepForeignKernelArgs = prepareGpuKernelArguments(
                fepForeignKernel, fepForeignConfig, adat, nbp, feplist, &nLambda);

        launchGpuKernel(
                fepForeignKernel, fepForeignConfig, deviceStream, nullptr, "k_calc_nb_fep_foreign", fepForeignKernelArgs);
    }

    if (GMX_NATIVE_WINDOWS)
    {
        /* Windows: force flushing WDDM queue */
        cudaStreamQuery(deviceStream.stream());
    }
}

} // namespace gmx
