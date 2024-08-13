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
/*! \internal \file
 *  \brief Define OpenCL implementation of nbnxm_gpu.h
 *
 *  \author Anca Hamuraru <anca@streamcomputing.eu>
 *  \author Teemu Virolainen <teemu@streamcomputing.eu>
 *  \author Dimitrios Karkoulis <dimitris.karkoulis@gmail.com>
 *  \author Szilárd Páll <pall.szilard@gmail.com>
 *  \ingroup module_nbnxm
 *
 *  TODO (psz):
 *  - Add a static const cl_uint c_pruneKernelWorkDim / c_nbnxnKernelWorkDim = 3;
 *  - Rework the copying of OCL data structures done before every invocation of both
 *    nb and prune kernels (using fillin_ocl_structures); also consider at the same
 *    time calling clSetKernelArg only on the updated parameters (if tracking changed
 *    parameters is feasible);
 *  - Consider using the event_wait_list argument to clEnqueueNDRangeKernel to mark
 *    dependencies on the kernel launched: e.g. the non-local nb kernel's dependency
 *    on the misc_ops_and_local_H2D_done event could be better expressed this way.
 *
 *  - Consider extracting common sections of the OpenCL and CUDA nbnxn logic, e.g:
 *    - in nbnxn_gpu_launch_kernel_pruneonly() the pre- and post-kernel launch logic
 *      is identical in the two implementations, so a 3-way split might allow sharing
 *      code;
 *    -
 *
 */
#include "gmxpre.h"

#include <cassert>
#include <climits>
#include <cstdlib>

#if defined(_MSVC)
#    include <limits>
#endif

#include "gromacs/gpu_utils/device_context.h"
#include "gromacs/gpu_utils/gputraits_ocl.h"
#include "gromacs/gpu_utils/oclutils.h"
#include "gromacs/hardware/device_information.h"
#include "gromacs/hardware/hw_info.h"
#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/nbnxm/atomdata.h"
#include "gromacs/nbnxm/gpu_common.h"
#include "gromacs/nbnxm/gpu_common_utils.h"
#include "gromacs/nbnxm/gpu_data_mgmt.h"
#include "gromacs/nbnxm/nbnxm.h"
#include "gromacs/nbnxm/nbnxm_gpu.h"
#include "gromacs/nbnxm/pairlist.h"
#include "gromacs/timing/gpu_timing.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"

#include "nbnxm_ocl_types.h"

namespace gmx
{

/*! \brief Validates the input global work size parameter.
 */
static inline void validate_global_work_size(const KernelLaunchConfig& config,
                                             int                       work_dim,
                                             const DeviceInformation*  dinfo)
{
    cl_uint device_size_t_size_bits;
    cl_uint host_size_t_size_bits;

    GMX_ASSERT(dinfo, "Need a valid device info object");

    size_t global_work_size[3];
    GMX_ASSERT(work_dim <= 3, "Not supporting hyper-grids just yet");
    for (int i = 0; i < work_dim; i++)
    {
        global_work_size[i] = config.blockSize[i] * config.gridSize[i];
    }

    /* Each component of a global_work_size must not exceed the range given by the
       sizeof(device size_t) for the device on which the kernel execution will
       be enqueued. See:
       https://www.khronos.org/registry/cl/sdk/1.0/docs/man/xhtml/clEnqueueNDRangeKernel.html
     */
    device_size_t_size_bits = dinfo->adress_bits;
    host_size_t_size_bits   = static_cast<cl_uint>(sizeof(size_t) * CHAR_BIT);

    /* If sizeof(host size_t) <= sizeof(device size_t)
            => global_work_size components will always be valid
       else
            => get device limit for global work size and
            compare it against each component of global_work_size.
     */
    if (host_size_t_size_bits > device_size_t_size_bits)
    {
        size_t device_limit;

        device_limit = (1ULL << device_size_t_size_bits) - 1;

        for (int i = 0; i < work_dim; i++)
        {
            if (global_work_size[i] > device_limit)
            {
                gmx_fatal(
                        FARGS,
                        "Watch out, the input system is too large to simulate!\n"
                        "The number of nonbonded work units (=number of super-clusters) exceeds the"
                        "device capabilities. Global work size limit exceeded (%zu > %zu)!",
                        global_work_size[i],
                        device_limit);
            }
        }
    }
}

/* Constant arrays listing non-bonded kernel function names. The arrays are
 * organized in 2-dim arrays by: electrostatics and VDW type.
 *
 *  Note that the row- and column-order of function pointers has to match the
 *  order of corresponding enumerated electrostatics and vdw types, resp.,
 *  defined in nbnxm_ocl_types.h.
 */

/*! \brief Force-only kernel function names. */
static const char* nb_kfunc_noener_noprune_ptr[c_numElecTypes][c_numVdwTypes] = {
    { "nbnxn_kernel_ElecCut_VdwLJ_F_opencl",
      "nbnxn_kernel_ElecCut_VdwLJCombGeom_F_opencl",
      "nbnxn_kernel_ElecCut_VdwLJCombLB_F_opencl",
      "nbnxn_kernel_ElecCut_VdwLJFsw_F_opencl",
      "nbnxn_kernel_ElecCut_VdwLJPsw_F_opencl",
      "nbnxn_kernel_ElecCut_VdwLJEwCombGeom_F_opencl",
      "nbnxn_kernel_ElecCut_VdwLJEwCombLB_F_opencl" },
    { "nbnxn_kernel_ElecRF_VdwLJ_F_opencl",
      "nbnxn_kernel_ElecRF_VdwLJCombGeom_F_opencl",
      "nbnxn_kernel_ElecRF_VdwLJCombLB_F_opencl",
      "nbnxn_kernel_ElecRF_VdwLJFsw_F_opencl",
      "nbnxn_kernel_ElecRF_VdwLJPsw_F_opencl",
      "nbnxn_kernel_ElecRF_VdwLJEwCombGeom_F_opencl",
      "nbnxn_kernel_ElecRF_VdwLJEwCombLB_F_opencl" },
    { "nbnxn_kernel_ElecEwQSTab_VdwLJ_F_opencl",
      "nbnxn_kernel_ElecEwQSTab_VdwLJCombGeom_F_opencl",
      "nbnxn_kernel_ElecEwQSTab_VdwLJCombLB_F_opencl",
      "nbnxn_kernel_ElecEwQSTab_VdwLJFsw_F_opencl",
      "nbnxn_kernel_ElecEwQSTab_VdwLJPsw_F_opencl",
      "nbnxn_kernel_ElecEwQSTab_VdwLJEwCombGeom_F_opencl",
      "nbnxn_kernel_ElecEwQSTab_VdwLJEwCombLB_F_opencl" },
    { "nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJ_F_opencl",
      "nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJCombGeom_F_opencl",
      "nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJCombLB_F_opencl",
      "nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJFsw_F_opencl",
      "nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJPsw_F_opencl",
      "nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJEwCombGeom_F_opencl",
      "nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJEwCombLB_F_opencl" },
    { "nbnxn_kernel_ElecEw_VdwLJ_F_opencl",
      "nbnxn_kernel_ElecEw_VdwLJCombGeom_F_opencl",
      "nbnxn_kernel_ElecEw_VdwLJCombLB_F_opencl",
      "nbnxn_kernel_ElecEw_VdwLJFsw_F_opencl",
      "nbnxn_kernel_ElecEw_VdwLJPsw_F_opencl",
      "nbnxn_kernel_ElecEw_VdwLJEwCombGeom_F_opencl",
      "nbnxn_kernel_ElecEw_VdwLJEwCombLB_F_opencl" },
    { "nbnxn_kernel_ElecEwTwinCut_VdwLJ_F_opencl",
      "nbnxn_kernel_ElecEwTwinCut_VdwLJCombGeom_F_opencl",
      "nbnxn_kernel_ElecEwTwinCut_VdwLJCombLB_F_opencl",
      "nbnxn_kernel_ElecEwTwinCut_VdwLJFsw_F_opencl",
      "nbnxn_kernel_ElecEwTwinCut_VdwLJPsw_F_opencl",
      "nbnxn_kernel_ElecEwTwinCut_VdwLJEwCombGeom_F_opencl",
      "nbnxn_kernel_ElecEwTwinCut_VdwLJEwCombLB_F_opencl" }
};

/*! \brief Force + energy kernel function pointers. */
static const char* nb_kfunc_ener_noprune_ptr[c_numElecTypes][c_numVdwTypes] = {
    { "nbnxn_kernel_ElecCut_VdwLJ_VF_opencl",
      "nbnxn_kernel_ElecCut_VdwLJCombGeom_VF_opencl",
      "nbnxn_kernel_ElecCut_VdwLJCombLB_VF_opencl",
      "nbnxn_kernel_ElecCut_VdwLJFsw_VF_opencl",
      "nbnxn_kernel_ElecCut_VdwLJPsw_VF_opencl",
      "nbnxn_kernel_ElecCut_VdwLJEwCombGeom_VF_opencl",
      "nbnxn_kernel_ElecCut_VdwLJEwCombLB_VF_opencl" },
    { "nbnxn_kernel_ElecRF_VdwLJ_VF_opencl",
      "nbnxn_kernel_ElecRF_VdwLJCombGeom_VF_opencl",
      "nbnxn_kernel_ElecRF_VdwLJCombLB_VF_opencl",
      "nbnxn_kernel_ElecRF_VdwLJFsw_VF_opencl",
      "nbnxn_kernel_ElecRF_VdwLJPsw_VF_opencl",
      "nbnxn_kernel_ElecRF_VdwLJEwCombGeom_VF_opencl",
      "nbnxn_kernel_ElecRF_VdwLJEwCombLB_VF_opencl" },
    { "nbnxn_kernel_ElecEwQSTab_VdwLJ_VF_opencl",
      "nbnxn_kernel_ElecEwQSTab_VdwLJCombGeom_VF_opencl",
      "nbnxn_kernel_ElecEwQSTab_VdwLJCombLB_VF_opencl",
      "nbnxn_kernel_ElecEwQSTab_VdwLJFsw_VF_opencl",
      "nbnxn_kernel_ElecEwQSTab_VdwLJPsw_VF_opencl",
      "nbnxn_kernel_ElecEwQSTab_VdwLJEwCombGeom_VF_opencl",
      "nbnxn_kernel_ElecEwQSTab_VdwLJEwCombLB_VF_opencl" },
    { "nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJ_VF_opencl",
      "nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJCombGeom_VF_opencl",
      "nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJCombLB_VF_opencl",
      "nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJFsw_VF_opencl",
      "nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJPsw_VF_opencl",
      "nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJEwCombGeom_VF_opencl",
      "nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJEwCombLB_VF_opencl" },
    { "nbnxn_kernel_ElecEw_VdwLJ_VF_opencl",
      "nbnxn_kernel_ElecEw_VdwLJCombGeom_VF_opencl",
      "nbnxn_kernel_ElecEw_VdwLJCombLB_VF_opencl",
      "nbnxn_kernel_ElecEw_VdwLJFsw_VF_opencl",
      "nbnxn_kernel_ElecEw_VdwLJPsw_VF_opencl",
      "nbnxn_kernel_ElecEw_VdwLJEwCombGeom_VF_opencl",
      "nbnxn_kernel_ElecEw_VdwLJEwCombLB_VF_opencl" },
    { "nbnxn_kernel_ElecEwTwinCut_VdwLJ_VF_opencl",
      "nbnxn_kernel_ElecEwTwinCut_VdwLJCombGeom_VF_opencl",
      "nbnxn_kernel_ElecEwTwinCut_VdwLJCombLB_VF_opencl",
      "nbnxn_kernel_ElecEwTwinCut_VdwLJFsw_VF_opencl",
      "nbnxn_kernel_ElecEwTwinCut_VdwLJPsw_VF_opencl",
      "nbnxn_kernel_ElecEwTwinCut_VdwLJEwCombGeom_VF_opencl",
      "nbnxn_kernel_ElecEwTwinCut_VdwLJEwCombLB_VF_opencl" }
};

/*! \brief Force + pruning kernel function pointers. */
static const char* nb_kfunc_noener_prune_ptr[c_numElecTypes][c_numVdwTypes] = {
    { "nbnxn_kernel_ElecCut_VdwLJ_F_prune_opencl",
      "nbnxn_kernel_ElecCut_VdwLJCombGeom_F_prune_opencl",
      "nbnxn_kernel_ElecCut_VdwLJCombLB_F_prune_opencl",
      "nbnxn_kernel_ElecCut_VdwLJFsw_F_prune_opencl",
      "nbnxn_kernel_ElecCut_VdwLJPsw_F_prune_opencl",
      "nbnxn_kernel_ElecCut_VdwLJEwCombGeom_F_prune_opencl",
      "nbnxn_kernel_ElecCut_VdwLJEwCombLB_F_prune_opencl" },
    { "nbnxn_kernel_ElecRF_VdwLJ_F_prune_opencl",
      "nbnxn_kernel_ElecRF_VdwLJCombGeom_F_prune_opencl",
      "nbnxn_kernel_ElecRF_VdwLJCombLB_F_prune_opencl",
      "nbnxn_kernel_ElecRF_VdwLJFsw_F_prune_opencl",
      "nbnxn_kernel_ElecRF_VdwLJPsw_F_prune_opencl",
      "nbnxn_kernel_ElecRF_VdwLJEwCombGeom_F_prune_opencl",
      "nbnxn_kernel_ElecRF_VdwLJEwCombLB_F_prune_opencl" },
    { "nbnxn_kernel_ElecEwQSTab_VdwLJ_F_prune_opencl",
      "nbnxn_kernel_ElecEwQSTab_VdwLJCombGeom_F_prune_opencl",
      "nbnxn_kernel_ElecEwQSTab_VdwLJCombLB_F_prune_opencl",
      "nbnxn_kernel_ElecEwQSTab_VdwLJFsw_F_prune_opencl",
      "nbnxn_kernel_ElecEwQSTab_VdwLJPsw_F_prune_opencl",
      "nbnxn_kernel_ElecEwQSTab_VdwLJEwCombGeom_F_prune_opencl",
      "nbnxn_kernel_ElecEwQSTab_VdwLJEwCombLB_F_prune_opencl" },
    { "nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJ_F_prune_opencl",
      "nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJCombGeom_F_prune_opencl",
      "nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJCombLB_F_prune_opencl",
      "nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJFsw_F_prune_opencl",
      "nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJPsw_F_prune_opencl",
      "nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJEwCombGeom_F_prune_opencl",
      "nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJEwCombLB_F_prune_opencl" },
    { "nbnxn_kernel_ElecEw_VdwLJ_F_prune_opencl",
      "nbnxn_kernel_ElecEw_VdwLJCombGeom_F_prune_opencl",
      "nbnxn_kernel_ElecEw_VdwLJCombLB_F_prune_opencl",
      "nbnxn_kernel_ElecEw_VdwLJFsw_F_prune_opencl",
      "nbnxn_kernel_ElecEw_VdwLJPsw_F_prune_opencl",
      "nbnxn_kernel_ElecEw_VdwLJEwCombGeom_F_prune_opencl",
      "nbnxn_kernel_ElecEw_VdwLJEwCombLB_F_prune_opencl" },
    { "nbnxn_kernel_ElecEwTwinCut_VdwLJ_F_prune_opencl",
      "nbnxn_kernel_ElecEwTwinCut_VdwLJCombGeom_F_prune_opencl",
      "nbnxn_kernel_ElecEwTwinCut_VdwLJCombLB_F_prune_opencl",
      "nbnxn_kernel_ElecEwTwinCut_VdwLJFsw_F_prune_opencl",
      "nbnxn_kernel_ElecEwTwinCut_VdwLJPsw_F_prune_opencl",
      "nbnxn_kernel_ElecEwTwinCut_VdwLJEwCombGeom_F_prune_opencl",
      "nbnxn_kernel_ElecEwTwinCut_VdwLJEwCombLB_F_prune_opencl" }
};

/*! \brief Force + energy + pruning kernel function pointers. */
static const char* nb_kfunc_ener_prune_ptr[c_numElecTypes][c_numVdwTypes] = {
    { "nbnxn_kernel_ElecCut_VdwLJ_VF_prune_opencl",
      "nbnxn_kernel_ElecCut_VdwLJCombGeom_VF_prune_opencl",
      "nbnxn_kernel_ElecCut_VdwLJCombLB_VF_prune_opencl",
      "nbnxn_kernel_ElecCut_VdwLJFsw_VF_prune_opencl",
      "nbnxn_kernel_ElecCut_VdwLJPsw_VF_prune_opencl",
      "nbnxn_kernel_ElecCut_VdwLJEwCombGeom_VF_prune_opencl",
      "nbnxn_kernel_ElecCut_VdwLJEwCombLB_VF_prune_opencl" },
    { "nbnxn_kernel_ElecRF_VdwLJ_VF_prune_opencl",
      "nbnxn_kernel_ElecRF_VdwLJCombGeom_VF_prune_opencl",
      "nbnxn_kernel_ElecRF_VdwLJCombLB_VF_prune_opencl",
      "nbnxn_kernel_ElecRF_VdwLJFsw_VF_prune_opencl",
      "nbnxn_kernel_ElecRF_VdwLJPsw_VF_prune_opencl",
      "nbnxn_kernel_ElecRF_VdwLJEwCombGeom_VF_prune_opencl",
      "nbnxn_kernel_ElecRF_VdwLJEwCombLB_VF_prune_opencl" },
    { "nbnxn_kernel_ElecEwQSTab_VdwLJ_VF_prune_opencl",
      "nbnxn_kernel_ElecEwQSTab_VdwLJCombGeom_VF_prune_opencl",
      "nbnxn_kernel_ElecEwQSTab_VdwLJCombLB_VF_prune_opencl",
      "nbnxn_kernel_ElecEwQSTab_VdwLJFsw_VF_prune_opencl",
      "nbnxn_kernel_ElecEwQSTab_VdwLJPsw_VF_prune_opencl",
      "nbnxn_kernel_ElecEwQSTab_VdwLJEwCombGeom_VF_prune_opencl",
      "nbnxn_kernel_ElecEwQSTab_VdwLJEwCombLB_VF_prune_opencl" },
    { "nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJ_VF_prune_opencl",
      "nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJCombGeom_VF_prune_opencl",
      "nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJCombLB_VF_prune_opencl",
      "nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJFsw_VF_prune_opencl",
      "nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJPsw_VF_prune_opencl",
      "nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJEwCombGeom_VF_prune_opencl",
      "nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJEwCombLB_VF_prune_opencl" },
    { "nbnxn_kernel_ElecEw_VdwLJ_VF_prune_opencl",
      "nbnxn_kernel_ElecEw_VdwLJCombGeom_VF_prune_opencl",
      "nbnxn_kernel_ElecEw_VdwLJCombLB_VF_prune_opencl",
      "nbnxn_kernel_ElecEw_VdwLJFsw_VF_prune_opencl",
      "nbnxn_kernel_ElecEw_VdwLJPsw_VF_prune_opencl",
      "nbnxn_kernel_ElecEw_VdwLJEwCombGeom_VF_prune_opencl",
      "nbnxn_kernel_ElecEw_VdwLJEwCombLB_VF_prune_opencl" },
    { "nbnxn_kernel_ElecEwTwinCut_VdwLJ_VF_prune_opencl",
      "nbnxn_kernel_ElecEwTwinCut_VdwLJCombGeom_VF_prune_opencl",
      "nbnxn_kernel_ElecEwTwinCut_VdwLJCombLB_VF_prune_opencl",
      "nbnxn_kernel_ElecEwTwinCut_VdwLJFsw_VF_prune_opencl",
      "nbnxn_kernel_ElecEwTwinCut_VdwLJPsw_VF_prune_opencl",
      "nbnxn_kernel_ElecEwTwinCut_VdwLJEwCombGeom_VF_prune_opencl",
      "nbnxn_kernel_ElecEwTwinCut_VdwLJEwCombLB_VF_prune_opencl" }
};

/*! \brief Return a pointer to the prune kernel version to be executed at the current invocation.
 *
 * \param[in] kernel_pruneonly  array of prune kernel objects
 * \param[in] firstPrunePass    true if the first pruning pass is being executed
 */
static inline cl_kernel selectPruneKernel(cl_kernel kernel_pruneonly[], bool firstPrunePass)
{
    cl_kernel* kernelPtr;

    if (firstPrunePass)
    {
        kernelPtr = &(kernel_pruneonly[epruneFirst]);
    }
    else
    {
        kernelPtr = &(kernel_pruneonly[epruneRolling]);
    }
    // TODO: consider creating the prune kernel object here to avoid a
    // clCreateKernel for the rolling prune kernel if this is not needed.
    return *kernelPtr;
}

/*! \brief Return a pointer to the kernel version to be executed at the current step.
 *  OpenCL kernel objects are cached in nb. If the requested kernel is not
 *  found in the cache, it will be created and the cache will be updated.
 */
static inline cl_kernel
select_nbnxn_kernel(NbnxmGpu* nb, enum ElecType elecType, enum VdwType vdwType, bool bDoEne, bool bDoPrune)
{
    const char* kernel_name_to_run;
    cl_kernel*  kernel_ptr;
    cl_int      cl_error;

    const int elecTypeIdx = static_cast<int>(elecType);
    const int vdwTypeIdx  = static_cast<int>(vdwType);

    GMX_ASSERT(elecTypeIdx < c_numElecTypes,
               "The electrostatics type requested is not implemented in the OpenCL kernels.");
    GMX_ASSERT(vdwTypeIdx < c_numVdwTypes,
               "The VdW type requested is not implemented in the OpenCL kernels.");

    if (bDoEne)
    {
        if (bDoPrune)
        {
            kernel_name_to_run = nb_kfunc_ener_prune_ptr[elecTypeIdx][vdwTypeIdx];
            kernel_ptr         = &(nb->kernel_ener_prune_ptr[elecTypeIdx][vdwTypeIdx]);
        }
        else
        {
            kernel_name_to_run = nb_kfunc_ener_noprune_ptr[elecTypeIdx][vdwTypeIdx];
            kernel_ptr         = &(nb->kernel_ener_noprune_ptr[elecTypeIdx][vdwTypeIdx]);
        }
    }
    else
    {
        if (bDoPrune)
        {
            kernel_name_to_run = nb_kfunc_noener_prune_ptr[elecTypeIdx][vdwTypeIdx];
            kernel_ptr         = &(nb->kernel_noener_prune_ptr[elecTypeIdx][vdwTypeIdx]);
        }
        else
        {
            kernel_name_to_run = nb_kfunc_noener_noprune_ptr[elecTypeIdx][vdwTypeIdx];
            kernel_ptr         = &(nb->kernel_noener_noprune_ptr[elecTypeIdx][vdwTypeIdx]);
        }
    }

    if (nullptr == kernel_ptr[0])
    {
        *kernel_ptr = clCreateKernel(nb->dev_rundata->program, kernel_name_to_run, &cl_error);
        GMX_ASSERT(cl_error == CL_SUCCESS,
                   ("clCreateKernel failed: " + ocl_get_error_string(cl_error)
                    + " for kernel named " + kernel_name_to_run)
                           .c_str());
    }

    return *kernel_ptr;
}

/*! \brief Calculates the amount of shared memory required by the nonbonded kernel in use.
 */
static inline int calc_shmem_required_nonbonded(enum VdwType vdwType, bool bPrefetchLjParam)
{
    int shmem;

    /* size of shmem (force-buffers/xq/atom type preloading) */
    /* NOTE: with the default kernel on sm3.0 we need shmem only for pre-loading */
    /* i-atom x+q in shared memory */
    shmem = c_nbnxnGpuNumClusterPerSupercluster * c_clSize * sizeof(float) * 4; /* xqib */
    /* cj in shared memory, for both warps separately
     * TODO: in the "nowarp kernels we load cj only once  so the factor 2 is not needed.
     */
    shmem += 2 * c_nbnxnGpuJgroupSize * sizeof(int); /* cjs  */
    if (bPrefetchLjParam)
    {
        if (useLjCombRule(vdwType))
        {
            /* i-atom LJ combination parameters in shared memory */
            shmem += c_nbnxnGpuNumClusterPerSupercluster * c_clSize * 2
                     * sizeof(float); /* atib abused for ljcp, float2 */
        }
        else
        {
            /* i-atom types in shared memory */
            shmem += c_nbnxnGpuNumClusterPerSupercluster * c_clSize * sizeof(int); /* atib */
        }
    }
    /* force reduction buffers in shared memory */
    shmem += c_clSize * c_clSize * 3 * sizeof(float); /* f_buf */
    /* Warp vote. In fact it must be * number of warps in block.. */
    shmem += sizeof(cl_uint) * 2; /* warp_any */
    return shmem;
}

/*! \brief Initializes data structures that are going to be sent to the OpenCL device.
 *
 *  The device can't use the same data structures as the host for two main reasons:
 *  - OpenCL restrictions (pointers are not accepted inside data structures)
 *  - some host side fields are not needed for the OpenCL kernels.
 *
 *  This function is called before the launch of both nbnxn and prune kernels.
 */
static void fillin_ocl_structures(NBParamGpu* nbp, cl_nbparam_params_t* nbparams_params)
{
    nbparams_params->coulomb_tab_scale = nbp->coulomb_tab_scale;
    nbparams_params->c_rf              = nbp->c_rf;
    nbparams_params->dispersion_shift  = nbp->dispersion_shift;
    nbparams_params->elecType          = nbp->elecType;
    nbparams_params->epsfac            = nbp->epsfac;
    nbparams_params->ewaldcoeff_lj     = nbp->ewaldcoeff_lj;
    nbparams_params->ewald_beta        = nbp->ewald_beta;
    nbparams_params->rcoulomb_sq       = nbp->rcoulomb_sq;
    nbparams_params->repulsion_shift   = nbp->repulsion_shift;
    nbparams_params->rlistOuter_sq     = nbp->rlistOuter_sq;
    nbparams_params->rvdw_sq           = nbp->rvdw_sq;
    nbparams_params->rlistInner_sq     = nbp->rlistInner_sq;
    nbparams_params->rvdw_switch       = nbp->rvdw_switch;
    nbparams_params->sh_ewald          = nbp->sh_ewald;
    nbparams_params->sh_lj_ewald       = nbp->sh_lj_ewald;
    nbparams_params->two_k_rf          = nbp->two_k_rf;
    nbparams_params->vdwType           = nbp->vdwType;
    nbparams_params->vdw_switch        = nbp->vdw_switch;
}

/*! \brief Launch GPU kernel

   As we execute nonbonded workload in separate queues, before launching
   the kernel we need to make sure that he following operations have completed:
   - atomdata allocation and related H2D transfers (every nstlist step);
   - pair list H2D transfer (every nstlist step);
   - shift vector H2D transfer (every nstlist step);
   - force (+shift force and energy) output clearing (every step).

   These operations are issued in the local queue at the beginning of the step
   and therefore always complete before the local kernel launch. The non-local
   kernel is launched after the local on the same device/context, so this is
   inherently scheduled after the operations in the local stream (including the
   above "misc_ops").
   However, for the sake of having a future-proof implementation, we use the
   misc_ops_done event to record the point in time when the above  operations
   are finished and synchronize with this event in the non-local stream.
 */
void gpu_launch_kernel(NbnxmGpu* nb, const gmx::StepWorkload& stepWork, const InteractionLocality iloc)
{
    NBAtomDataGpu*      adat         = nb->atdat;
    NBParamGpu*         nbp          = nb->nbparam;
    auto*               plist        = nb->plist[iloc].get();
    GpuTimers*          timers       = nb->timers;
    const DeviceStream& deviceStream = *nb->deviceStreams[iloc];

    bool bDoTime = nb->bDoTime;

    cl_nbparam_params_t nbparams_params;

    /* Don't launch the non-local kernel if there is no work to do.
       Doing the same for the local kernel is more complicated, since the
       local part of the force array also depends on the non-local kernel.
       So to avoid complicating the code and to reduce the risk of bugs,
       we always call the local kernel and later (not in
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
           (that's the way the timing accounting can distinguish between
           separate prune kernel and combined force+prune).
         */
        gpu_launch_kernel_pruneonly(nb, iloc, 1);
    }

    if (plist->numSci == 0)
    {
        /* Don't launch an empty local kernel (is not allowed with OpenCL).
         */
        return;
    }

    /* beginning of timed nonbonded calculation section */
    if (bDoTime)
    {
        timers->interaction[iloc].nb_k.openTimingRegion(deviceStream);
    }

    /* kernel launch config */

    KernelLaunchConfig config;
    config.sharedMemorySize = calc_shmem_required_nonbonded(nbp->vdwType, nb->bPrefetchLjParam);
    config.blockSize[0]     = c_clSize;
    config.blockSize[1]     = c_clSize;
    config.gridSize[0]      = plist->numSci;

    validate_global_work_size(config, 3, &nb->deviceContext_->deviceInfo());

    if (debug)
    {
        fprintf(debug,
                "Non-bonded GPU launch configuration:\n\tLocal work size: %zux%zux%zu\n\t"
                "Global work size : %zux%zu\n\t#Super-clusters/clusters: %d/%d (%d)\n",
                config.blockSize[0],
                config.blockSize[1],
                config.blockSize[2],
                config.blockSize[0] * config.gridSize[0],
                config.blockSize[1] * config.gridSize[1],
                plist->numSci * c_nbnxnGpuNumClusterPerSupercluster,
                c_nbnxnGpuNumClusterPerSupercluster,
                plist->numAtomsPerCluster);
    }

    fillin_ocl_structures(nbp, &nbparams_params);

    auto* timingEvent = bDoTime ? timers->interaction[iloc].nb_k.fetchNextEvent() : nullptr;
    constexpr char kernelName[] = "k_calc_nb";
    const auto     kernel =
            select_nbnxn_kernel(nb,
                                nbp->elecType,
                                nbp->vdwType,
                                stepWork.computeEnergy,
                                (plist->haveFreshList && !nb->timers->interaction[iloc].didPrune));


    // The OpenCL kernel takes int as second to last argument because bool is
    // not supported as a kernel argument type (sizeof(bool) is implementation defined).
    const int computeFshift = static_cast<int>(stepWork.computeVirial);
    if (useLjCombRule(nb->nbparam->vdwType))
    {
        const auto kernelArgs = prepareGpuKernelArguments(kernel,
                                                          config,
                                                          &nbparams_params,
                                                          &adat->xq,
                                                          &adat->f,
                                                          &adat->eLJ,
                                                          &adat->eElec,
                                                          &adat->fShift,
                                                          &adat->ljComb,
                                                          &adat->shiftVec,
                                                          &nbp->nbfp,
                                                          &nbp->nbfp_comb,
                                                          &nbp->coulomb_tab,
                                                          &plist->sci,
                                                          &plist->cjPacked,
                                                          &plist->excl,
                                                          &computeFshift);

        launchGpuKernel(kernel, config, deviceStream, timingEvent, kernelName, kernelArgs);
    }
    else
    {
        const auto kernelArgs = prepareGpuKernelArguments(kernel,
                                                          config,
                                                          &adat->numTypes,
                                                          &nbparams_params,
                                                          &adat->xq,
                                                          &adat->f,
                                                          &adat->eLJ,
                                                          &adat->eElec,
                                                          &adat->fShift,
                                                          &adat->atomTypes,
                                                          &adat->shiftVec,
                                                          &nbp->nbfp,
                                                          &nbp->nbfp_comb,
                                                          &nbp->coulomb_tab,
                                                          &plist->sci,
                                                          &plist->cjPacked,
                                                          &plist->excl,
                                                          &computeFshift);
        launchGpuKernel(kernel, config, deviceStream, timingEvent, kernelName, kernelArgs);
    }

    if (bDoTime)
    {
        timers->interaction[iloc].nb_k.closeTimingRegion(deviceStream);
    }
}


/*! \brief Calculates the amount of shared memory required by the prune kernel.
 *
 *  Note that for the sake of simplicity we use the CUDA terminology "shared memory"
 *  for OpenCL local memory.
 *
 * \param[in] num_threads_z cjPacked concurrency equal to the number of threads/work items in the
 * 3-rd dimension. \returns   the amount of local memory in bytes required by the pruning kernel
 */
static inline int calc_shmem_required_prune(const int num_threads_z)
{
    int shmem;

    /* i-atom x in shared memory (for convenience we load all 4 components including q) */
    shmem = c_nbnxnGpuNumClusterPerSupercluster * c_clSize * sizeof(float) * 4;
    /* cj in shared memory, for each warp separately
     * Note: only need to load once per wavefront, but to keep the code simple,
     * for now we load twice on AMD.
     */
    shmem += num_threads_z * c_nbnxnGpuClusterpairSplit * c_nbnxnGpuJgroupSize * sizeof(int);
    /* Warp vote, requires one uint per warp/32 threads per block. */
    shmem += sizeof(cl_uint) * 2 * num_threads_z;

    return shmem;
}

/*! \brief
 * Launch the pairlist prune only kernel for the given locality.
 * \p numParts tells in how many parts, i.e. calls the list will be pruned.
 */
void gpu_launch_kernel_pruneonly(NbnxmGpu* nb, const InteractionLocality iloc, const int numParts)
{
    NBAtomDataGpu*      adat         = nb->atdat;
    NBParamGpu*         nbp          = nb->nbparam;
    auto*               plist        = nb->plist[iloc].get();
    GpuTimers*          timers       = nb->timers;
    const DeviceStream& deviceStream = *nb->deviceStreams[iloc];
    bool                bDoTime      = nb->bDoTime;

    if (plist->haveFreshList)
    {
        GMX_ASSERT(numParts == 1, "With first pruning we expect 1 part");

        /* Set rollingPruningNumParts to signal that it is not set */
        plist->rollingPruningNumParts = 0;
        plist->rollingPruningPart     = 0;
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

    /* Use a local variable for part and update in plist, so we can return here
     * without duplicating the part increment code.
     */
    int part = plist->rollingPruningPart;

    plist->rollingPruningPart++;
    if (plist->rollingPruningPart >= plist->rollingPruningNumParts)
    {
        plist->rollingPruningPart = 0;
    }

    /* Compute the number of list entries to prune in this pass */
    int numSciInPart = (plist->numSci - part) / numParts;

    /* Don't launch the kernel if there is no work to do. */
    if (numSciInPart <= 0)
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
    /* kernel launch config */
    KernelLaunchConfig config;
    config.sharedMemorySize = calc_shmem_required_prune(num_threads_z);
    config.blockSize[0]     = c_clSize;
    config.blockSize[1]     = c_clSize;
    config.blockSize[2]     = num_threads_z;
    config.gridSize[0]      = numSciInPart;

    validate_global_work_size(config, 3, &nb->deviceContext_->deviceInfo());

    if (debug)
    {
        fprintf(debug,
                "Pruning GPU kernel launch configuration:\n\tLocal work size: %zux%zux%zu\n\t"
                "\tGlobal work size: %zux%zu\n\t#Super-clusters/clusters: %d/%d (%d)\n"
                "\tShMem: %zu\n",
                config.blockSize[0],
                config.blockSize[1],
                config.blockSize[2],
                config.blockSize[0] * config.gridSize[0],
                config.blockSize[1] * config.gridSize[1],
                plist->numSci * c_nbnxnGpuNumClusterPerSupercluster,
                c_nbnxnGpuNumClusterPerSupercluster,
                plist->numAtomsPerCluster,
                config.sharedMemorySize);
    }

    cl_nbparam_params_t nbparams_params;
    fillin_ocl_structures(nbp, &nbparams_params);

    auto*          timingEvent  = bDoTime ? timer->fetchNextEvent() : nullptr;
    constexpr char kernelName[] = "k_pruneonly";
    const auto     pruneKernel  = selectPruneKernel(nb->kernel_pruneonly, plist->haveFreshList);
    const auto     kernelArgs   = prepareGpuKernelArguments(pruneKernel,
                                                      config,
                                                      &nbparams_params,
                                                      &adat->xq,
                                                      &adat->shiftVec,
                                                      &plist->sci,
                                                      &plist->cjPacked,
                                                      &plist->imask,
                                                      &numParts,
                                                      &part);
    launchGpuKernel(pruneKernel, config, deviceStream, timingEvent, kernelName, kernelArgs);

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
}

} // namespace gmx
