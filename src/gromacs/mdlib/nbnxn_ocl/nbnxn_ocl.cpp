/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015,2016,2017,2018, by the GROMACS development team, led by
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
 *  \brief Define OpenCL implementation of nbnxn_gpu.h
 *
 *  \author Anca Hamuraru <anca@streamcomputing.eu>
 *  \author Teemu Virolainen <teemu@streamcomputing.eu>
 *  \author Dimitrios Karkoulis <dimitris.karkoulis@gmail.com>
 *  \author Szilárd Páll <pall.szilard@gmail.com>
 *  \ingroup module_mdlib
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

#include <assert.h>
#include <stdlib.h>

#if defined(_MSVC)
#include <limits>
#endif

#include "thread_mpi/atomic.h"

#include "gromacs/gpu_utils/gputraits_ocl.h"
#include "gromacs/gpu_utils/oclutils.h"
#include "gromacs/hardware/hw_info.h"
#include "gromacs/mdlib/force_flags.h"
#include "gromacs/mdlib/nb_verlet.h"
#include "gromacs/mdlib/nbnxn_consts.h"
#include "gromacs/mdlib/nbnxn_gpu.h"
#include "gromacs/mdlib/nbnxn_gpu_common.h"
#include "gromacs/mdlib/nbnxn_gpu_common_utils.h"
#include "gromacs/mdlib/nbnxn_gpu_data_mgmt.h"
#include "gromacs/mdlib/nbnxn_pairlist.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/timing/gpu_timing.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"

#include "nbnxn_ocl_internal.h"
#include "nbnxn_ocl_types.h"


/*! \brief Convenience constants */
//@{
static const int c_numClPerSupercl = c_nbnxnGpuNumClusterPerSupercluster;
static const int c_clSize          = c_nbnxnGpuClusterSize;
//@}


/*! \brief Validates the input global work size parameter.
 */
static inline void validate_global_work_size(const KernelLaunchConfig &config, int work_dim, const gmx_device_info_t *dinfo)
{
    cl_uint device_size_t_size_bits;
    cl_uint host_size_t_size_bits;

    assert(dinfo);

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
    host_size_t_size_bits   = static_cast<cl_uint>(sizeof(size_t) * 8);

    /* If sizeof(host size_t) <= sizeof(device size_t)
            => global_work_size components will always be valid
       else
            => get device limit for global work size and
            compare it against each component of global_work_size.
     */
    if (host_size_t_size_bits > device_size_t_size_bits)
    {
        size_t device_limit;

        device_limit = (1ull << device_size_t_size_bits) - 1;

        for (int i = 0; i < work_dim; i++)
        {
            if (global_work_size[i] > device_limit)
            {
                gmx_fatal(FARGS, "Watch out, the input system is too large to simulate!\n"
                          "The number of nonbonded work units (=number of super-clusters) exceeds the"
                          "device capabilities. Global work size limit exceeded (%zu > %zu)!",
                          global_work_size[i], device_limit);
            }
        }
    }
}

/* Constant arrays listing non-bonded kernel function names. The arrays are
 * organized in 2-dim arrays by: electrostatics and VDW type.
 *
 *  Note that the row- and column-order of function pointers has to match the
 *  order of corresponding enumerated electrostatics and vdw types, resp.,
 *  defined in nbnxn_cuda_types.h.
 */

/*! \brief Force-only kernel function names. */
static const char* nb_kfunc_noener_noprune_ptr[eelOclNR][evdwOclNR] =
{
    { "nbnxn_kernel_ElecCut_VdwLJ_F_opencl",            "nbnxn_kernel_ElecCut_VdwLJCombGeom_F_opencl",            "nbnxn_kernel_ElecCut_VdwLJCombLB_F_opencl",            "nbnxn_kernel_ElecCut_VdwLJFsw_F_opencl",            "nbnxn_kernel_ElecCut_VdwLJPsw_F_opencl",            "nbnxn_kernel_ElecCut_VdwLJEwCombGeom_F_opencl",            "nbnxn_kernel_ElecCut_VdwLJEwCombLB_F_opencl"            },
    { "nbnxn_kernel_ElecRF_VdwLJ_F_opencl",             "nbnxn_kernel_ElecRF_VdwLJCombGeom_F_opencl",             "nbnxn_kernel_ElecRF_VdwLJCombLB_F_opencl",             "nbnxn_kernel_ElecRF_VdwLJFsw_F_opencl",             "nbnxn_kernel_ElecRF_VdwLJPsw_F_opencl",             "nbnxn_kernel_ElecRF_VdwLJEwCombGeom_F_opencl",             "nbnxn_kernel_ElecRF_VdwLJEwCombLB_F_opencl"             },
    { "nbnxn_kernel_ElecEwQSTab_VdwLJ_F_opencl",        "nbnxn_kernel_ElecEwQSTab_VdwLJCombGeom_F_opencl",        "nbnxn_kernel_ElecEwQSTab_VdwLJCombLB_F_opencl",        "nbnxn_kernel_ElecEwQSTab_VdwLJFsw_F_opencl",        "nbnxn_kernel_ElecEwQSTab_VdwLJPsw_F_opencl",        "nbnxn_kernel_ElecEwQSTab_VdwLJEwCombGeom_F_opencl",        "nbnxn_kernel_ElecEwQSTab_VdwLJEwCombLB_F_opencl"        },
    { "nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJ_F_opencl", "nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJCombGeom_F_opencl", "nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJCombLB_F_opencl", "nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJFsw_F_opencl", "nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJPsw_F_opencl", "nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJEwCombGeom_F_opencl", "nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJEwCombLB_F_opencl" },
    { "nbnxn_kernel_ElecEw_VdwLJ_F_opencl",             "nbnxn_kernel_ElecEw_VdwLJCombGeom_F_opencl",             "nbnxn_kernel_ElecEw_VdwLJCombLB_F_opencl",             "nbnxn_kernel_ElecEw_VdwLJFsw_F_opencl",             "nbnxn_kernel_ElecEw_VdwLJPsw_F_opencl",             "nbnxn_kernel_ElecEw_VdwLJEwCombGeom_F_opencl",             "nbnxn_kernel_ElecEw_VdwLJEwCombLB_F_opencl"             },
    { "nbnxn_kernel_ElecEwTwinCut_VdwLJ_F_opencl",      "nbnxn_kernel_ElecEwTwinCut_VdwLJCombGeom_F_opencl",      "nbnxn_kernel_ElecEwTwinCut_VdwLJCombLB_F_opencl",      "nbnxn_kernel_ElecEwTwinCut_VdwLJFsw_F_opencl",      "nbnxn_kernel_ElecEwTwinCut_VdwLJPsw_F_opencl",      "nbnxn_kernel_ElecEwTwinCut_VdwLJEwCombGeom_F_opencl",      "nbnxn_kernel_ElecEwTwinCut_VdwLJEwCombLB_F_opencl"      }
};

/*! \brief Force + energy kernel function pointers. */
static const char* nb_kfunc_ener_noprune_ptr[eelOclNR][evdwOclNR] =
{
    { "nbnxn_kernel_ElecCut_VdwLJ_VF_opencl",            "nbnxn_kernel_ElecCut_VdwLJCombGeom_VF_opencl",            "nbnxn_kernel_ElecCut_VdwLJCombLB_VF_opencl",            "nbnxn_kernel_ElecCut_VdwLJFsw_VF_opencl",            "nbnxn_kernel_ElecCut_VdwLJPsw_VF_opencl",            "nbnxn_kernel_ElecCut_VdwLJEwCombGeom_VF_opencl",            "nbnxn_kernel_ElecCut_VdwLJEwCombLB_VF_opencl"            },
    { "nbnxn_kernel_ElecRF_VdwLJ_VF_opencl",             "nbnxn_kernel_ElecRF_VdwLJCombGeom_VF_opencl",             "nbnxn_kernel_ElecRF_VdwLJCombLB_VF_opencl",             "nbnxn_kernel_ElecRF_VdwLJFsw_VF_opencl",             "nbnxn_kernel_ElecRF_VdwLJPsw_VF_opencl",             "nbnxn_kernel_ElecRF_VdwLJEwCombGeom_VF_opencl",             "nbnxn_kernel_ElecRF_VdwLJEwCombLB_VF_opencl"             },
    { "nbnxn_kernel_ElecEwQSTab_VdwLJ_VF_opencl",        "nbnxn_kernel_ElecEwQSTab_VdwLJCombGeom_VF_opencl",        "nbnxn_kernel_ElecEwQSTab_VdwLJCombLB_VF_opencl",        "nbnxn_kernel_ElecEwQSTab_VdwLJFsw_VF_opencl",        "nbnxn_kernel_ElecEwQSTab_VdwLJPsw_VF_opencl",        "nbnxn_kernel_ElecEwQSTab_VdwLJEwCombGeom_VF_opencl",        "nbnxn_kernel_ElecEwQSTab_VdwLJEwCombLB_VF_opencl"        },
    { "nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJ_VF_opencl", "nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJCombGeom_VF_opencl", "nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJCombLB_VF_opencl", "nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJFsw_VF_opencl", "nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJPsw_VF_opencl", "nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJEwCombGeom_VF_opencl", "nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJEwCombLB_VF_opencl" },
    { "nbnxn_kernel_ElecEw_VdwLJ_VF_opencl",             "nbnxn_kernel_ElecEw_VdwLJCombGeom_VF_opencl",             "nbnxn_kernel_ElecEw_VdwLJCombLB_VF_opencl",             "nbnxn_kernel_ElecEw_VdwLJFsw_VF_opencl",             "nbnxn_kernel_ElecEw_VdwLJPsw_VF_opencl",             "nbnxn_kernel_ElecEw_VdwLJEwCombGeom_VF_opencl",             "nbnxn_kernel_ElecEw_VdwLJEwCombLB_VF_opencl"             },
    { "nbnxn_kernel_ElecEwTwinCut_VdwLJ_VF_opencl",      "nbnxn_kernel_ElecEwTwinCut_VdwLJCombGeom_VF_opencl",      "nbnxn_kernel_ElecEwTwinCut_VdwLJCombLB_VF_opencl",      "nbnxn_kernel_ElecEwTwinCut_VdwLJFsw_VF_opencl",      "nbnxn_kernel_ElecEwTwinCut_VdwLJPsw_VF_opencl",      "nbnxn_kernel_ElecEwTwinCut_VdwLJEwCombGeom_VF_opencl",      "nbnxn_kernel_ElecEwTwinCut_VdwLJEwCombLB_VF_opencl"      }
};

/*! \brief Force + pruning kernel function pointers. */
static const char* nb_kfunc_noener_prune_ptr[eelOclNR][evdwOclNR] =
{
    { "nbnxn_kernel_ElecCut_VdwLJ_F_prune_opencl",            "nbnxn_kernel_ElecCut_VdwLJCombGeom_F_prune_opencl",            "nbnxn_kernel_ElecCut_VdwLJCombLB_F_prune_opencl",            "nbnxn_kernel_ElecCut_VdwLJFsw_F_prune_opencl",            "nbnxn_kernel_ElecCut_VdwLJPsw_F_prune_opencl",            "nbnxn_kernel_ElecCut_VdwLJEwCombGeom_F_prune_opencl",            "nbnxn_kernel_ElecCut_VdwLJEwCombLB_F_prune_opencl"             },
    { "nbnxn_kernel_ElecRF_VdwLJ_F_prune_opencl",             "nbnxn_kernel_ElecRF_VdwLJCombGeom_F_prune_opencl",             "nbnxn_kernel_ElecRF_VdwLJCombLB_F_prune_opencl",             "nbnxn_kernel_ElecRF_VdwLJFsw_F_prune_opencl",             "nbnxn_kernel_ElecRF_VdwLJPsw_F_prune_opencl",             "nbnxn_kernel_ElecRF_VdwLJEwCombGeom_F_prune_opencl",             "nbnxn_kernel_ElecRF_VdwLJEwCombLB_F_prune_opencl"              },
    { "nbnxn_kernel_ElecEwQSTab_VdwLJ_F_prune_opencl",        "nbnxn_kernel_ElecEwQSTab_VdwLJCombGeom_F_prune_opencl",        "nbnxn_kernel_ElecEwQSTab_VdwLJCombLB_F_prune_opencl",        "nbnxn_kernel_ElecEwQSTab_VdwLJFsw_F_prune_opencl",        "nbnxn_kernel_ElecEwQSTab_VdwLJPsw_F_prune_opencl",        "nbnxn_kernel_ElecEwQSTab_VdwLJEwCombGeom_F_prune_opencl",        "nbnxn_kernel_ElecEwQSTab_VdwLJEwCombLB_F_prune_opencl"         },
    { "nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJ_F_prune_opencl", "nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJCombGeom_F_prune_opencl", "nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJCombLB_F_prune_opencl", "nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJFsw_F_prune_opencl", "nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJPsw_F_prune_opencl", "nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJEwCombGeom_F_prune_opencl", "nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJEwCombLB_F_prune_opencl"  },
    { "nbnxn_kernel_ElecEw_VdwLJ_F_prune_opencl",             "nbnxn_kernel_ElecEw_VdwLJCombGeom_F_prune_opencl",             "nbnxn_kernel_ElecEw_VdwLJCombLB_F_prune_opencl",             "nbnxn_kernel_ElecEw_VdwLJFsw_F_prune_opencl",             "nbnxn_kernel_ElecEw_VdwLJPsw_F_prune_opencl",             "nbnxn_kernel_ElecEw_VdwLJEwCombGeom_F_prune_opencl",             "nbnxn_kernel_ElecEw_VdwLJEwCombLB_F_prune_opencl"              },
    { "nbnxn_kernel_ElecEwTwinCut_VdwLJ_F_prune_opencl",      "nbnxn_kernel_ElecEwTwinCut_VdwLJCombGeom_F_prune_opencl",      "nbnxn_kernel_ElecEwTwinCut_VdwLJCombLB_F_prune_opencl",      "nbnxn_kernel_ElecEwTwinCut_VdwLJFsw_F_prune_opencl",      "nbnxn_kernel_ElecEwTwinCut_VdwLJPsw_F_prune_opencl",      "nbnxn_kernel_ElecEwTwinCut_VdwLJEwCombGeom_F_prune_opencl",      "nbnxn_kernel_ElecEwTwinCut_VdwLJEwCombLB_F_prune_opencl"       }
};

/*! \brief Force + energy + pruning kernel function pointers. */
static const char* nb_kfunc_ener_prune_ptr[eelOclNR][evdwOclNR] =
{
    { "nbnxn_kernel_ElecCut_VdwLJ_VF_prune_opencl",            "nbnxn_kernel_ElecCut_VdwLJCombGeom_VF_prune_opencl",            "nbnxn_kernel_ElecCut_VdwLJCombLB_VF_prune_opencl",            "nbnxn_kernel_ElecCut_VdwLJFsw_VF_prune_opencl",            "nbnxn_kernel_ElecCut_VdwLJPsw_VF_prune_opencl",            "nbnxn_kernel_ElecCut_VdwLJEwCombGeom_VF_prune_opencl",            "nbnxn_kernel_ElecCut_VdwLJEwCombLB_VF_prune_opencl"            },
    { "nbnxn_kernel_ElecRF_VdwLJ_VF_prune_opencl",             "nbnxn_kernel_ElecRF_VdwLJCombGeom_VF_prune_opencl",             "nbnxn_kernel_ElecRF_VdwLJCombLB_VF_prune_opencl",             "nbnxn_kernel_ElecRF_VdwLJFsw_VF_prune_opencl",             "nbnxn_kernel_ElecRF_VdwLJPsw_VF_prune_opencl",             "nbnxn_kernel_ElecRF_VdwLJEwCombGeom_VF_prune_opencl",             "nbnxn_kernel_ElecRF_VdwLJEwCombLB_VF_prune_opencl"             },
    { "nbnxn_kernel_ElecEwQSTab_VdwLJ_VF_prune_opencl",        "nbnxn_kernel_ElecEwQSTab_VdwLJCombGeom_VF_prune_opencl",        "nbnxn_kernel_ElecEwQSTab_VdwLJCombLB_VF_prune_opencl",        "nbnxn_kernel_ElecEwQSTab_VdwLJFsw_VF_prune_opencl",        "nbnxn_kernel_ElecEwQSTab_VdwLJPsw_VF_prune_opencl",        "nbnxn_kernel_ElecEwQSTab_VdwLJEwCombGeom_VF_prune_opencl",        "nbnxn_kernel_ElecEwQSTab_VdwLJEwCombLB_VF_prune_opencl"        },
    { "nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJ_VF_prune_opencl", "nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJCombGeom_VF_prune_opencl", "nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJCombLB_VF_prune_opencl", "nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJFsw_VF_prune_opencl", "nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJPsw_VF_prune_opencl", "nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJEwCombGeom_VF_prune_opencl", "nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJEwCombLB_VF_prune_opencl" },
    { "nbnxn_kernel_ElecEw_VdwLJ_VF_prune_opencl",             "nbnxn_kernel_ElecEw_VdwLJCombGeom_VF_prune_opencl",             "nbnxn_kernel_ElecEw_VdwLJCombLB_VF_prune_opencl",             "nbnxn_kernel_ElecEw_VdwLJFsw_VF_prune_opencl",             "nbnxn_kernel_ElecEw_VdwLJPsw_VF_prune_opencl",             "nbnxn_kernel_ElecEw_VdwLJEwCombGeom_VF_prune_opencl",             "nbnxn_kernel_ElecEw_VdwLJEwCombLB_VF_prune_opencl"             },
    { "nbnxn_kernel_ElecEwTwinCut_VdwLJ_VF_prune_opencl",      "nbnxn_kernel_ElecEwTwinCut_VdwLJCombGeom_VF_prune_opencl",      "nbnxn_kernel_ElecEwTwinCut_VdwLJCombLB_VF_prune_opencl",      "nbnxn_kernel_ElecEwTwinCut_VdwLJFsw_VF_prune_opencl",      "nbnxn_kernel_ElecEwTwinCut_VdwLJPsw_VF_prune_opencl",      "nbnxn_kernel_ElecEwTwinCut_VdwLJEwCombGeom_VF_prune_opencl",      "nbnxn_kernel_ElecEwTwinCut_VdwLJEwCombLB_VF_prune_opencl"      }
};

/*! \brief Return a pointer to the prune kernel version to be executed at the current invocation.
 *
 * \param[in] kernel_pruneonly  array of prune kernel objects
 * \param[in] firstPrunePass    true if the first pruning pass is being executed
 */
static inline cl_kernel selectPruneKernel(cl_kernel kernel_pruneonly[],
                                          bool      firstPrunePass)
{
    cl_kernel  *kernelPtr;

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
static inline cl_kernel select_nbnxn_kernel(gmx_nbnxn_ocl_t   *nb,
                                            int                eeltype,
                                            int                evdwtype,
                                            bool               bDoEne,
                                            bool               bDoPrune)
{
    const char* kernel_name_to_run;
    cl_kernel  *kernel_ptr;
    cl_int      cl_error;

    assert(eeltype  < eelOclNR);
    assert(evdwtype < evdwOclNR);

    if (bDoEne)
    {
        if (bDoPrune)
        {
            kernel_name_to_run = nb_kfunc_ener_prune_ptr[eeltype][evdwtype];
            kernel_ptr         = &(nb->kernel_ener_prune_ptr[eeltype][evdwtype]);
        }
        else
        {
            kernel_name_to_run = nb_kfunc_ener_noprune_ptr[eeltype][evdwtype];
            kernel_ptr         = &(nb->kernel_ener_noprune_ptr[eeltype][evdwtype]);
        }
    }
    else
    {
        if (bDoPrune)
        {
            kernel_name_to_run = nb_kfunc_noener_prune_ptr[eeltype][evdwtype];
            kernel_ptr         = &(nb->kernel_noener_prune_ptr[eeltype][evdwtype]);
        }
        else
        {
            kernel_name_to_run = nb_kfunc_noener_noprune_ptr[eeltype][evdwtype];
            kernel_ptr         = &(nb->kernel_noener_noprune_ptr[eeltype][evdwtype]);
        }
    }

    if (nullptr == kernel_ptr[0])
    {
        *kernel_ptr = clCreateKernel(nb->dev_rundata->program, kernel_name_to_run, &cl_error);
        assert(cl_error == CL_SUCCESS);
    }
    // TODO: handle errors

    return *kernel_ptr;
}

/*! \brief Calculates the amount of shared memory required by the nonbonded kernel in use.
 */
static inline int calc_shmem_required_nonbonded(int  vdwType,
                                                bool bPrefetchLjParam)
{
    int shmem;

    /* size of shmem (force-buffers/xq/atom type preloading) */
    /* NOTE: with the default kernel on sm3.0 we need shmem only for pre-loading */
    /* i-atom x+q in shared memory */
    shmem  = c_numClPerSupercl * c_clSize * sizeof(float) * 4; /* xqib */
    /* cj in shared memory, for both warps separately
     * TODO: in the "nowarp kernels we load cj only once  so the factor 2 is not needed.
     */
    shmem += 2 * c_nbnxnGpuJgroupSize * sizeof(int);           /* cjs  */
    if (bPrefetchLjParam)
    {
        if (useLjCombRule(vdwType))
        {
            /* i-atom LJ combination parameters in shared memory */
            shmem += c_numClPerSupercl * c_clSize * 2*sizeof(float); /* atib abused for ljcp, float2 */
        }
        else
        {
            /* i-atom types in shared memory */
            shmem += c_numClPerSupercl * c_clSize * sizeof(int); /* atib */
        }
    }
    /* force reduction buffers in shared memory */
    shmem += c_clSize * c_clSize * 3 * sizeof(float);    /* f_buf */
    /* Warp vote. In fact it must be * number of warps in block.. */
    shmem += sizeof(cl_uint) * 2;                        /* warp_any */
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
static void fillin_ocl_structures(cl_nbparam_t        *nbp,
                                  cl_nbparam_params_t *nbparams_params)
{
    nbparams_params->coulomb_tab_scale = nbp->coulomb_tab_scale;
    nbparams_params->c_rf              = nbp->c_rf;
    nbparams_params->dispersion_shift  = nbp->dispersion_shift;
    nbparams_params->eeltype           = nbp->eeltype;
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
    nbparams_params->vdwtype           = nbp->vdwtype;
    nbparams_params->vdw_switch        = nbp->vdw_switch;
}

/*! \brief Enqueues a wait for event completion.
 *
 * Then it releases the event and sets it to 0.
 * Don't use this function when more than one wait will be issued for the event.
 * Equivalent to Cuda Stream Sync. */
static void sync_ocl_event(cl_command_queue stream, cl_event *ocl_event)
{
    cl_int gmx_unused cl_error;

    /* Enqueue wait */
    cl_error = clEnqueueBarrierWithWaitList(stream, 1, ocl_event, nullptr);
    GMX_RELEASE_ASSERT(CL_SUCCESS == cl_error, ocl_get_error_string(cl_error).c_str());

    /* Release event and reset it to 0. It is ok to release it as enqueuewaitforevents performs implicit retain for events. */
    cl_error = clReleaseEvent(*ocl_event);
    assert(CL_SUCCESS == cl_error);
    *ocl_event = nullptr;
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
void nbnxn_gpu_launch_kernel(gmx_nbnxn_ocl_t               *nb,
                             const struct nbnxn_atomdata_t *nbatom,
                             int                            flags,
                             int                            iloc)
{
    int                  adat_begin, adat_len; /* local/nonlocal offset and length used for xq and f */
    /* OpenCL kernel launch-related stuff */
    cl_kernel            nb_kernel = nullptr;  /* fn pointer to the nonbonded kernel */

    cl_atomdata_t       *adat    = nb->atdat;
    cl_nbparam_t        *nbp     = nb->nbparam;
    cl_plist_t          *plist   = nb->plist[iloc];
    cl_timers_t         *t       = nb->timers;
    cl_command_queue     stream  = nb->stream[iloc];

    bool                 bCalcEner   = (flags & GMX_FORCE_ENERGY) != 0;
    int                  bCalcFshift = flags & GMX_FORCE_VIRIAL;
    bool                 bDoTime     = (nb->bDoTime) != 0;

    cl_nbparam_params_t  nbparams_params;

    /* Don't launch the non-local kernel if there is no work to do.
       Doing the same for the local kernel is more complicated, since the
       local part of the force array also depends on the non-local kernel.
       So to avoid complicating the code and to reduce the risk of bugs,
       we always call the local kernel, the local x+q copy and later (not in
       this function) the stream wait, local f copyback and the f buffer
       clearing. All these operations, except for the local interaction kernel,
       are needed for the non-local interactions. The skip of the local kernel
       call is taken care of later in this function. */
    if (canSkipWork(nb, iloc))
    {
        plist->haveFreshList = false;

        return;
    }

    /* calculate the atom data index range based on locality */
    if (LOCAL_I(iloc))
    {
        adat_begin  = 0;
        adat_len    = adat->natoms_local;
    }
    else
    {
        adat_begin  = adat->natoms_local;
        adat_len    = adat->natoms - adat->natoms_local;
    }

    /* beginning of timed HtoD section */
    if (bDoTime)
    {
        t->nb_h2d[iloc].openTimingRegion(stream);
    }

    /* HtoD x, q */
    ocl_copy_H2D_async(adat->xq, nbatom->x + adat_begin * 4, adat_begin*sizeof(float)*4,
                       adat_len * sizeof(float) * 4, stream, bDoTime ? t->nb_h2d[iloc].fetchNextEvent() : nullptr);

    if (bDoTime)
    {
        t->nb_h2d[iloc].closeTimingRegion(stream);
    }

    /* When we get here all misc operations issues in the local stream as well as
       the local xq H2D are done,
       so we record that in the local stream and wait for it in the nonlocal one. */
    if (nb->bUseTwoStreams)
    {
        if (iloc == eintLocal)
        {
            cl_int gmx_used_in_debug cl_error = clEnqueueMarkerWithWaitList(stream, 0, nullptr, &(nb->misc_ops_and_local_H2D_done));
            assert(CL_SUCCESS == cl_error);

            /* Based on the v1.2 section 5.13 of the OpenCL spec, a flush is needed
             * in the local stream in order to be able to sync with the above event
             * from the non-local stream.
             */
            cl_error = clFlush(stream);
            assert(CL_SUCCESS == cl_error);
        }
        else
        {
            sync_ocl_event(stream, &(nb->misc_ops_and_local_H2D_done));
        }
    }

    if (nbp->useDynamicPruning && plist->haveFreshList)
    {
        /* Prunes for rlistOuter and rlistInner, sets plist->haveFreshList=false
           (that's the way the timing accounting can distinguish between
           separate prune kernel and combined force+prune).
         */
        nbnxn_gpu_launch_kernel_pruneonly(nb, iloc, 1);
    }

    if (plist->nsci == 0)
    {
        /* Don't launch an empty local kernel (is not allowed with OpenCL).
         * TODO: Separate H2D and kernel launch into separate functions.
         */
        return;
    }

    /* beginning of timed nonbonded calculation section */
    if (bDoTime)
    {
        t->nb_k[iloc].openTimingRegion(stream);
    }

    /* get the pointer to the kernel flavor we need to use */
    nb_kernel = select_nbnxn_kernel(nb,
                                    nbp->eeltype,
                                    nbp->vdwtype,
                                    bCalcEner,
                                    (plist->haveFreshList && !nb->timers->didPrune[iloc]));

    /* kernel launch config */

    KernelLaunchConfig config;
    config.sharedMemorySize = calc_shmem_required_nonbonded(nbp->vdwtype, nb->bPrefetchLjParam);
    config.stream           = stream;
    config.blockSize[0]     = c_clSize;
    config.blockSize[1]     = c_clSize;
    config.gridSize[0]      = plist->nsci;

    validate_global_work_size(config, 3, nb->dev_info);

    if (debug)
    {
        fprintf(debug, "Non-bonded GPU launch configuration:\n\tLocal work size: %zux%zux%zu\n\t"
                "Global work size : %zux%zu\n\t#Super-clusters/clusters: %d/%d (%d)\n",
                config.blockSize[0], config.blockSize[1], config.blockSize[2],
                config.blockSize[0] * config.gridSize[0], config.blockSize[1] * config.gridSize[1], plist->nsci*c_numClPerSupercl,
                c_numClPerSupercl, plist->na_c);
    }

    fillin_ocl_structures(nbp, &nbparams_params);

    auto          *timingEvent  = bDoTime ? t->nb_k[iloc].fetchNextEvent() : nullptr;
    constexpr char kernelName[] = "k_calc_nb";
    if (useLjCombRule(nb->nbparam->vdwtype))
    {
        const auto kernelArgs = prepareGpuKernelArguments(nb_kernel, config,
                                                          &nbparams_params, &adat->xq, &adat->f, &adat->e_lj, &adat->e_el, &adat->fshift,
                                                          &adat->lj_comb,
                                                          &adat->shift_vec, &nbp->nbfp_climg2d, &nbp->nbfp_comb_climg2d, &nbp->coulomb_tab_climg2d,
                                                          &plist->sci, &plist->cj4, &plist->excl, &bCalcFshift);

        launchGpuKernel(nb_kernel, config, timingEvent, kernelName, kernelArgs);
    }
    else
    {
        const auto kernelArgs = prepareGpuKernelArguments(nb_kernel, config,
                                                          &adat->ntypes,
                                                          &nbparams_params, &adat->xq, &adat->f, &adat->e_lj, &adat->e_el, &adat->fshift,
                                                          &adat->atom_types,
                                                          &adat->shift_vec, &nbp->nbfp_climg2d, &nbp->nbfp_comb_climg2d, &nbp->coulomb_tab_climg2d,
                                                          &plist->sci, &plist->cj4, &plist->excl, &bCalcFshift);
        launchGpuKernel(nb_kernel, config, timingEvent, kernelName, kernelArgs);
    }

    if (bDoTime)
    {
        t->nb_k[iloc].closeTimingRegion(stream);
    }
}


/*! \brief Calculates the amount of shared memory required by the prune kernel.
 *
 *  Note that for the sake of simplicity we use the CUDA terminology "shared memory"
 *  for OpenCL local memory.
 *
 * \param[in] num_threads_z cj4 concurrency equal to the number of threads/work items in the 3-rd dimension.
 * \returns   the amount of local memory in bytes required by the pruning kernel
 */
static inline int calc_shmem_required_prune(const int num_threads_z)
{
    int shmem;

    /* i-atom x in shared memory (for convenience we load all 4 components including q) */
    shmem  = c_numClPerSupercl * c_clSize * sizeof(float)*4;
    /* cj in shared memory, for each warp separately
     * Note: only need to load once per wavefront, but to keep the code simple,
     * for now we load twice on AMD.
     */
    shmem += num_threads_z * c_nbnxnGpuClusterpairSplit * c_nbnxnGpuJgroupSize * sizeof(int);
    /* Warp vote, requires one uint per warp/32 threads per block. */
    shmem += sizeof(cl_uint) * 2*num_threads_z;

    return shmem;
}

void nbnxn_gpu_launch_kernel_pruneonly(gmx_nbnxn_gpu_t       *nb,
                                       int                    iloc,
                                       int                    numParts)
{
    cl_atomdata_t       *adat    = nb->atdat;
    cl_nbparam_t        *nbp     = nb->nbparam;
    cl_plist_t          *plist   = nb->plist[iloc];
    cl_timers_t         *t       = nb->timers;
    cl_command_queue     stream  = nb->stream[iloc];
    bool                 bDoTime = nb->bDoTime == CL_TRUE;

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
            GMX_ASSERT(numParts == plist->rollingPruningNumParts, "It is not allowed to change numParts in between list generation steps");
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
    int numSciInPart = (plist->nsci - part)/numParts;

    /* Don't launch the kernel if there is no work to do. */
    if (numSciInPart <= 0)
    {
        plist->haveFreshList = false;

        return;
    }

    GpuRegionTimer *timer = nullptr;
    if (bDoTime)
    {
        timer = &(plist->haveFreshList ? t->prune_k[iloc] : t->rollingPrune_k[iloc]);
    }

    /* beginning of timed prune calculation section */
    if (bDoTime)
    {
        timer->openTimingRegion(stream);
    }

    /* Kernel launch config:
     * - The thread block dimensions match the size of i-clusters, j-clusters,
     *   and j-cluster concurrency, in x, y, and z, respectively.
     * - The 1D block-grid contains as many blocks as super-clusters.
     */
    int       num_threads_z = getOclPruneKernelJ4Concurrency(nb->dev_info->vendor_e);
    cl_kernel pruneKernel   = selectPruneKernel(nb->kernel_pruneonly, plist->haveFreshList);

    /* kernel launch config */
    KernelLaunchConfig config;
    config.sharedMemorySize = calc_shmem_required_prune(num_threads_z);
    config.stream           = stream;
    config.blockSize[0]     = c_clSize;
    config.blockSize[1]     = c_clSize;
    config.blockSize[2]     = num_threads_z;
    config.gridSize[0]      = numSciInPart;

    validate_global_work_size(config, 3, nb->dev_info);

    if (debug)
    {
        fprintf(debug, "Pruning GPU kernel launch configuration:\n\tLocal work size: %zux%zux%zu\n\t"
                "\tGlobal work size: %zux%zu\n\t#Super-clusters/clusters: %d/%d (%d)\n"
                "\tShMem: %zu\n",
                config.blockSize[0], config.blockSize[1], config.blockSize[2],
                config.blockSize[0] * config.gridSize[0], config.blockSize[1] * config.gridSize[1], plist->nsci*c_numClPerSupercl,
                c_numClPerSupercl, plist->na_c, config.sharedMemorySize);
    }

    cl_nbparam_params_t  nbparams_params;
    fillin_ocl_structures(nbp, &nbparams_params);

    auto          *timingEvent  = bDoTime ? timer->fetchNextEvent() : nullptr;
    constexpr char kernelName[] = "k_pruneonly";
    const auto     kernelArgs   = prepareGpuKernelArguments(pruneKernel, config,
                                                            &nbparams_params, &adat->xq, &adat->shift_vec,
                                                            &plist->sci, &plist->cj4, &plist->imask, &numParts, &part);
    launchGpuKernel(pruneKernel, config, timingEvent, kernelName, kernelArgs);

    if (plist->haveFreshList)
    {
        plist->haveFreshList         = false;
        /* Mark that pruning has been done */
        nb->timers->didPrune[iloc] = true;
    }
    else
    {
        /* Mark that rolling pruning has been done */
        nb->timers->didRollingPrune[iloc] = true;
    }

    if (bDoTime)
    {
        timer->closeTimingRegion(stream);
    }
}

/*! \brief
 * Launch asynchronously the download of nonbonded forces from the GPU
 * (and energies/shift forces if required).
 */
void nbnxn_gpu_launch_cpyback(gmx_nbnxn_ocl_t               *nb,
                              const struct nbnxn_atomdata_t *nbatom,
                              int                            flags,
                              int                            aloc)
{
    cl_int gmx_unused cl_error;
    int               adat_begin, adat_len; /* local/nonlocal offset and length used for xq and f */

    /* determine interaction locality from atom locality */
    int              iloc = gpuAtomToInteractionLocality(aloc);

    cl_atomdata_t   *adat    = nb->atdat;
    cl_timers_t     *t       = nb->timers;
    bool             bDoTime = nb->bDoTime == CL_TRUE;
    cl_command_queue stream  = nb->stream[iloc];

    bool             bCalcEner   = (flags & GMX_FORCE_ENERGY) != 0;
    int              bCalcFshift = flags & GMX_FORCE_VIRIAL;


    /* don't launch non-local copy-back if there was no non-local work to do */
    if (canSkipWork(nb, iloc))
    {
        /* TODO An alternative way to signal that non-local work is
           complete is to use a clEnqueueMarker+clEnqueueBarrier
           pair. However, the use of bNonLocalStreamActive has the
           advantage of being local to the host, so probably minimizes
           overhead. Curiously, for NVIDIA OpenCL with an empty-domain
           test case, overall simulation performance was higher with
           the API calls, but this has not been tested on AMD OpenCL,
           so could be worth considering in future. */
        nb->bNonLocalStreamActive = CL_FALSE;
        return;
    }

    getGpuAtomRange(adat, aloc, &adat_begin, &adat_len);

    /* beginning of timed D2H section */
    if (bDoTime)
    {
        t->nb_d2h[iloc].openTimingRegion(stream);
    }

    /* With DD the local D2H transfer can only start after the non-local
       has been launched. */
    if (iloc == eintLocal && nb->bNonLocalStreamActive)
    {
        sync_ocl_event(stream, &(nb->nonlocal_done));
    }

    /* DtoH f */
    ocl_copy_D2H_async(nbatom->out[0].f + adat_begin * 3, adat->f, adat_begin*3*sizeof(float),
                       (adat_len)* adat->f_elem_size, stream, bDoTime ? t->nb_d2h[iloc].fetchNextEvent() : nullptr);

    /* kick off work */
    cl_error = clFlush(stream);
    assert(CL_SUCCESS == cl_error);

    /* After the non-local D2H is launched the nonlocal_done event can be
       recorded which signals that the local D2H can proceed. This event is not
       placed after the non-local kernel because we first need the non-local
       data back first. */
    if (iloc == eintNonlocal)
    {
        cl_error = clEnqueueMarkerWithWaitList(stream, 0, nullptr, &(nb->nonlocal_done));
        assert(CL_SUCCESS == cl_error);
        nb->bNonLocalStreamActive = CL_TRUE;
    }

    /* only transfer energies in the local stream */
    if (LOCAL_I(iloc))
    {
        /* DtoH fshift */
        if (bCalcFshift)
        {
            ocl_copy_D2H_async(nb->nbst.fshift, adat->fshift, 0,
                               SHIFTS * adat->fshift_elem_size, stream, bDoTime ? t->nb_d2h[iloc].fetchNextEvent() : nullptr);
        }

        /* DtoH energies */
        if (bCalcEner)
        {
            ocl_copy_D2H_async(nb->nbst.e_lj, adat->e_lj, 0,
                               sizeof(float), stream, bDoTime ? t->nb_d2h[iloc].fetchNextEvent() : nullptr);

            ocl_copy_D2H_async(nb->nbst.e_el, adat->e_el, 0,
                               sizeof(float), stream, bDoTime ? t->nb_d2h[iloc].fetchNextEvent() : nullptr);
        }
    }

    if (bDoTime)
    {
        t->nb_d2h[iloc].closeTimingRegion(stream);
    }
}


/*! \brief Selects the Ewald kernel type, analytical or tabulated, single or twin cut-off. */
int nbnxn_gpu_pick_ewald_kernel_type(bool bTwinCut)
{
    bool bUseAnalyticalEwald, bForceAnalyticalEwald, bForceTabulatedEwald;
    int  kernel_type;

    /* Benchmarking/development environment variables to force the use of
       analytical or tabulated Ewald kernel. */
    bForceAnalyticalEwald = (getenv("GMX_OCL_NB_ANA_EWALD") != nullptr);
    bForceTabulatedEwald  = (getenv("GMX_OCL_NB_TAB_EWALD") != nullptr);

    if (bForceAnalyticalEwald && bForceTabulatedEwald)
    {
        gmx_incons("Both analytical and tabulated Ewald OpenCL non-bonded kernels "
                   "requested through environment variables.");
    }

    /* OpenCL: By default, use analytical Ewald
     * TODO: tabulated does not work, it needs fixing, see init_nbparam() in nbnxn_ocl_data_mgmt.cpp
     *
     * TODO: decide if dev_info parameter should be added to recognize NVIDIA CC>=3.0 devices.
     *
     */
    /* By default use analytical Ewald. */
    bUseAnalyticalEwald = true;
    if (bForceAnalyticalEwald)
    {
        if (debug)
        {
            fprintf(debug, "Using analytical Ewald OpenCL kernels\n");
        }
    }
    else if (bForceTabulatedEwald)
    {
        bUseAnalyticalEwald = false;

        if (debug)
        {
            fprintf(debug, "Using tabulated Ewald OpenCL kernels\n");
        }
    }

    /* Use twin cut-off kernels if requested by bTwinCut or the env. var.
       forces it (use it for debugging/benchmarking only). */
    if (!bTwinCut && (getenv("GMX_OCL_NB_EWALD_TWINCUT") == nullptr))
    {
        kernel_type = bUseAnalyticalEwald ? eelOclEWALD_ANA : eelOclEWALD_TAB;
    }
    else
    {
        kernel_type = bUseAnalyticalEwald ? eelOclEWALD_ANA_TWIN : eelOclEWALD_TAB_TWIN;
    }

    return kernel_type;
}
