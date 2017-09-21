/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015,2016,2017, by the GROMACS development team, led by
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

#include "gromacs/gpu_utils/oclutils.h"
#include "gromacs/hardware/hw_info.h"
#include "gromacs/mdlib/force_flags.h"
#include "gromacs/mdlib/nb_verlet.h"
#include "gromacs/mdlib/nbnxn_consts.h"
#include "gromacs/mdlib/nbnxn_gpu.h"
#include "gromacs/mdlib/nbnxn_gpu_common.h"
#include "gromacs/mdlib/nbnxn_gpu_data_mgmt.h"
#include "gromacs/mdlib/nbnxn_pairlist.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/timing/gpu_timing.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"

#include "nbnxn_ocl_internal.h"
#include "nbnxn_ocl_types.h"

#if defined TEXOBJ_SUPPORTED && __CUDA_ARCH__ >= 300
#define USE_TEXOBJ
#endif

/*! \brief Convenience constants */
//@{
static const int c_numClPerSupercl = c_nbnxnGpuNumClusterPerSupercluster;
static const int c_clSize          = c_nbnxnGpuClusterSize;
//@}

/*! \brief Always/never run the energy/pruning kernels -- only for benchmarking purposes */
//@{
static bool always_ener  = (getenv("GMX_GPU_ALWAYS_ENER") != NULL);
static bool never_ener   = (getenv("GMX_GPU_NEVER_ENER") != NULL);
static bool always_prune = (getenv("GMX_GPU_ALWAYS_PRUNE") != NULL);
//@}

/* Uncomment this define to enable kernel debugging */
//#define DEBUG_OCL

/*! \brief Specifies which kernel run to debug */
#define DEBUG_RUN_STEP 2

/*! \brief Validates the input global work size parameter.
 */
static inline void validate_global_work_size(size_t *global_work_size, int work_dim, const gmx_device_info_t *dinfo)
{
    cl_uint device_size_t_size_bits;
    cl_uint host_size_t_size_bits;

    assert(dinfo);

    /* Each component of a global_work_size must not exceed the range given by the
       sizeof(device size_t) for the device on which the kernel execution will
       be enqueued. See:
       https://www.khronos.org/registry/cl/sdk/1.0/docs/man/xhtml/clEnqueueNDRangeKernel.html
     */
    device_size_t_size_bits = dinfo->adress_bits;
    host_size_t_size_bits   = (cl_uint)(sizeof(size_t) * 8);

    /* If sizeof(host size_t) <= sizeof(device size_t)
            => global_work_size components will always be valid
       else
            => get device limit for global work size and
            compare it against each component of global_work_size.
     */
    if (host_size_t_size_bits > device_size_t_size_bits)
    {
        size_t device_limit;

        device_limit = (((size_t)1) << device_size_t_size_bits) - 1;

        for (int i = 0; i < work_dim; i++)
        {
            if (global_work_size[i] > device_limit)
            {
                gmx_fatal(FARGS, "Watch out, the input system is too large to simulate!\n"
                          "The number of nonbonded work units (=number of super-clusters) exceeds the"
                          "device capabilities. Global work size limit exceeded (%d > %d)!",
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

    if (NULL == kernel_ptr[0])
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
    /* cj in shared memory, for both warps separately */
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
#ifdef CL_VERSION_1_2
    cl_error = clEnqueueBarrierWithWaitList(stream, 1, ocl_event, NULL);
#else
    cl_error = clEnqueueWaitForEvents(stream, 1, ocl_event);
#endif

    GMX_RELEASE_ASSERT(CL_SUCCESS == cl_error, ocl_get_error_string(cl_error).c_str());

    /* Release event and reset it to 0. It is ok to release it as enqueuewaitforevents performs implicit retain for events. */
    cl_error = clReleaseEvent(*ocl_event);
    assert(CL_SUCCESS == cl_error);
    *ocl_event = 0;
}

/*! \brief Returns the duration in milliseconds for the command associated with the event.
 *
 * It then releases the event and sets it to 0.
 * Before calling this function, make sure the command has finished either by
 * calling clFinish or clWaitForEvents.
 * The function returns 0.0 if the input event, *ocl_event, is 0.
 * Don't use this function when more than one wait will be issued for the event.
 */
static double ocl_event_elapsed_ms(cl_event *ocl_event)
{
    cl_int gmx_unused cl_error;
    cl_ulong          start_ns, end_ns;
    double            elapsed_ms;

    elapsed_ms = 0.0;
    assert(NULL != ocl_event);

    if (*ocl_event)
    {
        cl_error = clGetEventProfilingInfo(*ocl_event, CL_PROFILING_COMMAND_START,
                                           sizeof(cl_ulong), &start_ns, NULL);
        assert(CL_SUCCESS == cl_error);

        cl_error = clGetEventProfilingInfo(*ocl_event, CL_PROFILING_COMMAND_END,
                                           sizeof(cl_ulong), &end_ns, NULL);
        assert(CL_SUCCESS == cl_error);

        clReleaseEvent(*ocl_event);
        *ocl_event = 0;

        elapsed_ms = (end_ns - start_ns) / 1000000.0;
    }

    return elapsed_ms;
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
    cl_int               cl_error;
    int                  adat_begin, adat_len; /* local/nonlocal offset and length used for xq and f */
    /* OpenCL kernel launch-related stuff */
    int                  shmem;
    size_t               local_work_size[3], global_work_size[3];
    cl_kernel            nb_kernel = NULL; /* fn pointer to the nonbonded kernel */

    cl_atomdata_t       *adat    = nb->atdat;
    cl_nbparam_t        *nbp     = nb->nbparam;
    cl_plist_t          *plist   = nb->plist[iloc];
    cl_timers_t         *t       = nb->timers;
    cl_command_queue     stream  = nb->stream[iloc];

    bool                 bCalcEner   = flags & GMX_FORCE_ENERGY;
    int                  bCalcFshift = flags & GMX_FORCE_VIRIAL;
    bool                 bDoTime     = nb->bDoTime;
    cl_uint              arg_no;

    cl_nbparam_params_t  nbparams_params;
#ifdef DEBUG_OCL
    float              * debug_buffer_h;
    size_t               debug_buffer_size;
#endif

    /* turn energy calculation always on/off (for debugging/testing only) */
    bCalcEner = (bCalcEner || always_ener) && !never_ener;

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

    /* HtoD x, q */
    ocl_copy_H2D_async(adat->xq, nbatom->x + adat_begin * 4, adat_begin*sizeof(float)*4,
                       adat_len * sizeof(float) * 4, stream, bDoTime ? (&(t->nb_h2d[iloc])) : NULL);

    /* When we get here all misc operations issues in the local stream as well as
       the local xq H2D are done,
       so we record that in the local stream and wait for it in the nonlocal one. */
    if (nb->bUseTwoStreams)
    {
        if (iloc == eintLocal)
        {
#ifdef CL_VERSION_1_2
            cl_error = clEnqueueMarkerWithWaitList(stream, 0, NULL, &(nb->misc_ops_and_local_H2D_done));
#else
            cl_error = clEnqueueMarker(stream, &(nb->misc_ops_and_local_H2D_done));
#endif
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

    /* get the pointer to the kernel flavor we need to use */
    nb_kernel = select_nbnxn_kernel(nb,
                                    nbp->eeltype,
                                    nbp->vdwtype,
                                    bCalcEner,
                                    (plist->haveFreshList && !nb->timers->didPrune[iloc]) || always_prune);

    /* kernel launch config */
    local_work_size[0] = c_clSize;
    local_work_size[1] = c_clSize;
    local_work_size[2] = 1;

    global_work_size[0] = plist->nsci * local_work_size[0];
    global_work_size[1] = 1 * local_work_size[1];
    global_work_size[2] = 1 * local_work_size[2];

    validate_global_work_size(global_work_size, 3, nb->dev_info);

    shmem     = calc_shmem_required_nonbonded(nbp->vdwtype, nb->bPrefetchLjParam);

#ifdef DEBUG_OCL
    {
        static int run_step = 1;

        if (DEBUG_RUN_STEP == run_step)
        {
            debug_buffer_size = global_work_size[0] * global_work_size[1] * global_work_size[2] * sizeof(float);
            debug_buffer_h    = (float*)calloc(1, debug_buffer_size);
            assert(NULL != debug_buffer_h);

            if (NULL == nb->debug_buffer)
            {
                nb->debug_buffer = clCreateBuffer(nb->dev_rundata->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                                                  debug_buffer_size, debug_buffer_h, &cl_error);

                assert(CL_SUCCESS == cl_error);
            }
        }

        run_step++;
    }
#endif
    if (debug)
    {
        fprintf(debug, "Non-bonded GPU launch configuration:\n\tLocal work size: %dx%dx%d\n\t"
                "Global work size : %dx%d\n\t#Super-clusters/clusters: %d/%d (%d)\n",
                (int)(local_work_size[0]), (int)(local_work_size[1]), (int)(local_work_size[2]),
                (int)(global_work_size[0]), (int)(global_work_size[1]), plist->nsci*c_numClPerSupercl,
                c_numClPerSupercl, plist->na_c);
    }

    fillin_ocl_structures(nbp, &nbparams_params);

    arg_no    = 0;
    cl_error  = CL_SUCCESS;
    if (!useLjCombRule(nb->nbparam->vdwtype))
    {
        cl_error  = clSetKernelArg(nb_kernel, arg_no++, sizeof(int), &(adat->ntypes));
    }
    cl_error |= clSetKernelArg(nb_kernel, arg_no++, sizeof(nbparams_params), &(nbparams_params));
    cl_error |= clSetKernelArg(nb_kernel, arg_no++, sizeof(cl_mem), &(adat->xq));
    cl_error |= clSetKernelArg(nb_kernel, arg_no++, sizeof(cl_mem), &(adat->f));
    cl_error |= clSetKernelArg(nb_kernel, arg_no++, sizeof(cl_mem), &(adat->e_lj));
    cl_error |= clSetKernelArg(nb_kernel, arg_no++, sizeof(cl_mem), &(adat->e_el));
    cl_error |= clSetKernelArg(nb_kernel, arg_no++, sizeof(cl_mem), &(adat->fshift));
    if (useLjCombRule(nb->nbparam->vdwtype))
    {
        cl_error |= clSetKernelArg(nb_kernel, arg_no++, sizeof(cl_mem), &(adat->lj_comb));
    }
    else
    {
        cl_error |= clSetKernelArg(nb_kernel, arg_no++, sizeof(cl_mem), &(adat->atom_types));
    }
    cl_error |= clSetKernelArg(nb_kernel, arg_no++, sizeof(cl_mem), &(adat->shift_vec));
    cl_error |= clSetKernelArg(nb_kernel, arg_no++, sizeof(cl_mem), &(nbp->nbfp_climg2d));
    cl_error |= clSetKernelArg(nb_kernel, arg_no++, sizeof(cl_mem), &(nbp->nbfp_comb_climg2d));
    cl_error |= clSetKernelArg(nb_kernel, arg_no++, sizeof(cl_mem), &(nbp->coulomb_tab_climg2d));
    cl_error |= clSetKernelArg(nb_kernel, arg_no++, sizeof(cl_mem), &(plist->sci));
    cl_error |= clSetKernelArg(nb_kernel, arg_no++, sizeof(cl_mem), &(plist->cj4));
    cl_error |= clSetKernelArg(nb_kernel, arg_no++, sizeof(cl_mem), &(plist->excl));
    cl_error |= clSetKernelArg(nb_kernel, arg_no++, sizeof(int), &bCalcFshift);
    cl_error |= clSetKernelArg(nb_kernel, arg_no++, shmem, NULL);
    cl_error |= clSetKernelArg(nb_kernel, arg_no++, sizeof(cl_mem), &(nb->debug_buffer));

    assert(cl_error == CL_SUCCESS);

    if (cl_error)
    {
        printf("OpenCL error: %s\n", ocl_get_error_string(cl_error).c_str());
    }
    cl_error = clEnqueueNDRangeKernel(stream, nb_kernel, 3, NULL, global_work_size, local_work_size, 0, NULL, bDoTime ? &(t->nb_k[iloc]) : NULL);
    assert(cl_error == CL_SUCCESS);

#ifdef DEBUG_OCL
    {
        static int run_step = 1;

        if (DEBUG_RUN_STEP == run_step)
        {
            FILE *pf;
            char  file_name[256] = {0};

            ocl_copy_D2H_async(debug_buffer_h, nb->debug_buffer, 0,
                               debug_buffer_size, stream, NULL);

            // Make sure all data has been transfered back from device
            clFinish(stream);

            printf("\nWriting debug_buffer to debug_buffer_ocl.txt...");

            sprintf(file_name, "debug_buffer_ocl_%d.txt", DEBUG_RUN_STEP);
            pf = fopen(file_name, "wt");
            assert(pf != NULL);

            fprintf(pf, "%20s", "");
            for (int j = 0; j < global_work_size[0]; j++)
            {
                char label[20];
                sprintf(label, "(wIdx=%2d thIdx=%2d)", j / local_work_size[0], j % local_work_size[0]);
                fprintf(pf, "%20s", label);
            }

            for (int i = 0; i < global_work_size[1]; i++)
            {
                char label[20];
                sprintf(label, "(wIdy=%2d thIdy=%2d)", i / local_work_size[1], i % local_work_size[1]);
                fprintf(pf, "\n%20s", label);

                for (int j = 0; j < global_work_size[0]; j++)
                {
                    fprintf(pf, "%20.5f", debug_buffer_h[i * global_work_size[0] + j]);
                }

                //fprintf(pf, "\n");
            }

            fclose(pf);

            printf(" done.\n");


            free(debug_buffer_h);
            debug_buffer_h = NULL;
        }

        run_step++;
    }
#endif
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
    /* cj in shared memory, for each warp separately */
    shmem += num_threads_z * c_nbnxnGpuClusterpairSplit * c_nbnxnGpuJgroupSize * sizeof(int);
    /* Warp vote, requires one uint per warp/32 threads per block. */
    shmem += sizeof(cl_uint) * 2*num_threads_z;

    return shmem;
}

void nbnxn_gpu_launch_kernel_pruneonly(gmx_nbnxn_gpu_t       *nb,
                                       int                    iloc,
                                       int                    numParts)
{
    cl_int               cl_error;

    cl_atomdata_t       *adat    = nb->atdat;
    cl_nbparam_t        *nbp     = nb->nbparam;
    cl_plist_t          *plist   = nb->plist[iloc];
    cl_timers_t         *t       = nb->timers;
    cl_command_queue     stream  = nb->stream[iloc];
    bool                 bDoTime = nb->bDoTime;

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

    /* Kernel launch config:
     * - The thread block dimensions match the size of i-clusters, j-clusters,
     *   and j-cluster concurrency, in x, y, and z, respectively.
     * - The 1D block-grid contains as many blocks as super-clusters.
     */
    int       num_threads_z = getOclPruneKernelJ4Concurrency(nb->dev_info->vendor_e);
    cl_kernel pruneKernel   = selectPruneKernel(nb->kernel_pruneonly, plist->haveFreshList);

    /* kernel launch config */
    size_t  local_work_size[3], global_work_size[3];
    local_work_size[0] = c_clSize;
    local_work_size[1] = c_clSize;
    local_work_size[2] = num_threads_z;

    global_work_size[0] = numSciInPart * local_work_size[0];
    global_work_size[1] = 1 * local_work_size[1];
    global_work_size[2] = 1 * local_work_size[2];

    validate_global_work_size(global_work_size, 3, nb->dev_info);

    int shmem = calc_shmem_required_prune(num_threads_z);

    if (debug)
    {
        fprintf(debug, "Pruning GPU kernel launch configuration:\n\tLocal work size: %dx%dx%d\n\t"
                "\tGlobal work size: %dx%d\n\t#Super-clusters/clusters: %d/%d (%d)\n"
                "\tShMem: %d\n",
                (int)(local_work_size[0]), (int)(local_work_size[1]), (int)(local_work_size[2]),
                (int)(global_work_size[0]), (int)(global_work_size[1]), plist->nsci*c_numClPerSupercl,
                c_numClPerSupercl, plist->na_c, shmem);
    }

    cl_nbparam_params_t  nbparams_params;
    fillin_ocl_structures(nbp, &nbparams_params);

    cl_uint  arg_no = 0;
    cl_error = CL_SUCCESS;

    cl_error |= clSetKernelArg(pruneKernel, arg_no++, sizeof(nbparams_params), &(nbparams_params));
    cl_error |= clSetKernelArg(pruneKernel, arg_no++, sizeof(cl_mem), &(adat->xq));
    cl_error |= clSetKernelArg(pruneKernel, arg_no++, sizeof(cl_mem), &(adat->shift_vec));
    cl_error |= clSetKernelArg(pruneKernel, arg_no++, sizeof(cl_mem), &(plist->sci));
    cl_error |= clSetKernelArg(pruneKernel, arg_no++, sizeof(cl_mem), &(plist->cj4));
    cl_error |= clSetKernelArg(pruneKernel, arg_no++, sizeof(cl_mem), &(plist->imask));
    cl_error |= clSetKernelArg(pruneKernel, arg_no++, sizeof(int), &(numParts));
    cl_error |= clSetKernelArg(pruneKernel, arg_no++, sizeof(int), &(part));
    cl_error |= clSetKernelArg(pruneKernel, arg_no++, shmem, nullptr);
    assert(cl_error == CL_SUCCESS);

    cl_event *pruneEventPtr = nullptr;
    if (bDoTime)
    {
        pruneEventPtr = plist->haveFreshList ? &t->prune_k[iloc] : &t->rollingPrune_k[iloc];
    }

    cl_error = clEnqueueNDRangeKernel(stream, pruneKernel, 3,
                                      nullptr, global_work_size, local_work_size,
                                      0, nullptr, pruneEventPtr);
    GMX_RELEASE_ASSERT(CL_SUCCESS == cl_error, ocl_get_error_string(cl_error).c_str());

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
    int               iloc = -1;

    /* determine interaction locality from atom locality */
    if (LOCAL_A(aloc))
    {
        iloc = eintLocal;
    }
    else if (NONLOCAL_A(aloc))
    {
        iloc = eintNonlocal;
    }
    else
    {
        char stmp[STRLEN];
        sprintf(stmp, "Invalid atom locality passed (%d); valid here is only "
                "local (%d) or nonlocal (%d)", aloc, eatLocal, eatNonlocal);

        gmx_incons(stmp);
    }

    cl_atomdata_t   *adat    = nb->atdat;
    cl_timers_t     *t       = nb->timers;
    bool             bDoTime = nb->bDoTime;
    cl_command_queue stream  = nb->stream[iloc];

    bool             bCalcEner   = flags & GMX_FORCE_ENERGY;
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
        nb->bNonLocalStreamActive = false;
        return;
    }

    /* calculate the atom data index range based on locality */
    if (LOCAL_A(aloc))
    {
        adat_begin  = 0;
        adat_len    = adat->natoms_local;
    }
    else
    {
        adat_begin  = adat->natoms_local;
        adat_len    = adat->natoms - adat->natoms_local;
    }

    /* beginning of timed D2H section */

    /* With DD the local D2H transfer can only start after the non-local
       has been launched. */
    if (iloc == eintLocal && nb->bNonLocalStreamActive)
    {
        sync_ocl_event(stream, &(nb->nonlocal_done));
    }

    /* DtoH f */
    ocl_copy_D2H_async(nbatom->out[0].f + adat_begin * 3, adat->f, adat_begin*3*sizeof(float),
                       (adat_len)* adat->f_elem_size, stream, bDoTime ? &(t->nb_d2h_f[iloc]) : NULL);

    /* kick off work */
    cl_error = clFlush(stream);
    assert(CL_SUCCESS == cl_error);

    /* After the non-local D2H is launched the nonlocal_done event can be
       recorded which signals that the local D2H can proceed. This event is not
       placed after the non-local kernel because we first need the non-local
       data back first. */
    if (iloc == eintNonlocal)
    {
#ifdef CL_VERSION_1_2
        cl_error = clEnqueueMarkerWithWaitList(stream, 0, NULL, &(nb->nonlocal_done));
#else
        cl_error = clEnqueueMarker(stream, &(nb->nonlocal_done));
#endif
        assert(CL_SUCCESS == cl_error);
        nb->bNonLocalStreamActive = true;
    }

    /* only transfer energies in the local stream */
    if (LOCAL_I(iloc))
    {
        /* DtoH fshift */
        if (bCalcFshift)
        {
            ocl_copy_D2H_async(nb->nbst.fshift, adat->fshift, 0,
                               SHIFTS * adat->fshift_elem_size, stream, bDoTime ? &(t->nb_d2h_fshift[iloc]) : NULL);
        }

        /* DtoH energies */
        if (bCalcEner)
        {
            ocl_copy_D2H_async(nb->nbst.e_lj, adat->e_lj, 0,
                               sizeof(float), stream, bDoTime ? &(t->nb_d2h_e_lj[iloc]) : NULL);

            ocl_copy_D2H_async(nb->nbst.e_el, adat->e_el, 0,
                               sizeof(float), stream, bDoTime ? &(t->nb_d2h_e_el[iloc]) : NULL);
        }
    }
}

/*! \brief Count pruning kernel time if either kernel has been triggered
 *
 *  We do the accounting for either of the two pruning kernel flavors:
 *   - 1st pass prune: ran during the current step (prior to the force kernel);
 *   - rolling prune:  ran at the end of the previous step (prior to the current step H2D xq);
 *
 * Note that the resetting of cu_timers_t::didPrune and cu_timers_t::didRollingPrune should happen
 * after calling this function.
 *
 * \param[inout] timers   structs with OCL timer objects
 * \param[inout] timings  GPU task timing data
 * \param[in] iloc        interaction locality
 */
static void countPruneKernelTime(cl_timers_t         *timers,
                                 gmx_wallclock_gpu_t *timings,
                                 const int            iloc)
{
    // We might have not done any pruning (e.g. if we skipped with empty domains).
    if (!timers->didPrune[iloc] && !timers->didRollingPrune[iloc])
    {
        return;
    }

    if (timers->didPrune[iloc])
    {
        timings->pruneTime.c++;
        timings->pruneTime.t += ocl_event_elapsed_ms(&timers->prune_k[iloc]);
    }

    if (timers->didRollingPrune[iloc])
    {
        timings->dynamicPruneTime.c++;
        timings->dynamicPruneTime.t += ocl_event_elapsed_ms(&timers->rollingPrune_k[iloc]);
    }
}

/*! \brief
 * Wait for the asynchronously launched nonbonded calculations and data
 * transfers to finish.
 */
void nbnxn_gpu_wait_for_gpu(gmx_nbnxn_ocl_t *nb,
                            int flags, int aloc,
                            real *e_lj, real *e_el, rvec *fshift)
{
    /* NOTE:  only implemented for single-precision at this time */
    cl_int gmx_unused      cl_error;
    int                    iloc = -1;

    /* determine interaction locality from atom locality */
    if (LOCAL_A(aloc))
    {
        iloc = eintLocal;
    }
    else if (NONLOCAL_A(aloc))
    {
        iloc = eintNonlocal;
    }
    else
    {
        char stmp[STRLEN];
        sprintf(stmp, "Invalid atom locality passed (%d); valid here is only "
                "local (%d) or nonlocal (%d)", aloc, eatLocal, eatNonlocal);
        gmx_incons(stmp);
    }

    cl_plist_t                 *plist    = nb->plist[iloc];
    cl_timers_t                *timers   = nb->timers;
    struct gmx_wallclock_gpu_t *timings  = nb->timings;
    cl_nb_staging               nbst     = nb->nbst;

    bool                        bCalcEner   = flags & GMX_FORCE_ENERGY;
    int                         bCalcFshift = flags & GMX_FORCE_VIRIAL;

    /* turn energy calculation always on/off (for debugging/testing only) */
    bCalcEner = (bCalcEner || always_ener) && !never_ener;

    /* Launch wait/update timers & counters and do reduction into staging buffers
       BUT skip it when during the non-local phase there was actually no work to do.
       This is consistent with nbnxn_gpu_launch_kernel.

       NOTE: if timing with multiple GPUs (streams) becomes possible, the
       counters could end up being inconsistent due to not being incremented
       on some of the nodes! */
    if (!canSkipWork(nb, iloc))
    {
        /* Actual sync point. Waits for everything to be finished in the command queue. TODO: Find out if a more fine grained solution is needed */
        cl_error = clFinish(nb->stream[iloc]);
        assert(CL_SUCCESS == cl_error);

        /* timing data accumulation */
        if (nb->bDoTime)
        {
            /* only increase counter once (at local F wait) */
            if (LOCAL_I(iloc))
            {
                timings->nb_c++;
                timings->ktime[plist->haveFreshList ? 1 : 0][bCalcEner ? 1 : 0].c += 1;
            }

            /* kernel timings */

            timings->ktime[plist->haveFreshList ? 1 : 0][bCalcEner ? 1 : 0].t +=
                ocl_event_elapsed_ms(timers->nb_k + iloc);

            /* X/q H2D and F D2H timings */
            timings->nb_h2d_t += ocl_event_elapsed_ms(timers->nb_h2d        + iloc);
            timings->nb_d2h_t += ocl_event_elapsed_ms(timers->nb_d2h_f      + iloc);
            timings->nb_d2h_t += ocl_event_elapsed_ms(timers->nb_d2h_fshift + iloc);
            timings->nb_d2h_t += ocl_event_elapsed_ms(timers->nb_d2h_e_el   + iloc);
            timings->nb_d2h_t += ocl_event_elapsed_ms(timers->nb_d2h_e_lj   + iloc);

            /* Count the pruning kernel times for both cases:1st pass (at search step)
               and rolling pruning (if called at the previous step).
               We do the accounting here as this is the only sync point where we
               know (without checking or additional sync-ing) that prune tasks in
               in the current stream have completed (having just blocking-waited
               for the force D2H). */
            countPruneKernelTime(timers, timings, iloc);

            /* only count atdat and pair-list H2D at pair-search step */
            if (timers->didPairlistH2D[iloc])
            {
                /* atdat transfer timing (add only once, at local F wait) */
                if (LOCAL_A(aloc))
                {
                    timings->pl_h2d_c++;
                    timings->pl_h2d_t += ocl_event_elapsed_ms(&(timers->atdat));
                }

                timings->pl_h2d_t +=
                    ocl_event_elapsed_ms(timers->pl_h2d_sci     + iloc) +
                    ocl_event_elapsed_ms(timers->pl_h2d_cj4     + iloc) +
                    ocl_event_elapsed_ms(timers->pl_h2d_excl    + iloc);

                /* Clear the timing flag for the next step */
                timers->didPairlistH2D[iloc] = false;
            }
        }

        /* add up energies and shift forces (only once at local F wait) */
        if (LOCAL_I(iloc))
        {
            if (bCalcEner)
            {
                *e_lj += *nbst.e_lj;
                *e_el += *nbst.e_el;
            }

            if (bCalcFshift)
            {
                for (int i = 0; i < SHIFTS; i++)
                {
                    fshift[i][0] += (nbst.fshift)[i][0];
                    fshift[i][1] += (nbst.fshift)[i][1];
                    fshift[i][2] += (nbst.fshift)[i][2];
                }
            }
        }
    }

    /* Always reset both pruning flags (doesn't hurt doing it even when timing is off). */
    timers->didPrune[iloc] = timers->didRollingPrune[iloc] = false;

    /* Turn off initial list pruning (doesn't hurt if this is not pair-search step). */
    plist->haveFreshList = false;
}

/*! \brief Selects the Ewald kernel type, analytical or tabulated, single or twin cut-off. */
int nbnxn_gpu_pick_ewald_kernel_type(bool bTwinCut)
{
    bool bUseAnalyticalEwald, bForceAnalyticalEwald, bForceTabulatedEwald;
    int  kernel_type;

    /* Benchmarking/development environment variables to force the use of
       analytical or tabulated Ewald kernel. */
    bForceAnalyticalEwald = (getenv("GMX_OCL_NB_ANA_EWALD") != NULL);
    bForceTabulatedEwald  = (getenv("GMX_OCL_NB_TAB_EWALD") != NULL);

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
    //if ((dev_info->prop.major >= 3 || bForceAnalyticalEwald) && !bForceTabulatedEwald)
    if ((1                         || bForceAnalyticalEwald) && !bForceTabulatedEwald)
    {
        bUseAnalyticalEwald = true;

        if (debug)
        {
            fprintf(debug, "Using analytical Ewald OpenCL kernels\n");
        }
    }
    else
    {
        bUseAnalyticalEwald = false;

        if (debug)
        {
            fprintf(debug, "Using tabulated Ewald OpenCL kernels\n");
        }
    }

    /* Use twin cut-off kernels if requested by bTwinCut or the env. var.
       forces it (use it for debugging/benchmarking only). */
    if (!bTwinCut && (getenv("GMX_OCL_NB_EWALD_TWINCUT") == NULL))
    {
        kernel_type = bUseAnalyticalEwald ? eelOclEWALD_ANA : eelOclEWALD_TAB;
    }
    else
    {
        kernel_type = bUseAnalyticalEwald ? eelOclEWALD_ANA_TWIN : eelOclEWALD_TAB_TWIN;
    }

    return kernel_type;
}
