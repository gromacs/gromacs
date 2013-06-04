/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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

#ifndef HWINFO_H
#define HWINFO_H

#include "simple.h"
#include "nbnxn_cuda_types_ext.h"
#include "../gmx_cpuid.h"

#ifdef __cplusplus
extern "C" {
#endif
#if 0
} /* fixes auto-indentation problems */
#endif

/* Possible results of the GPU detection/check.
 *
 * The egpuInsane value means that during the sanity checks an error
 * occurred that indicates malfunctioning of the device, driver, or
 * incompatible driver/runtime. */
typedef enum
{
    egpuCompatible = 0,  egpuNonexistent,  egpuIncompatible, egpuInsane
} e_gpu_detect_res_t;

/* Textual names of the GPU detection/check results (see e_gpu_detect_res_t). */
static const char * const gpu_detect_res_str[] =
{
    "compatible", "inexistent", "incompatible", "insane"
};

/* GPU device information -- for now with only CUDA devices.
 * The gmx_hardware_detect module initializes it. */
typedef struct
{
    gmx_bool             bUserSet;      /* true if the GPUs in cuda_dev_use are manually provided by the user */

    int                  ncuda_dev_use; /* number of devices selected to be used */
    int                 *cuda_dev_use;  /* index of the devices selected to be used */
    int                  ncuda_dev;     /* total number of devices detected */
    cuda_dev_info_ptr_t  cuda_dev;      /* devices detected in the system (per node) */
} gmx_gpu_info_t;

/* Hardware information structure with CPU and GPU information.
 * It is initialized by gmx_detect_hardware().
 * NOTE: this structure may only contain structures that are globally valid
 *       (i.e. must be able to be shared among all threads) */
typedef struct
{
    gmx_bool        bCanUseGPU;          /* True if compatible GPUs are detected during hardware detection */
    gmx_gpu_info_t  gpu_info;            /* Information about GPUs detected in the system */

    gmx_cpuid_t     cpuid_info;          /* CPUID information about CPU detected;
                                            NOTE: this will only detect the CPU thread 0 of the
                                            current process runs on. */
    int             nthreads_hw_avail;   /* Number of hardware threads available; this number
                                            is based on the number of CPUs reported as available
                                            by the OS at the time of detection. */
    gmx_bool        bConsistencyChecked; /* whether
                                            gmx_check_hw_runconf_consistency()
                                            has been run with this hw_info */
} gmx_hw_info_t;

#ifdef __cplusplus
}
#endif

#endif /* HWINFO_H */
