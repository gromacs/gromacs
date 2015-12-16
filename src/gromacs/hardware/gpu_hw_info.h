/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015, by the GROMACS development team, led by
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
#ifndef GMX_HARDWARE_GPU_HW_INFO_H
#define GMX_HARDWARE_GPU_HW_INFO_H

#include "gromacs/utility/basedefinitions.h"

#ifdef __cplusplus
extern "C" {
#endif
#if 0
} /* fixes auto-indentation problems */
#endif

struct gmx_device_info_t;

/* Possible results of the GPU detection/check.
 *
 * The egpuInsane value means that during the sanity checks an error
 * occurred that indicates malfunctioning of the device, driver, or
 * incompatible driver/runtime. */
typedef enum
{
    egpuCompatible = 0,  egpuNonexistent,  egpuIncompatible, egpuInsane, egpuNR
} e_gpu_detect_res_t;

/* Names of the GPU detection/check results */
extern const char * const gpu_detect_res_str[egpuNR];

/* GPU device information -- includes either CUDA or OpenCL devices.
 * The gmx_hardware_detect module initializes it. */
struct gmx_gpu_info_t
{
    gmx_bool                  bDetectGPUs;      /* Did we try to detect GPUs? */
    int                       n_dev;            /* total number of GPU devices detected */
    struct gmx_device_info_t *gpu_dev;          /* GPU devices detected in the system (per node) */
    int                       n_dev_compatible; /* number of compatible GPUs */
};

/* GPU device selection information -- includes either CUDA or OpenCL devices */
typedef struct gmx_gpu_opt_t
{
    char     *gpu_id;           /* GPU id's to use, each specified as chars */
    gmx_bool  bUserSet;         /* true if the GPUs in dev_use are manually provided by the user */

    int       n_dev_compatible; /* number of compatible GPU devices that could be used */
    int      *dev_compatible;   /* array of compatible GPU device IDs, from which automatic selection occurs */
    int       n_dev_use;        /* number of GPU devices selected to be used, either by the user or automatically */
    int      *dev_use;          /* array mapping from PP rank index to GPU device ID; GPU IDs can be listed multiple times when ranks share them */
} gmx_gpu_opt_t;

#ifdef __cplusplus
}
#endif

#endif
