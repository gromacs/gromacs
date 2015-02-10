/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014, by the GROMACS development team, led by
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

#ifndef HWINFO_H
#define HWINFO_H

#include "gromacs/legacyheaders/gmx_cpuid.h"
#include "gromacs/legacyheaders/types/simple.h"

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
    egpuCompatible = 0,  egpuNonexistent,  egpuIncompatible, egpuInsane
} e_gpu_detect_res_t;

/* Textual names of the GPU detection/check results (see e_gpu_detect_res_t). */
static const char * const gpu_detect_res_str[] =
{
    "compatible", "inexistent", "incompatible", "insane"
};

/* GPU device information -- for now with only CUDA devices
 * The gmx_hardware_detect module initializes it. */
struct gmx_gpu_info_t
{
    gmx_bool                  bDetectGPUs;      /* Did we try to detect GPUs? */
    int                       n_dev;            /* total number of GPU devices detected */
    struct gmx_device_info_t *gpu_dev;          /* GPU devices detected in the system (per node) */
    int                       n_dev_compatible; /* number of compatible GPUs */
};

/* Hardware information structure with CPU and GPU information.
 * It is initialized by gmx_detect_hardware().
 * NOTE: this structure may only contain structures that are globally valid
 *       (i.e. must be able to be shared among all threads) */
typedef struct
{
    struct gmx_gpu_info_t gpu_info;      /* Information about GPUs detected in the system */

    gmx_cpuid_t           cpuid_info;    /* CPUID information about CPU detected;
                                            NOTE: this will only detect the CPU thread 0 of the
                                            current process runs on. */
    int             nthreads_hw_avail;   /* Number of hardware threads available; this number
                                            is based on the number of CPUs reported as available
                                            by the OS at the time of detection. */
} gmx_hw_info_t;


/* The options for the thread affinity setting, default: auto */
enum {
    threadaffSEL, threadaffAUTO, threadaffON, threadaffOFF, threadaffNR
};

/* GPU device selection information -- for now with only CUDA devices */
typedef struct
{
    char     *gpu_id;           /* GPU id's to use, each specified as chars */
    gmx_bool  bUserSet;         /* true if the GPUs in dev_use are manually provided by the user */

    int       n_dev_compatible; /* number of compatible GPU devices that could be used */
    int      *dev_compatible;   /* array of compatible GPU device IDs, from which automatic selection occurs */
    int       n_dev_use;        /* number of GPU devices selected to be used, either by the user or automatically */
    int      *dev_use;          /* array mapping from PP rank index to GPU device ID; GPU IDs can be listed multiple times when ranks share them */
} gmx_gpu_opt_t;

/* Threading and GPU options, can be set automatically or by the user */
typedef struct {
    int           nthreads_tot;        /* Total number of threads requested (TMPI) */
    int           nthreads_tmpi;       /* Number of TMPI threads requested         */
    int           nthreads_omp;        /* Number of OpenMP threads requested       */
    int           nthreads_omp_pme;    /* As nthreads_omp, but for PME only nodes  */
    int           thread_affinity;     /* Thread affinity switch, see enum above   */
    int           core_pinning_stride; /* Logical core pinning stride              */
    int           core_pinning_offset; /* Logical core pinning offset              */

    gmx_gpu_opt_t gpu_opt;             /* The GPU options                          */
} gmx_hw_opt_t;

#ifdef __cplusplus
}
#endif

#endif /* HWINFO_H */
