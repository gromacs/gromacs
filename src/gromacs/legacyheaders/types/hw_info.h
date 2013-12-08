/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
 *
 * This file is part of GROMACS.
 * Copyright (c) 2012-
 *
 * Written by the Gromacs development team under coordination of
 * David van der Spoel, Berk Hess, and Erik Lindahl.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org
 *
 * And Hey:
 * Gromacs Runs On Most of All Computer Systems
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
    gmx_bool             bDetectGPUs;          /* Did we try to detect GPUs? */
    int                  ncuda_dev;            /* total number of devices detected */
    cuda_dev_info_ptr_t  cuda_dev;             /* devices detected in the system (per node) */
    int                  ncuda_dev_compatible; /* number of compatible GPUs */
} gmx_gpu_info_t;

/* Hardware information structure with CPU and GPU information.
 * It is initialized by gmx_detect_hardware().
 * NOTE: this structure may only contain structures that are globally valid
 *       (i.e. must be able to be shared among all threads) */
typedef struct
{
    gmx_gpu_info_t  gpu_info;            /* Information about GPUs detected in the system */

    gmx_cpuid_t     cpuid_info;          /* CPUID information about CPU detected;
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
    char     *gpu_id;        /* GPU id's to use, each specified as chars */
    gmx_bool  bUserSet;      /* true if the GPUs in cuda_dev_use are manually provided by the user */

    int       ncuda_dev_use; /* number of device (IDs) selected to be used */
    int      *cuda_dev_use;  /* device index list providing GPU to PP rank mapping, GPUs can be listed multiple times when ranks share them */
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
