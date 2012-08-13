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

#include "types/simple.h"
#include "types/nbnxn_cuda_types_ext.h"
#include "gmx_detectcpu.h"

#ifdef __cplusplus
extern "C" {
#endif
#if 0
} /* fixes auto-indentation problems */
#endif

typedef enum
{
    egpuCompatible = 0,  egpuNonexistent,  egpuIncompatible, egpuInsane
} e_gpu_detect_res_t;


static const char * const gpu_detect_res_str[] =
{
    "compatible", "inexistent", "incompatible", "insane"
};

/* GPU device information -- for now with only CUDA devices.
 * The gmx_hardware_detect module initializes it. */
typedef struct 
{
    gmx_bool            bUserSet;       /* true if the GPUs in cuda_dev_use are manually provided by the user */

    int                 ncuda_dev_use;  /* number of devices selected to be used */
    int                 *cuda_dev_use;  /* index of the devices selected to be used */
    int                 ncuda_dev;      /* total number of devices detected */
    cuda_dev_info_ptr_t cuda_dev;       /* devices detected in the system (per node) */
} gmx_gpu_info_t;

typedef struct
{
    gmx_bool        bCanUseGPU; /* TODO ATM this is false if -nb cpu is passed. Do we want that? */
    gmx_gpu_info_t  gpu_info;   /* TODO ATM no GPUs are detected if -nb cpu is passed. Do we wan that? */
    gmx_detectcpu_t cpu_info;
} gmx_hwinfo_t;

/* FIXME either rename gmx_hwinfo_t to gmx_hw_info_t or gmx_gpu_info_t to gmx_gpuinfo_t */

#ifdef __cplusplus
}
#endif

#endif /* HWINFO_H */
