/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2010, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
 */

#ifndef _GPU_UTILS_H_
#define _GPU_UTILS_H_

#include "types/simple.h"
#include "types/hw_info.h"

#ifdef GMX_GPU
#define FUNC_TERM_INT ;
#define FUNC_TERM_VOID ;
#define FUNC_QUALIFIER
#else
#define FUNC_TERM_INT {return -1;}
#define FUNC_TERM_VOID {}
#define FUNC_QUALIFIER static
#endif

#ifdef __cplusplus
extern "C" {
#endif

FUNC_QUALIFIER
int do_quick_memtest(int dev_id) FUNC_TERM_INT

FUNC_QUALIFIER
int do_full_memtest(int dev_id) FUNC_TERM_INT

FUNC_QUALIFIER
int do_timed_memtest(int dev_id, int time_limit) FUNC_TERM_INT

FUNC_QUALIFIER
gmx_bool is_gmx_openmm_supported_gpu(int dev_id, char *gpu_name) FUNC_TERM_INT

FUNC_QUALIFIER
void detect_cuda_gpus(gmx_gpu_info_t *gpu_info) FUNC_TERM_VOID

FUNC_QUALIFIER
void pick_compatible_gpus(gmx_gpu_info_t *gpu_info) FUNC_TERM_VOID

FUNC_QUALIFIER
gmx_bool check_select_cuda_gpus(int *checkres, gmx_gpu_info_t *gpu_info,
                                const int *requested_devs, int count) FUNC_TERM_INT

FUNC_QUALIFIER
void free_gpu_info(const gmx_gpu_info_t *gpu_info) FUNC_TERM_VOID

FUNC_QUALIFIER
gmx_bool init_gpu(int mygpu, char *result_str, const gmx_gpu_info_t *gpu_info) FUNC_TERM_INT

FUNC_QUALIFIER
gmx_bool free_gpu(char *result_str) FUNC_TERM_INT

/*! \brief Returns the device ID of the GPU currently in use.*/
FUNC_QUALIFIER
int get_current_gpu_device_id(void) FUNC_TERM_INT

FUNC_QUALIFIER
int get_gpu_device_id(const gmx_gpu_info_t *gpu_info, int index) FUNC_TERM_INT

FUNC_QUALIFIER
void get_gpu_device_info_string(char *s, const gmx_gpu_info_t *gpu_info, int index) FUNC_TERM_VOID

#ifdef __cplusplus
}
#endif

#undef FUNC_TERM_INT
#undef FUNC_TERM_VOID
#undef FUNC_QUALIFIER

#endif /* _GPU_UTILS_H_ */
