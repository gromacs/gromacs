/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2010, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
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

#ifndef _GPU_UTILS_H_
#define _GPU_UTILS_H_

#include "types/simple.h"
#include "types/hw_info.h"

#ifdef GMX_GPU
#define FUNC_TERM_INT ;
#define FUNC_TERM_VOID ;
#define FUNC_QUALIFIER
#else
#define FUNC_TERM_INT {return -1; }
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
int detect_cuda_gpus(gmx_gpu_info_t *gpu_info, char *err_str) FUNC_TERM_INT

FUNC_QUALIFIER
void pick_compatible_gpus(const gmx_gpu_info_t *gpu_info,
                          gmx_gpu_opt_t *gpu_opt) FUNC_TERM_VOID

FUNC_QUALIFIER
gmx_bool check_selected_cuda_gpus(int *checkres,
                                  const gmx_gpu_info_t *gpu_info,
                                  gmx_gpu_opt_t *gpu_opt) FUNC_TERM_INT

FUNC_QUALIFIER
void free_gpu_info(const gmx_gpu_info_t *gpu_info) FUNC_TERM_VOID

FUNC_QUALIFIER
gmx_bool init_gpu(int mygpu, char *result_str,
                  const gmx_gpu_info_t *gpu_info,
                  const gmx_gpu_opt_t *gpu_opt) FUNC_TERM_INT

FUNC_QUALIFIER
gmx_bool free_gpu(char *result_str) FUNC_TERM_INT

/*! \brief Returns the device ID of the GPU currently in use.*/
FUNC_QUALIFIER
int get_current_gpu_device_id(void) FUNC_TERM_INT

FUNC_QUALIFIER
int get_gpu_device_id(const gmx_gpu_info_t *gpu_info,
                      const gmx_gpu_opt_t *gpu_opt,
                      int index) FUNC_TERM_INT

FUNC_QUALIFIER
void get_gpu_device_info_string(char *s, const gmx_gpu_info_t *gpu_info, int index) FUNC_TERM_VOID

FUNC_QUALIFIER
size_t sizeof_cuda_dev_info(void) FUNC_TERM_INT

#ifdef __cplusplus
}
#endif

#undef FUNC_TERM_INT
#undef FUNC_TERM_VOID
#undef FUNC_QUALIFIER

#endif /* _GPU_UTILS_H_ */
