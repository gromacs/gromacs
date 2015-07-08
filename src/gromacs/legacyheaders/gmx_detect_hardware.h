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

#ifndef GMX_HARDWARE_DETECT_H
#define GMX_HARDWARE_DETECT_H

#include "gromacs/legacyheaders/typedefs.h"

#ifdef __cplusplus
extern "C" {
#endif
#if 0
} /* fixes auto-indentation problems */
#endif

/*! \brief Return whether mdrun can use more than one GPU per node
 *
 * The OpenCL implementation cannot use more than one GPU per node,
 * for example. */
gmx_bool gmx_multiple_gpu_per_node_supported();

/*! \brief Return whether PP ranks can share a GPU
 *
 * The OpenCL implementation cannot share a GPU between ranks, for
 * example. */
gmx_bool gmx_gpu_sharing_supported();

/* the init and consistency functions depend on commrec that may not be
   consistent in cuda because MPI types don't exist there.  */
#ifndef __CUDACC__
/* Construct the global hwinfo structure and return a pointer to
   it. Caller is responsible for freeing this pointer. */
gmx_hw_info_t *gmx_detect_hardware(FILE *fplog, const t_commrec *cr,
                                   gmx_bool bDetectGPUs);

/* Print information about the detected hardware to fplog (if != NULL)
 * and to stderr the master rank.
 */
void gmx_print_detected_hardware(FILE *fplog, const t_commrec *cr,
                                 const gmx_hw_info_t *hwinfo);

void gmx_hardware_info_free(gmx_hw_info_t *hwinfo);

void gmx_parse_gpu_ids(gmx_gpu_opt_t *gpu_opt);

void gmx_select_gpu_ids(FILE *fplog, const t_commrec *cr,
                        const gmx_gpu_info_t *gpu_info,
                        gmx_bool bForceUseGPU,
                        gmx_gpu_opt_t *gpu_opt);

/* Check the consistency of hw_opt with hwinfo.
   This function should be called once on each MPI rank. */
void gmx_check_hw_runconf_consistency(FILE                *fplog,
                                      const gmx_hw_info_t *hwinfo,
                                      const t_commrec     *cr,
                                      const gmx_hw_opt_t  *hw_opt,
                                      gmx_bool             bUseGPU);
#endif


/* Check whether a GPU is shared among ranks, and return the number of shared
   gpus

   gpu_opt       = the gpu options struct

   returns: The number of GPUs shared among ranks, or 0 */
int gmx_count_gpu_dev_shared(const gmx_gpu_opt_t *gpu_opt);


#ifdef __cplusplus
}
#endif


#endif /* GMX_HARDWARE_DETECT_H */
