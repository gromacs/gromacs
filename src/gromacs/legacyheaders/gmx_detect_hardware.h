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
 * GROup of MAchos and Cynical Suckers
 */

#ifndef GMX_HARDWARE_DETECT_H
#define GMX_HARDWARE_DETECT_H

#include "types/hw_info.h"

#ifdef __cplusplus
extern "C" {
#endif
#if 0
} /* fixes auto-indentation problems */
#endif

/* the init and consistency functions depend on commrec that may not be
   consistent in cuda because MPI types don't exist there.  */
#ifndef __CUDACC__
#include "types/commrec.h"
/* return a pointer to a global hwinfo structure. */
gmx_hw_info_t *gmx_detect_hardware(FILE *fplog, const t_commrec *cr,
                                   gmx_bool bDetectGPUs);

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
