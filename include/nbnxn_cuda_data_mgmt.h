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
 * Copyright (c) 2001-2012, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 *
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

#ifndef NBNXN_CUDA_DATA_MGMT_H
#define NBNXN_CUDA_DATA_MGMT_H

#include "types/simple.h"
#include "types/interaction_const.h"
#include "types/nbnxn_cuda_types_ext.h"
#include "types/hwinfo.h"

#ifdef GMX_GPU
#define FUNC_TERM ;
#define FUNC_QUALIFIER
#else
#define FUNC_TERM {}
#define FUNC_QUALIFIER static
#endif

#ifdef __cplusplus
extern "C" {
#endif

/*! Initializes the data structures related to CUDA nonbonded calculations. */
FUNC_QUALIFIER
void nbnxn_cuda_init(FILE *fplog,
                     nbnxn_cuda_ptr_t *p_cu_nb,
                     gmx_gpu_info_t *gpu_info, int my_gpu_index,
                     /* true of both local and non-local are don on GPU */
                     gmx_bool bLocalAndNonlocal) FUNC_TERM

/*! Initilizes simulation constant data. */
FUNC_QUALIFIER
void nbnxn_cuda_init_const(nbnxn_cuda_ptr_t p_cu_nb,
                           const interaction_const_t *ic,
                           const nonbonded_verlet_t *nbv) FUNC_TERM

/*! Initilizes pair-list data for GPU, called at every pair search step. */
FUNC_QUALIFIER
void nbnxn_cuda_init_pairlist(nbnxn_cuda_ptr_t cu_nb,
                              const nbnxn_pairlist_t *h_nblist,
                              int iloc) FUNC_TERM

/*! Initilizes atom-data on the GPU, called at every pair search step. */
FUNC_QUALIFIER
void nbnxn_cuda_init_atomdata(nbnxn_cuda_ptr_t cu_nb,
                              const nbnxn_atomdata_t *atomdata) FUNC_TERM

/*! Re-generates the GPU Ewald force table and resets rlist - used with PME auto-tuning. */
FUNC_QUALIFIER
void reset_gpu_rlist_ewaldtab(nbnxn_cuda_ptr_t cu_nb,
                              const interaction_const_t *ic) FUNC_TERM

/*! Uploads shift vector to the GPU if the box is dynamic (otherwise just returns). */
FUNC_QUALIFIER
void nbnxn_cuda_upload_shiftvec(nbnxn_cuda_ptr_t cu_nb,
                                const nbnxn_atomdata_t *nbatom) FUNC_TERM

/*! Clears GPU outputs: nonbonded force, shift force and energy. */
FUNC_QUALIFIER
void nbnxn_cuda_clear_outputs(nbnxn_cuda_ptr_t cu_nb,
                              int flags) FUNC_TERM

/*! Frees all GPU resources used for the nonbonded calculations. */
FUNC_QUALIFIER
void nbnxn_cuda_free(FILE *fplog,
                     nbnxn_cuda_ptr_t cu_nb) FUNC_TERM

/*! Returns the GPU timings structure or NULL if GPU is not used or timing is off. */
FUNC_QUALIFIER
wallclock_gpu_t * nbnxn_cuda_get_timings(nbnxn_cuda_ptr_t cu_nb)
#ifdef GMX_GPU
;
#else
{ return NULL; }
#endif

/*! Resets nonbonded GPU timings. */
FUNC_QUALIFIER
void nbnxn_cuda_reset_timings(nbnxn_cuda_ptr_t cu_nb) FUNC_TERM

/*! Calculates the minimum size of proximity lists to improve SM load balance 
    with CUDA non-bonded kernels. */
FUNC_QUALIFIER
int nbnxn_cuda_min_ci_balanced(nbnxn_cuda_ptr_t cu_nb)
#ifdef GMX_GPU
;
#else
{ return -1; }
#endif

#ifdef __cplusplus
}
#endif

#undef FUNC_TERM
#undef FUNC_QUALIFIER

#endif /* NBNXN_CUDA_DATA_MGMT_H */
