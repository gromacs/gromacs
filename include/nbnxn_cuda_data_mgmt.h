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

#ifdef __cplusplus
extern "C" {
#endif

/*! Initializes the data structures related to CUDA nonbonded calculations. */
void nbnxn_cuda_init(FILE * /*fplog*/,
                     nbnxn_cuda_ptr_t * /*p_cu_nb*/,
                     gmx_bool /*bDomDec*/);

/*! Initilizes simulation constant data. */
void nbnxn_cuda_init_const(nbnxn_cuda_ptr_t /*p_cu_nb*/,
                           const interaction_const_t * /*ic*/,
                           const nonbonded_verlet_t * /*nbv*/);

/*! Initilizes pair-list data for GPU, called at every pair search step. */
void nbnxn_cuda_init_pairlist(nbnxn_cuda_ptr_t /*cu_nb*/,
                              const nbnxn_pairlist_t * /*h_nblist*/,
                              int /*iloc*/);

/*! Initilizes atom-data on the GPU, called at every pair search step. */
void nbnxn_cuda_init_atomdata(nbnxn_cuda_ptr_t /*cu_nb*/,
                              const nbnxn_atomdata_t * /*atomdata*/);

/*! Re-generates the GPU Ewald force table and resets rlist - used with PME auto-tuning. */
void reset_gpu_rlist_ewaldtab(nbnxn_cuda_ptr_t /*cu_nb*/,
                              const interaction_const_t * /*ic*/);

/*! Uploads shift vector to the GPU if the box is dynamic (otherwise just returns). */
void nbnxn_cuda_upload_shiftvec(nbnxn_cuda_ptr_t /*cu_nb*/,
                                const nbnxn_atomdata_t * /*nbatom*/);

/*! Clears GPU outputs: nonbonded force, shift force and energy. */
void nbnxn_cuda_clear_outputs(nbnxn_cuda_ptr_t /*cu_nb*/,
                              int /*flags*/);

/*! Frees all GPU resources used for the nonbonded calculations. */
void nbnxn_cuda_free(FILE * /*fplog*/,
                     nbnxn_cuda_ptr_t /*cu_nb*/,
                     gmx_bool /*bDomDec*/);

/*! Returns the GPU timings structure or NULL if GPU is not used or timing is off. */
wallclock_gpu_t * nbnxn_cuda_get_timings(nbnxn_cuda_ptr_t /*cu_nb*/);

/*! Resets nonbonded GPU timings. */
void nbnxn_cuda_reset_timings(nbnxn_cuda_ptr_t /*cu_nb*/);

/*! Calculates the minimum size of proximity lists to improve SM load balance 
    with CUDA non-bonded kernels. */
int nbnxn_cuda_min_ci_balanced(nbnxn_cuda_ptr_t /*cu_nb*/);

#ifdef __cplusplus
}
#endif

#endif /* NBNXN_CUDA_DATA_MGMT_H */
