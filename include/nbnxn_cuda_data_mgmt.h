/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2012, The GROMACS development team,
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

#ifndef NBNXN_CUDA_DATA_MGMT_H
#define NBNXN_CUDA_DATA_MGMT_H

#include "types/simple.h"
#include "types/interaction_const.h"
#include "types/nbnxn_cuda_types_ext.h"
#include "types/hw_info.h"

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
                     const gmx_gpu_info_t *gpu_info,
                     const gmx_gpu_opt_t *gpu_opt,
                     int my_gpu_index,
                     /* true of both local and non-local are don on GPU */
                     gmx_bool bLocalAndNonlocal) FUNC_TERM

/*! Initializes simulation constant data. */
FUNC_QUALIFIER
void nbnxn_cuda_init_const(nbnxn_cuda_ptr_t                cu_nb,
                           const interaction_const_t      *ic,
                           const nonbonded_verlet_group_t *nbv_group) FUNC_TERM

/*! Initializes pair-list data for GPU, called at every pair search step. */
FUNC_QUALIFIER
void nbnxn_cuda_init_pairlist(nbnxn_cuda_ptr_t        cu_nb,
                              const nbnxn_pairlist_t *h_nblist,
                              int                     iloc) FUNC_TERM

/*! Initializes atom-data on the GPU, called at every pair search step. */
FUNC_QUALIFIER
void nbnxn_cuda_init_atomdata(nbnxn_cuda_ptr_t        cu_nb,
                              const nbnxn_atomdata_t *atomdata) FUNC_TERM

/*! \brief Update parameters during PP-PME load balancing. */
FUNC_QUALIFIER
void nbnxn_cuda_pme_loadbal_update_param(nbnxn_cuda_ptr_t           cu_nb,
                                         const interaction_const_t *ic) FUNC_TERM

/*! Uploads shift vector to the GPU if the box is dynamic (otherwise just returns). */
FUNC_QUALIFIER
void nbnxn_cuda_upload_shiftvec(nbnxn_cuda_ptr_t        cu_nb,
                                const nbnxn_atomdata_t *nbatom) FUNC_TERM

/*! Clears GPU outputs: nonbonded force, shift force and energy. */
FUNC_QUALIFIER
void nbnxn_cuda_clear_outputs(nbnxn_cuda_ptr_t cu_nb,
                              int              flags) FUNC_TERM

/*! Frees all GPU resources used for the nonbonded calculations. */
FUNC_QUALIFIER
void nbnxn_cuda_free(FILE            *fplog,
                     nbnxn_cuda_ptr_t cu_nb) FUNC_TERM

/*! Returns the GPU timings structure or NULL if GPU is not used or timing is off. */
FUNC_QUALIFIER
wallclock_gpu_t * nbnxn_cuda_get_timings(nbnxn_cuda_ptr_t cu_nb)
#ifdef GMX_GPU
;
#else
{
    return NULL;
}
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
{
    return -1;
}
#endif

/*! Returns if analytical Ewald CUDA kernels are used. */
FUNC_QUALIFIER
gmx_bool nbnxn_cuda_is_kernel_ewald_analytical(const nbnxn_cuda_ptr_t cu_nb)
#ifdef GMX_GPU
;
#else
{
    return FALSE;
}
#endif

#ifdef __cplusplus
}
#endif

#undef FUNC_TERM
#undef FUNC_QUALIFIER

#endif /* NBNXN_CUDA_DATA_MGMT_H */
