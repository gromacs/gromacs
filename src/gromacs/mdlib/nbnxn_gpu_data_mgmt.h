/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015,2017, by the GROMACS development team, led by
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
/*! \libinternal \file
 *  \brief Declare interface for GPU data transfer for NBNXN module
 *
 *  \author Mark Abraham <mark.j.abraham@gmail.com>
 *  \ingroup module_mdlib
 *  \inlibraryapi
 */

#ifndef NBNXN_GPU_DATA_MGMT_H
#define NBNXN_GPU_DATA_MGMT_H

#include "gromacs/gpu_utils/gpu_macros.h"
#include "gromacs/mdlib/nbnxn_gpu_types.h"
#include "gromacs/mdtypes/interaction_const.h"

#ifdef __cplusplus
extern "C" {
#endif

struct nonbonded_verlet_group_t;
struct nbnxn_pairlist_t;
struct nbnxn_atomdata_t;
struct NbnxnListParameters;
struct gmx_wallclock_gpu_nbnxn_t;
struct gmx_gpu_info_t;
struct gmx_device_info_t;

/** Initializes the data structures related to GPU nonbonded calculations. */
GPU_FUNC_QUALIFIER
void nbnxn_gpu_init(gmx_nbnxn_gpu_t gmx_unused            **p_nb,
                    const gmx_device_info_t gmx_unused     *deviceInfo,
                    const interaction_const_t gmx_unused   *ic,
                    const NbnxnListParameters gmx_unused   *listParams,
                    const nbnxn_atomdata_t gmx_unused      *nbat,
                    int gmx_unused                          rank,
                    /* true if both local and non-local are done on GPU */
                    gmx_bool gmx_unused                     bLocalAndNonlocal) GPU_FUNC_TERM

/** Initializes pair-list data for GPU, called at every pair search step. */
GPU_FUNC_QUALIFIER
void nbnxn_gpu_init_pairlist(gmx_nbnxn_gpu_t gmx_unused               *nb,
                             const struct nbnxn_pairlist_t gmx_unused *h_nblist,
                             int                    gmx_unused         iloc) GPU_FUNC_TERM

/** Initializes atom-data on the GPU, called at every pair search step. */
GPU_FUNC_QUALIFIER
void nbnxn_gpu_init_atomdata(gmx_nbnxn_gpu_t gmx_unused               *nb,
                             const nbnxn_atomdata_t gmx_unused        *nbat) GPU_FUNC_TERM

/*! \brief Re-generate the GPU Ewald force table, resets rlist, and update the
 *  electrostatic type switching to twin cut-off (or back) if needed.
 */
GPU_FUNC_QUALIFIER
void nbnxn_gpu_pme_loadbal_update_param(const struct nonbonded_verlet_t gmx_unused *nbv,
                                        const interaction_const_t gmx_unused       *ic,
                                        const NbnxnListParameters gmx_unused       *listParams) GPU_FUNC_TERM

/** Uploads shift vector to the GPU if the box is dynamic (otherwise just returns). */
GPU_FUNC_QUALIFIER
void nbnxn_gpu_upload_shiftvec(gmx_nbnxn_gpu_t gmx_unused               *nb,
                               const nbnxn_atomdata_t gmx_unused        *nbatom) GPU_FUNC_TERM

/** Clears GPU outputs: nonbonded force, shift force and energy. */
GPU_FUNC_QUALIFIER
void nbnxn_gpu_clear_outputs(gmx_nbnxn_gpu_t gmx_unused *nb,
                             int              gmx_unused flags) GPU_FUNC_TERM

/** Frees all GPU resources used for the nonbonded calculations. */
GPU_FUNC_QUALIFIER
void nbnxn_gpu_free(gmx_nbnxn_gpu_t gmx_unused *nb) GPU_FUNC_TERM

/** Returns the GPU timings structure or NULL if GPU is not used or timing is off. */
GPU_FUNC_QUALIFIER
struct gmx_wallclock_gpu_nbnxn_t *nbnxn_gpu_get_timings(gmx_nbnxn_gpu_t gmx_unused *nb) GPU_FUNC_TERM_WITH_RETURN(nullptr)

/** Resets nonbonded GPU timings. */
GPU_FUNC_QUALIFIER
void nbnxn_gpu_reset_timings(struct nonbonded_verlet_t gmx_unused *nbv) GPU_FUNC_TERM

/** Calculates the minimum size of proximity lists to improve SM load balance
 *  with GPU non-bonded kernels. */
GPU_FUNC_QUALIFIER
int nbnxn_gpu_min_ci_balanced(gmx_nbnxn_gpu_t gmx_unused *nb) GPU_FUNC_TERM_WITH_RETURN(-1)

/** Returns if analytical Ewald GPU kernels are used. */
GPU_FUNC_QUALIFIER
gmx_bool nbnxn_gpu_is_kernel_ewald_analytical(const gmx_nbnxn_gpu_t gmx_unused *nb) GPU_FUNC_TERM_WITH_RETURN(FALSE)

#ifdef __cplusplus
}
#endif

#endif
