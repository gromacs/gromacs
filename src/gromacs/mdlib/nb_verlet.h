/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015,2017, by the GROMACS development team, led by
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

#ifndef NB_VERLET_H
#define NB_VERLET_H

#include "gromacs/mdlib/nbnxn_gpu_types.h"
#include "gromacs/mdlib/nbnxn_pairlist.h"

#ifdef __cplusplus
extern "C" {
#endif


/** Nonbonded NxN kernel types: plain C, CPU SIMD, GPU, GPU emulation */
typedef enum
{
    nbnxnkNotSet = 0,
    nbnxnk4x4_PlainC,
    nbnxnk4xN_SIMD_4xN,
    nbnxnk4xN_SIMD_2xNN,
    nbnxnk8x8x8_GPU,
    nbnxnk8x8x8_PlainC,
    nbnxnkNR
} nbnxn_kernel_type;

/** Return a string indentifying the kernel type */
const char *lookup_nbnxn_kernel_name(int kernel_type);

enum {
    ewaldexclTable, ewaldexclAnalytical
};

/* Atom locality indicator: local, non-local, all, used for calls to:
   gridding, pair-search, force calculation, x/f buffer operations */
enum {
    eatLocal = 0, eatNonlocal = 1, eatAll
};

#define LOCAL_A(x)               ((x) == eatLocal)
#define NONLOCAL_A(x)            ((x) == eatNonlocal)
#define LOCAL_OR_NONLOCAL_A(x)   (LOCAL_A(x) || NONLOCAL_A(x))

/* Interaction locality indicator (used in pair-list search/calculations):
    - local interactions require local atom data and affect local output only;
    - non-local interactions require both local and non-local atom data and
      affect both local- and non-local output. */
enum {
    eintLocal = 0, eintNonlocal = 1
};

#define LOCAL_I(x)               ((x) == eintLocal)
#define NONLOCAL_I(x)            ((x) == eintNonlocal)

enum {
    enbvClearFNo, enbvClearFYes
};

typedef struct nonbonded_verlet_group_t {
    nbnxn_pairlist_set_t  nbl_lists;   /* pair list(s)                       */
    nbnxn_atomdata_t     *nbat;        /* atom data                          */
    int                   kernel_type; /* non-bonded kernel - see enum above */
    int                   ewald_excl;  /* Ewald exclusion - see enum above   */
} nonbonded_verlet_group_t;

/* non-bonded data structure with Verlet-type cut-off */
typedef struct nonbonded_verlet_t {
    nbnxn_search_t           nbs;             /* n vs n atom pair searching data       */
    int                      ngrp;            /* number of interaction groups          */
    nonbonded_verlet_group_t grp[2];          /* local and non-local interaction group */

    gmx_bool                 bUseGPU;         /* TRUE when GPU acceleration is used */
    gmx_nbnxn_gpu_t         *gpu_nbv;         /* pointer to GPU nb verlet data     */
    int                      min_ci_balanced; /* pair list balancing parameter
                                                 used for the 8x8x8 GPU kernels    */
} nonbonded_verlet_t;

/* The minimum nstlist for dynamic pair list pruning on a CPU or a GPU.
 * In most cases going lower than 4 will lead to a too high pruning cost
 */
static const int c_nbnxnDynamicListPruningMinLifetime = 4;

/* The interval in steps at which we perform dynamic, rolling pruning on a GPU.
 * Ideally we should auto-tune this value.
 * Not considering overheads, 1 would be the ideal value. But 2 seems a reasonabe
 * compromise that reduces GPU kernel launch overheads and also avoids inefficiency
 * due to scaling on large GPUs when pruning small workloads.
 * Note that with domain decomposition we alternate between pruning the local
 * and non-local pair-lists, so 2 is a convenient value.
 */
static const int c_nbnxnGpuRollingListPruningInterval = 2;

/*! \brief Getter for bUseGPU */
gmx_bool
usingGpu(nonbonded_verlet_t *nbv);

#ifdef __cplusplus
}
#endif

#endif /* NB_VERLET_H */
