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

#ifndef NB_VERLET_H
#define NB_VERLET_H

#include "nbnxn_pairlist.h"
#include "nbnxn_cuda_types_ext.h"

#ifdef __cplusplus
extern "C" {
#endif


/* For testing the reference plain-C SIMD kernels, uncomment the next lines,
 * as well as the GMX_SIMD_REFERENCE_PLAIN_C define in gmx_simd_macros.h
 * The actual SIMD width is set in gmx_simd_macros.h
 * The 4xN reference kernels support 2-, 4- and 8-way SIMD.
 * The 2x(N+N) reference kernels support 8- and 16-way SIMD.
 */
/* #define GMX_NBNXN_SIMD */
/* #define GMX_NBNXN_SIMD_4XN */
/* #define GMX_NBNXN_SIMD_2XNN */


#ifdef GMX_X86_SSE2
/* Use SIMD accelerated nbnxn search and kernels */
#define GMX_NBNXN_SIMD

/* Uncomment the next line to use, slower, 128-bit SIMD with AVX-256 */
/* #define GMX_NBNXN_HALF_WIDTH_SIMD */

/* The nbnxn SIMD 4xN and 2x(N+N) kernels can be added independently.
 * Currently the 2xNN SIMD kernels only make sense with:
 *  8-way SIMD: 4x4 setup, works with AVX-256 in single precision
 * 16-way SIMD: 4x8 setup, not used, but most of the kernel code is there
 */
#define GMX_NBNXN_SIMD_4XN
#if defined GMX_X86_AVX_256 && !(defined GMX_DOUBLE || defined GMX_NBNXN_HALF_WIDTH_SIMD)
#define GMX_NBNXN_SIMD_2XNN
#endif

#endif


/*! Nonbonded NxN kernel types: plain C, CPU SIMD, GPU CUDA, GPU emulation */
typedef enum
{
    nbnxnkNotSet = 0,
    nbnxnk4x4_PlainC,
    nbnxnk4xN_SIMD_4xN,
    nbnxnk4xN_SIMD_2xNN,
    nbnxnk8x8x8_CUDA,
    nbnxnk8x8x8_PlainC,
    nbnxnkNR
} nbnxn_kernel_type;

/*! Return a string indentifying the kernel type */
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

typedef struct {
    nbnxn_pairlist_set_t  nbl_lists;   /* pair list(s)                       */
    nbnxn_atomdata_t     *nbat;        /* atom data                          */
    int                   kernel_type; /* non-bonded kernel - see enum above */
    int                   ewald_excl;  /* Ewald exclusion - see enum above   */
} nonbonded_verlet_group_t;

/* non-bonded data structure with Verlet-type cut-off */
typedef struct {
    nbnxn_search_t           nbs;             /* n vs n atom pair searching data       */
    int                      ngrp;            /* number of interaction groups          */
    nonbonded_verlet_group_t grp[2];          /* local and non-local interaction group */

    gmx_bool                 bUseGPU;         /* TRUE when GPU acceleration is used */
    nbnxn_cuda_ptr_t         cu_nbv;          /* pointer to CUDA nb verlet data     */
    int                      min_ci_balanced; /* pair list balancing parameter
                                                 used for the 8x8x8 CUDA kernels    */
} nonbonded_verlet_t;

#ifdef __cplusplus
}
#endif

#endif /* NB_VERLET_H */
