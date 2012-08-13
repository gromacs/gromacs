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

#ifndef NB_VERLET_H
#define NB_VERLET_H

#include "nbnxn_pairlist.h"
#include "nbnxn_cuda_types_ext.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! Nonbonded NxN kernel types: plain C, SSE/AVX, GPU CUDA, GPU emulation, etc */
enum { nbkNotSet = 0, 
       nbk4x4_PlainC, 
       nbk4xN_X86_SIMD128,
       nbk4xN_X86_SIMD256,
       nbk8x8x8_CUDA,
       nbk8x8x8_PlainC };

/* Note that _mm_... intrinsics can be converted to either SSE or AVX
 * depending on compiler flags.
 * For gcc we check for __AVX__
 * At least a check for icc should be added (if there is a macro)
 */
static const char *nbk_name[] =
  { "not set", "plain C 4x4",
#if !(defined GMX_X86_AVX_256 || defined GMX_X86_AVX128_FMA || defined __AVX__)
#ifndef GMX_X86_SSE4_1
#ifndef GMX_DOUBLE
    "SSE2 4x4",
#else
    "SSE2 4x2",
#endif
#else
#ifndef GMX_DOUBLE
    "SSE4.1 4x4",
#else
    "SSE4.1 4x2",
#endif
#endif
#else
#ifndef GMX_DOUBLE
    "AVX-128 4x4",
#else
    "AVX-128 4x2",
#endif
#endif
#ifndef GMX_DOUBLE
    "AVX-256 4x8",
#else
    "AVX-256 4x4",
#endif
    "CUDA 8x8x8", "plain C 8x8x8" };

/* Atom locality indicator: local, non-local, all, used for calls to:
   gridding, pair-search, force calculation, x/f buffer operations */
enum { eatLocal = 0, eatNonlocal = 1, eatAll  };

#define LOCAL_A(x)               ((x) == eatLocal)
#define NONLOCAL_A(x)            ((x) == eatNonlocal)
#define LOCAL_OR_NONLOCAL_A(x)   (LOCAL_A(x) || NONLOCAL_A(x))

/* Interaction locality indicator (used in pair-list search/calculations):
    - local interactions require local atom data and affect local output only;
    - non-local interactions require both local and non-local atom data and
      affect both local- and non-local output. */
enum { eintLocal = 0, eintNonlocal = 1 };

#define LOCAL_I(x)               ((x) == eintLocal)
#define NONLOCAL_I(x)            ((x) == eintNonlocal)

enum { enbvClearFNo, enbvClearFYes };

typedef struct {
    nbnxn_pairlist_set_t nbl_lists;   /* pair list(s)                       */
    nbnxn_atomdata_t     *nbat;       /* atom data                          */
    int                  kernel_type; /* non-bonded kernel - see enum above */
} nonbonded_verlet_group_t;

/* non-bonded data structure with Verlet-type cut-off */
typedef struct {
    nbnxn_search_t           nbs;   /* n vs n atom pair searching data          */
    int                      ngrp;  /* number of interaction groups             */
    nonbonded_verlet_group_t grp[2];/* local and non-local interaction group    */

    gmx_bool         bUseGPU;          /* TRUE when GPU acceleration is used */
    nbnxn_cuda_ptr_t cu_nbv;           /* pointer to CUDA nb verlet data     */
    int              min_ci_balanced;  /* pair list balancing parameter
                                          used for the 8x8x8 CUDA kernels    */
} nonbonded_verlet_t;

#ifdef __cplusplus
}
#endif

#endif /* NB_VERLET_H */
