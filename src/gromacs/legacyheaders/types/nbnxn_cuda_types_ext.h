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

#ifndef NBNXN_CUDA_TYPES_EXT_H
#define NBNXN_CUDA_TYPES_EXT_H

#ifdef __cplusplus
extern "C" {
#endif

/* Abstract types */
/* CUDA nonbonded structure */
typedef struct nbnxn_cuda *nbnxn_cuda_ptr_t;
/* CUDA GPU device info */
typedef struct cuda_dev_info *cuda_dev_info_ptr_t;

/* Types defined for the structs below. */
typedef struct wallclock_gpu wallclock_gpu_t;
typedef struct nbnxn_cuda_ktime nbnxn_cuda_ktime_t;

/* Nonbonded kernel time and call count. */
struct nbnxn_cuda_ktime
{
    double  t;
    int     c;
};

/* GPU timings for kernels and H2d/D2H transfers. */
struct wallclock_gpu
{
    nbnxn_cuda_ktime_t ktime[2][2]; /* table containing the timings of the four
                                       version of the nonbonded kernels: force-only,
                                       force+energy, force+pruning, and force+energy+pruning */
    double  nb_h2d_t;               /* host to device transfer time in nb calculation  */
    double  nb_d2h_t;               /* device to host transfer time in nb calculation */
    int     nb_c;                   /* total call count of the nonbonded gpu operations */
    double  pl_h2d_t;               /* pair search step host to device transfer time */
    int     pl_h2d_c;               /* pair search step  host to device transfer call count */
};

#ifdef __cplusplus
}
#endif

#endif /* NBNXN_CUDA_TYPES_EXT_H */
