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

/*! \file
 *  This header has the sole purporse of generating kernels for the different
 *  type of electrostatics supported: Cut-off, Reaction-Field, and Ewald/PME.
 *  (No include fence as it is meant to be included multiple times.)
 */

/* Cut-Off */
#define EL_CUTOFF
#define NB_KERNEL_FUNC_NAME(x,...) x##_cutoff##__VA_ARGS__
#include "nbnxn_cuda_kernel_old.cuh"
#include "nbnxn_cuda_kernel_legacy.cuh"
#include "nbnxn_cuda_kernel.cuh"
#undef EL_CUTOFF
#undef NB_KERNEL_FUNC_NAME

/* Reaction-Field */
#define EL_RF
#define NB_KERNEL_FUNC_NAME(x,...) x##_rf##__VA_ARGS__
#include "nbnxn_cuda_kernel_old.cuh"
#include "nbnxn_cuda_kernel_legacy.cuh"
#include "nbnxn_cuda_kernel.cuh"
#undef EL_RF
#undef NB_KERNEL_FUNC_NAME

/* Ewald */
#define EL_EWALD
#define NB_KERNEL_FUNC_NAME(x,...) x##_ewald##__VA_ARGS__
#include "nbnxn_cuda_kernel_old.cuh"
#include "nbnxn_cuda_kernel_legacy.cuh"
#include "nbnxn_cuda_kernel.cuh"
#undef EL_EWALD
#undef NB_KERNEL_FUNC_NAME
