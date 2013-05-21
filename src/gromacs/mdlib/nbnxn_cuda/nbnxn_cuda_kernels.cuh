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
 *  This header has the sole purpose of generating kernels for the supported
 *  electrostatics types: cut-off, reaction-field, Ewald, and tabulated Ewald.
 *
 *  The Ewald kernels have twin-range cut-off versions with rcoul != rvdw which
 *  require an extra distance check to enable  PP-PME load balancing
 *  (otherwise, by default rcoul == rvdw).
 *
 *  NOTE: No include fence as it is meant to be included multiple times.
 */

/* Analytical plain cut-off kernels */
#define EL_CUTOFF
#define NB_KERNEL_FUNC_NAME(x,...) x##_cutoff##__VA_ARGS__
#include "nbnxn_cuda_kernel_legacy.cuh"
#include "nbnxn_cuda_kernel.cuh"
#undef EL_CUTOFF
#undef NB_KERNEL_FUNC_NAME

/* Analytical reaction-field kernels */
#define EL_RF
#define NB_KERNEL_FUNC_NAME(x,...) x##_rf##__VA_ARGS__
#include "nbnxn_cuda_kernel_legacy.cuh"
#include "nbnxn_cuda_kernel.cuh"
#undef EL_RF
#undef NB_KERNEL_FUNC_NAME

/* Analytical Ewald interaction kernels
 * NOTE: no legacy kernels with analytical Ewald.
 */
#define EL_EWALD_ANA
#define NB_KERNEL_FUNC_NAME(x,...) x##_ewald##__VA_ARGS__
#include "nbnxn_cuda_kernel.cuh"
#undef EL_EWALD_ANA
#undef NB_KERNEL_FUNC_NAME

/* Analytical Ewald interaction kernels with twin-range cut-off
 * NOTE: no legacy kernels with analytical Ewald.
 */
#define EL_EWALD_ANA
#define VDW_CUTOFF_CHECK
#define NB_KERNEL_FUNC_NAME(x,...) x##_ewald_twin##__VA_ARGS__
#include "nbnxn_cuda_kernel.cuh"
#undef EL_EWALD_ANA
#undef VDW_CUTOFF_CHECK
#undef NB_KERNEL_FUNC_NAME

/* Tabulated Ewald interaction kernels */
#define EL_EWALD_TAB
#define NB_KERNEL_FUNC_NAME(x,...) x##_ewald_tab##__VA_ARGS__
#include "nbnxn_cuda_kernel_legacy.cuh"
#include "nbnxn_cuda_kernel.cuh"
#undef EL_EWALD_TAB
#undef NB_KERNEL_FUNC_NAME

/* Tabulated Ewald interaction kernels with twin-range cut-off */
#define EL_EWALD_TAB
#define VDW_CUTOFF_CHECK
#define NB_KERNEL_FUNC_NAME(x,...) x##_ewald_tab_twin##__VA_ARGS__
#include "nbnxn_cuda_kernel_legacy.cuh"
#include "nbnxn_cuda_kernel.cuh"
#undef EL_EWALD_TAB
#undef VDW_CUTOFF_CHECK
#undef NB_KERNEL_FUNC_NAME
