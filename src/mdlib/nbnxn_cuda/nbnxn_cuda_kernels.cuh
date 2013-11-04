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
