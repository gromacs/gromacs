/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014, by the GROMACS development team, led by
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

/*! \internal \file
 *  This header has the sole purpose of generating kernels for the combinations of
 *  supported electrostatics types (cut-off, reaction-field, analytical and
 *  tabulated Ewald) and VDW types ( V shift, F switch, V swtich).
 *
 *  The Ewald kernels have twin-range cut-off versions with rcoul != rvdw which
 *  require an extra distance check to enable  PP-PME load balancing
 *  (otherwise, by default rcoul == rvdw).
 *
 *  NOTE: No include fence as it is meant to be included multiple times.
 */

/* Analytical plain cut-off electrostatics kernels
 */
#define EL_CUTOFF

/* V shift */
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecCut_VdwLJ ## __VA_ARGS__
#include "nbnxn_cuda_kernel.cuh"
#undef NB_KERNEL_FUNC_NAME
/* F switch */
#define VDW_FORCE_SWITCH
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecCut_VdwLJFsw ## __VA_ARGS__
#include "nbnxn_cuda_kernel.cuh"
#undef VDW_FORCE_SWITCH
#undef NB_KERNEL_FUNC_NAME
/* V switch */
#define VDW_POT_SWITCH
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecCut_VdwLJPsw ## __VA_ARGS__
#include "nbnxn_cuda_kernel.cuh"
#undef VDW_POT_SWITCH
#undef NB_KERNEL_FUNC_NAME

#undef EL_CUTOFF

/* Analytical reaction-field kernels
 */
#define EL_RF

/* V shift */
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecRF_VdwLJ ## __VA_ARGS__
#include "nbnxn_cuda_kernel.cuh"
#undef NB_KERNEL_FUNC_NAME
/* F switch */
#define VDW_FORCE_SWITCH
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecRF_VdwLJFsw ## __VA_ARGS__
#include "nbnxn_cuda_kernel.cuh"
#undef VDW_FORCE_SWITCH
#undef NB_KERNEL_FUNC_NAME
/* V switch */
#define VDW_POT_SWITCH
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecRF_VdwLJPsw ## __VA_ARGS__
#include "nbnxn_cuda_kernel.cuh"
#undef VDW_POT_SWITCH
#undef NB_KERNEL_FUNC_NAME

#undef EL_RF


/* Analytical Ewald interaction kernels
 */
#define EL_EWALD_ANA

/* V shift */
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecEw_VdwLJ ## __VA_ARGS__
#include "nbnxn_cuda_kernel.cuh"
#undef NB_KERNEL_FUNC_NAME
/* F switch */
#define VDW_FORCE_SWITCH
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecEw_VdwLJFsw ## __VA_ARGS__
#include "nbnxn_cuda_kernel.cuh"
#undef VDW_FORCE_SWITCH
#undef NB_KERNEL_FUNC_NAME
/* V switch */
#define VDW_POT_SWITCH
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecEw_VdwLJPsw ## __VA_ARGS__
#include "nbnxn_cuda_kernel.cuh"
#undef VDW_POT_SWITCH
#undef NB_KERNEL_FUNC_NAME

#undef EL_EWALD_ANA



/* Analytical Ewald interaction kernels with twin-range cut-off
 */
#define EL_EWALD_ANA
#define VDW_CUTOFF_CHECK

/* V shift */
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecEwTwinCut_VdwLJ ## __VA_ARGS__
#include "nbnxn_cuda_kernel.cuh"
#undef NB_KERNEL_FUNC_NAME
/* F switch */
#define VDW_FORCE_SWITCH
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecEwTwinCut_VdwLJFsw ## __VA_ARGS__
#include "nbnxn_cuda_kernel.cuh"
#undef VDW_FORCE_SWITCH
#undef NB_KERNEL_FUNC_NAME
/* V switch */
#define VDW_POT_SWITCH
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecEwTwinCut_VdwLJPsw ## __VA_ARGS__
#include "nbnxn_cuda_kernel.cuh"
#undef VDW_POT_SWITCH
#undef NB_KERNEL_FUNC_NAME

#undef EL_EWALD_ANA
#undef VDW_CUTOFF_CHECK



/* Tabulated Ewald interaction kernels */
#define EL_EWALD_TAB

/* V shift */
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecEwQSTab_VdwLJ ## __VA_ARGS__
#include "nbnxn_cuda_kernel.cuh"
#undef NB_KERNEL_FUNC_NAME
/* F switch */
#define VDW_FORCE_SWITCH
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecEwQSTab_VdwLJFsw ## __VA_ARGS__
#include "nbnxn_cuda_kernel.cuh"
#undef VDW_FORCE_SWITCH
#undef NB_KERNEL_FUNC_NAME
/* V switch */
#define VDW_POT_SWITCH
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecEwQSTab_VdwLJPsw ## __VA_ARGS__
#include "nbnxn_cuda_kernel.cuh"
#undef VDW_POT_SWITCH
#undef NB_KERNEL_FUNC_NAME

#undef EL_EWALD_TAB


/* Tabulated Ewald interaction kernels with twin-range cut-off */
#define EL_EWALD_TAB
#define VDW_CUTOFF_CHECK

/* V shift */
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecEwQSTabTwinCut_VdwLJ ## __VA_ARGS__
#include "nbnxn_cuda_kernel.cuh"
#undef NB_KERNEL_FUNC_NAME
/* F switch */
#define VDW_FORCE_SWITCH
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecEwQSTabTwinCut_VdwLJFsw ## __VA_ARGS__
#include "nbnxn_cuda_kernel.cuh"
#undef VDW_FORCE_SWITCH
#undef NB_KERNEL_FUNC_NAME
/* V switch */
#define VDW_POT_SWITCH
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecEwQSTabTwinCut_VdwLJPsw ## __VA_ARGS__
#include "nbnxn_cuda_kernel.cuh"
#undef VDW_POT_SWITCH
#undef NB_KERNEL_FUNC_NAME

#undef EL_EWALD_TAB
#undef VDW_CUTOFF_CHECK
