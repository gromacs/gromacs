/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015, by the GROMACS development team, led by
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
 *  supported electrostatics types (cut-off, reaction-field, analytical,
 *  tabulated Ewald, tabulated GENERIC, No electrostatics) and VDW types
 * (cut-off + V shift, LJ-Ewald with geometric or Lorentz-Berthelot combination
 * rule, F switch, V switch, Tabulated LJ6+12, Tabulated GENERIC).
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

/* cut-off + V shift LJ */
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecCut_VdwLJ ## __VA_ARGS__
#include "gromacs/mdlib/nbnxn_cuda/nbnxn_cuda_kernel.cuh"
#undef NB_KERNEL_FUNC_NAME
/* LJ-Ewald w geometric combination rules */
#define LJ_EWALD_COMB_GEOM
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecCut_VdwLJEwCombGeom ## __VA_ARGS__
#include "gromacs/mdlib/nbnxn_cuda/nbnxn_cuda_kernel.cuh"
#undef LJ_EWALD_COMB_GEOM
#undef NB_KERNEL_FUNC_NAME
/* LJ-Ewald w LB combination rules */
#define LJ_EWALD_COMB_LB
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecCut_VdwLJEwCombLB ## __VA_ARGS__
#include "gromacs/mdlib/nbnxn_cuda/nbnxn_cuda_kernel.cuh"
#undef LJ_EWALD_COMB_LB
#undef NB_KERNEL_FUNC_NAME
/* F switch LJ */
#define LJ_FORCE_SWITCH
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecCut_VdwLJFsw ## __VA_ARGS__
#include "gromacs/mdlib/nbnxn_cuda/nbnxn_cuda_kernel.cuh"
#undef LJ_FORCE_SWITCH
#undef NB_KERNEL_FUNC_NAME
/* V switch LJ */
#define LJ_POT_SWITCH
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecCut_VdwLJPsw ## __VA_ARGS__
#include "gromacs/mdlib/nbnxn_cuda/nbnxn_cuda_kernel.cuh"
#undef LJ_POT_SWITCH
#undef NB_KERNEL_FUNC_NAME
/* Vdw USER tables LJ6+LJ12 */
#define VDW_USER
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecCut_Vdw_USER ## __VA_ARGS__
#include "gromacs/mdlib/nbnxn_cuda/nbnxn_cuda_kernel.cuh"
#undef VDW_USER
#undef NB_KERNEL_FUNC_NAME
/* Vdw USER tables GENERIC */
#define VDW_USER_GENERIC
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecCut_Vdw_USER_GENERIC ## __VA_ARGS__
#include "gromacs/mdlib/nbnxn_cuda/nbnxn_cuda_kernel.cuh"
#undef VDW_USER_GENERIC
#undef NB_KERNEL_FUNC_NAME

#undef EL_CUTOFF


/* Analytical reaction-field kernels
 */
#define EL_RF

/* cut-off + V shift LJ */
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecRF_VdwLJ ## __VA_ARGS__
#include "gromacs/mdlib/nbnxn_cuda/nbnxn_cuda_kernel.cuh"
#undef NB_KERNEL_FUNC_NAME
/* LJ-Ewald w geometric combination rules */
#define LJ_EWALD_COMB_GEOM
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecRF_VdwLJEwCombGeom ## __VA_ARGS__
#include "gromacs/mdlib/nbnxn_cuda/nbnxn_cuda_kernel.cuh"
#undef LJ_EWALD_COMB_GEOM
#undef NB_KERNEL_FUNC_NAME
/* LJ-Ewald w LB combination rules */
#define LJ_EWALD_COMB_LB
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecRF_VdwLJEwCombLB ## __VA_ARGS__
#include "gromacs/mdlib/nbnxn_cuda/nbnxn_cuda_kernel.cuh"
#undef LJ_EWALD_COMB_LB
#undef NB_KERNEL_FUNC_NAME
/* F switch LJ */
#define LJ_FORCE_SWITCH
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecRF_VdwLJFsw ## __VA_ARGS__
#include "gromacs/mdlib/nbnxn_cuda/nbnxn_cuda_kernel.cuh"
#undef LJ_FORCE_SWITCH
#undef NB_KERNEL_FUNC_NAME
/* V switch LJ */
#define LJ_POT_SWITCH
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecRF_VdwLJPsw ## __VA_ARGS__
#include "gromacs/mdlib/nbnxn_cuda/nbnxn_cuda_kernel.cuh"
#undef LJ_POT_SWITCH
#undef NB_KERNEL_FUNC_NAME
/* Vdw USER tables LJ6+LJ12 */
#define VDW_USER
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecRF_Vdw_USER ## __VA_ARGS__
#include "gromacs/mdlib/nbnxn_cuda/nbnxn_cuda_kernel.cuh"
#undef VDW_USER
#undef NB_KERNEL_FUNC_NAME
/* Vdw USER tables GENERIC */
#define VDW_USER_GENERIC
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecRF_Vdw_USER_GENERIC ## __VA_ARGS__
#include "gromacs/mdlib/nbnxn_cuda/nbnxn_cuda_kernel.cuh"
#undef VDW_USER_GENERIC
#undef NB_KERNEL_FUNC_NAME

#undef EL_RF


/* Analytical Ewald interaction kernels
 */
#define EL_EWALD_ANA

/* cut-off + V shift LJ */
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecEw_VdwLJ ## __VA_ARGS__
#include "gromacs/mdlib/nbnxn_cuda/nbnxn_cuda_kernel.cuh"
#undef NB_KERNEL_FUNC_NAME
/* LJ-Ewald w geometric combination rules */
#define LJ_EWALD_COMB_GEOM
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecEw_VdwLJEwCombGeom ## __VA_ARGS__
#include "gromacs/mdlib/nbnxn_cuda/nbnxn_cuda_kernel.cuh"
#undef LJ_EWALD_COMB_GEOM
#undef NB_KERNEL_FUNC_NAME
/* LJ-Ewald w LB combination rules */
#define LJ_EWALD_COMB_LB
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecEw_VdwLJEwCombLB ## __VA_ARGS__
#include "gromacs/mdlib/nbnxn_cuda/nbnxn_cuda_kernel.cuh"
#undef LJ_EWALD_COMB_LB
#undef NB_KERNEL_FUNC_NAME
/* F switch LJ */
#define LJ_FORCE_SWITCH
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecEw_VdwLJFsw ## __VA_ARGS__
#include "gromacs/mdlib/nbnxn_cuda/nbnxn_cuda_kernel.cuh"
#undef LJ_FORCE_SWITCH
#undef NB_KERNEL_FUNC_NAME
/* V switch LJ */
#define LJ_POT_SWITCH
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecEw_VdwLJPsw ## __VA_ARGS__
#include "gromacs/mdlib/nbnxn_cuda/nbnxn_cuda_kernel.cuh"
#undef LJ_POT_SWITCH
#undef NB_KERNEL_FUNC_NAME
/* Vdw USER tables LJ6+LJ12 */
#define VDW_USER
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecEw_Vdw_USER ## __VA_ARGS__
#include "gromacs/mdlib/nbnxn_cuda/nbnxn_cuda_kernel.cuh"
#undef VDW_USER
#undef NB_KERNEL_FUNC_NAME
/* Vdw USER tables GENERIC */
#define VDW_USER_GENERIC
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecEw_Vdw_USER_GENERIC ## __VA_ARGS__
#include "gromacs/mdlib/nbnxn_cuda/nbnxn_cuda_kernel.cuh"
#undef VDW_USER_GENERIC
#undef NB_KERNEL_FUNC_NAME

#undef EL_EWALD_ANA


/* Analytical Ewald interaction kernels with twin-range cut-off
 */
#define EL_EWALD_ANA
#define VDW_CUTOFF_CHECK

/* cut-off + V shift LJ */
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecEwTwinCut_VdwLJ ## __VA_ARGS__
#include "gromacs/mdlib/nbnxn_cuda/nbnxn_cuda_kernel.cuh"
#undef NB_KERNEL_FUNC_NAME
/* LJ-Ewald w geometric combination rules */
#define LJ_EWALD_COMB_GEOM
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecEwTwinCut_VdwLJEwCombGeom ## __VA_ARGS__
#include "gromacs/mdlib/nbnxn_cuda/nbnxn_cuda_kernel.cuh"
#undef LJ_EWALD_COMB_GEOM
#undef NB_KERNEL_FUNC_NAME
/* LJ-Ewald w LB combination rules */
#define LJ_EWALD_COMB_LB
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecEwTwinCut_VdwLJEwCombLB ## __VA_ARGS__
#include "gromacs/mdlib/nbnxn_cuda/nbnxn_cuda_kernel.cuh"
#undef LJ_EWALD_COMB_LB
#undef NB_KERNEL_FUNC_NAME
/* F switch LJ */
#define LJ_FORCE_SWITCH
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecEwTwinCut_VdwLJFsw ## __VA_ARGS__
#include "gromacs/mdlib/nbnxn_cuda/nbnxn_cuda_kernel.cuh"
#undef LJ_FORCE_SWITCH
#undef NB_KERNEL_FUNC_NAME
/* V switch LJ */
#define LJ_POT_SWITCH
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecEwTwinCut_VdwLJPsw ## __VA_ARGS__
#include "gromacs/mdlib/nbnxn_cuda/nbnxn_cuda_kernel.cuh"
#undef LJ_POT_SWITCH
#undef NB_KERNEL_FUNC_NAME
/* Vdw USER tables LJ6+LJ12 */
#define VDW_USER
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecEwTwinCut_Vdw_USER ## __VA_ARGS__
#include "gromacs/mdlib/nbnxn_cuda/nbnxn_cuda_kernel.cuh"
#undef VDW_USER
#undef NB_KERNEL_FUNC_NAME
/* Vdw USER tables GENERIC */
#define VDW_USER_GENERIC
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecEwTwinCut_Vdw_USER_GENERIC ## __VA_ARGS__
#include "gromacs/mdlib/nbnxn_cuda/nbnxn_cuda_kernel.cuh"
#undef VDW_USER_GENERIC
#undef NB_KERNEL_FUNC_NAME

#undef EL_EWALD_ANA
#undef VDW_CUTOFF_CHECK


/* Tabulated Ewald interaction kernels */
#define EL_EWALD_TAB

/* cut-off + V shift LJ */
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecEwQSTab_VdwLJ ## __VA_ARGS__
#include "gromacs/mdlib/nbnxn_cuda/nbnxn_cuda_kernel.cuh"
#undef NB_KERNEL_FUNC_NAME
/* LJ-Ewald w geometric combination rules */
#define LJ_EWALD_COMB_GEOM
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecEwQSTab_VdwLJEwCombGeom ## __VA_ARGS__
#include "gromacs/mdlib/nbnxn_cuda/nbnxn_cuda_kernel.cuh"
#undef LJ_EWALD_COMB_GEOM
#undef NB_KERNEL_FUNC_NAME
/* LJ-Ewald w LB combination rules */
#define LJ_EWALD_COMB_LB
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecEwQSTab_VdwLJEwCombLB ## __VA_ARGS__
#include "gromacs/mdlib/nbnxn_cuda/nbnxn_cuda_kernel.cuh"
#undef LJ_EWALD_COMB_LB
#undef NB_KERNEL_FUNC_NAME
/* F switch LJ */
#define LJ_FORCE_SWITCH
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecEwQSTab_VdwLJFsw ## __VA_ARGS__
#include "gromacs/mdlib/nbnxn_cuda/nbnxn_cuda_kernel.cuh"
#undef LJ_FORCE_SWITCH
#undef NB_KERNEL_FUNC_NAME
/* V switch LJ */
#define LJ_POT_SWITCH
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecEwQSTab_VdwLJPsw ## __VA_ARGS__
#include "gromacs/mdlib/nbnxn_cuda/nbnxn_cuda_kernel.cuh"
#undef LJ_POT_SWITCH
#undef NB_KERNEL_FUNC_NAME
/* Vdw USER tables LJ6+LJ12 */
#define VDW_USER
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecEwQSTab_Vdw_USER ## __VA_ARGS__
#include "gromacs/mdlib/nbnxn_cuda/nbnxn_cuda_kernel.cuh"
#undef VDW_USER
#undef NB_KERNEL_FUNC_NAME
/* Vdw USER tables GENERIC */
#define VDW_USER_GENERIC
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecEwQSTab_Vdw_USER_GENERIC ## __VA_ARGS__
#include "gromacs/mdlib/nbnxn_cuda/nbnxn_cuda_kernel.cuh"
#undef VDW_USER_GENERIC
#undef NB_KERNEL_FUNC_NAME

#undef EL_EWALD_TAB


/* Tabulated Ewald interaction kernels with twin-range cut-off */
#define EL_EWALD_TAB
#define VDW_CUTOFF_CHECK

/* cut-off + V shift LJ */
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecEwQSTabTwinCut_VdwLJ ## __VA_ARGS__
#include "gromacs/mdlib/nbnxn_cuda/nbnxn_cuda_kernel.cuh"
#undef NB_KERNEL_FUNC_NAME
/* LJ-Ewald w geometric combination rules */
#define LJ_EWALD_COMB_GEOM
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecEwQSTabTwinCut_VdwLJEwCombGeom ## __VA_ARGS__
#include "gromacs/mdlib/nbnxn_cuda/nbnxn_cuda_kernel.cuh"
#undef LJ_EWALD_COMB_GEOM
#undef NB_KERNEL_FUNC_NAME
/* LJ-Ewald w LB combination rules */
#define LJ_EWALD_COMB_LB
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecEwQSTabTwinCut_VdwLJEwCombLB ## __VA_ARGS__
#include "gromacs/mdlib/nbnxn_cuda/nbnxn_cuda_kernel.cuh"
#undef LJ_EWALD_COMB_LB
#undef NB_KERNEL_FUNC_NAME
/* F switch LJ */
#define LJ_FORCE_SWITCH
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecEwQSTabTwinCut_VdwLJFsw ## __VA_ARGS__
#include "gromacs/mdlib/nbnxn_cuda/nbnxn_cuda_kernel.cuh"
#undef LJ_FORCE_SWITCH
#undef NB_KERNEL_FUNC_NAME
/* V switch LJ */
#define LJ_POT_SWITCH
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecEwQSTabTwinCut_VdwLJPsw ## __VA_ARGS__
#include "gromacs/mdlib/nbnxn_cuda/nbnxn_cuda_kernel.cuh"
#undef LJ_POT_SWITCH
#undef NB_KERNEL_FUNC_NAME
/* Vdw USER tables LJ6+LJ12 */
#define VDW_USER
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecEwQSTabTwinCut_Vdw_USER ## __VA_ARGS__
#include "gromacs/mdlib/nbnxn_cuda/nbnxn_cuda_kernel.cuh"
#undef VDW_USER
#undef NB_KERNEL_FUNC_NAME
/* Vdw USER tables GENERIC */
#define VDW_USER_GENERIC
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecEwQSTabTwinCut_Vdw_USER_GENERIC ## __VA_ARGS__
#include "gromacs/mdlib/nbnxn_cuda/nbnxn_cuda_kernel.cuh"
#undef VDW_USER_GENERIC
#undef NB_KERNEL_FUNC_NAME

#undef EL_EWALD_TAB
#undef VDW_CUTOFF_CHECK



/* Non-Bonded Tabulated pair potential kernels
 */

#define EL_NONE

#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecNone_VdwLJ ## __VA_ARGS__
#include "gromacs/mdlib/nbnxn_cuda/nbnxn_cuda_kernel.cuh"
#undef NB_KERNEL_FUNC_NAME
/* LJ-Ewald w geometric combination rules */
#define LJ_EWALD_COMB_GEOM
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecNone_VdwLJEwCombGeom ## __VA_ARGS__
#include "gromacs/mdlib/nbnxn_cuda/nbnxn_cuda_kernel.cuh"
#undef LJ_EWALD_COMB_GEOM
#undef NB_KERNEL_FUNC_NAME
/* LJ-Ewald w LB combination rules */
#define LJ_EWALD_COMB_LB
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecNone_VdwLJEwCombLB ## __VA_ARGS__
#include "gromacs/mdlib/nbnxn_cuda/nbnxn_cuda_kernel.cuh"
#undef LJ_EWALD_COMB_LB
#undef NB_KERNEL_FUNC_NAME
/* F switch LJ */
#define LJ_FORCE_SWITCH
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecNone_VdwLJFsw ## __VA_ARGS__
#include "gromacs/mdlib/nbnxn_cuda/nbnxn_cuda_kernel.cuh"
#undef LJ_FORCE_SWITCH
#undef NB_KERNEL_FUNC_NAME
/* V switch LJ */
#define LJ_POT_SWITCH
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecNone_VdwLJPsw ## __VA_ARGS__
#include "gromacs/mdlib/nbnxn_cuda/nbnxn_cuda_kernel.cuh"
#undef LJ_POT_SWITCH
#undef NB_KERNEL_FUNC_NAME
/* Vdw USER tables LJ6+LJ12 */
#define VDW_USER
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecNone_Vdw_USER ## __VA_ARGS__
#include "gromacs/mdlib/nbnxn_cuda/nbnxn_cuda_kernel.cuh"
#undef VDW_USER
#undef NB_KERNEL_FUNC_NAME
/* Vdw USER tables GENERIC */
#define VDW_USER_GENERIC
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecNone_Vdw_USER_GENERIC ## __VA_ARGS__
#include "gromacs/mdlib/nbnxn_cuda/nbnxn_cuda_kernel.cuh"
#undef VDW_USER_GENERIC
#undef NB_KERNEL_FUNC_NAME

#undef EL_NONE

#define EL_USER

#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecUSER_VdwLJ ## __VA_ARGS__
#include "gromacs/mdlib/nbnxn_cuda/nbnxn_cuda_kernel.cuh"
#undef NB_KERNEL_FUNC_NAME
/* LJ-Ewald w geometric combination rules */
#define LJ_EWALD_COMB_GEOM
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecUSER_VdwLJEwCombGeom ## __VA_ARGS__
#include "gromacs/mdlib/nbnxn_cuda/nbnxn_cuda_kernel.cuh"
#undef LJ_EWALD_COMB_GEOM
#undef NB_KERNEL_FUNC_NAME
/* LJ-Ewald w LB combination rules */
#define LJ_EWALD_COMB_LB
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecUSER_VdwLJEwCombLB ## __VA_ARGS__
#include "gromacs/mdlib/nbnxn_cuda/nbnxn_cuda_kernel.cuh"
#undef LJ_EWALD_COMB_LB
#undef NB_KERNEL_FUNC_NAME
/* F switch LJ */
#define LJ_FORCE_SWITCH
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecUSER_VdwLJFsw ## __VA_ARGS__
#include "gromacs/mdlib/nbnxn_cuda/nbnxn_cuda_kernel.cuh"
#undef LJ_FORCE_SWITCH
#undef NB_KERNEL_FUNC_NAME
/* V switch LJ */
#define LJ_POT_SWITCH
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecUSER_VdwLJPsw ## __VA_ARGS__
#include "gromacs/mdlib/nbnxn_cuda/nbnxn_cuda_kernel.cuh"
#undef LJ_POT_SWITCH
#undef NB_KERNEL_FUNC_NAME
/* Vdw USER tables LJ6+LJ12 */
#define VDW_USER
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecUSER_Vdw_USER ## __VA_ARGS__
#include "gromacs/mdlib/nbnxn_cuda/nbnxn_cuda_kernel.cuh"
#undef VDW_USER
#undef NB_KERNEL_FUNC_NAME
/* Vdw USER tables GENERIC */
#define VDW_USER_GENERIC
#define NB_KERNEL_FUNC_NAME(x, ...) x ## _ElecUSER_Vdw_USER_GENERIC ## __VA_ARGS__
#include "gromacs/mdlib/nbnxn_cuda/nbnxn_cuda_kernel.cuh"
#undef VDW_USER_GENERIC
#undef NB_KERNEL_FUNC_NAME

#undef EL_USER
