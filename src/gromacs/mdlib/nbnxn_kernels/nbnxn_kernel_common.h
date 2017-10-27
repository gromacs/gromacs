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

/*! \internal \file
 *
 * \brief
 * Declares the nbnxn pair interaction kernel function types and kind counts, also declares utility functions used in nbnxn_kernel.cpp.
 *
 * \author Berk Hess <hess@kth.se>
 */

#ifndef _nbnxn_kernel_common_h
#define _nbnxn_kernel_common_h

#include "gromacs/fda/FDA.h"
#include "gromacs/math/vectypes.h"
/* nbnxn_atomdata_t and nbnxn_pairlist_t could be forward declared, but that requires modifications in all SIMD kernel files */
#include "gromacs/mdlib/nbnxn_atomdata.h"
#include "gromacs/mdlib/nbnxn_pairlist.h"
#include "gromacs/utility/real.h"

struct interaction_const_t;

/*! \brief Pair-interaction kernel type that also calculates energies.
 */
typedef void (nbk_func_ener)(const nbnxn_pairlist_t     *nbl,
                             const nbnxn_atomdata_t     *nbat,
                             const interaction_const_t  *ic,
                             rvec                       *shift_vec,
                             real                       *f,
                             real                       *fshift,
                             real                       *Vvdw,
                             real                       *Vc
#ifdef BUILD_WITH_FDA
							 ,
                             FDA                        *fda,
                             int                        *cellInv
#endif
                            );

/*! \brief Pointer to \p nbk_func_ener.
 */
typedef nbk_func_ener *p_nbk_func_ener;

/*! \brief Pair-interaction kernel type that does not calculates energies.
 */
typedef void (nbk_func_noener)(const nbnxn_pairlist_t     *nbl,
                               const nbnxn_atomdata_t     *nbat,
                               const interaction_const_t  *ic,
                               rvec                       *shift_vec,
                               real                       *f,
                               real                       *fshift);

/*! \brief Pointer to \p nbk_func_noener.
 */
typedef nbk_func_noener *p_nbk_func_noener;

/*! \brief Kinds of electrostatic treatments in SIMD Verlet kernels
 */
enum {
    coulktRF, coulktTAB, coulktTAB_TWIN, coulktEWALD, coulktEWALD_TWIN, coulktNR
};

/*! \brief Kinds of Van der Waals treatments in SIMD Verlet kernels
 *
 * The \p LJCUT_COMB refers to the LJ combination rule for the short range.
 * The \p EWALDCOMB refers to the combination rule for the grid part.
 * \p vdwktNR is the number of VdW treatments for the SIMD kernels.
 * \p vdwktNR_ref is the number of VdW treatments for the C reference kernels.
 * These two numbers differ, because currently only the reference kernels
 * support LB combination rules for the LJ-Ewald grid part.
 */
enum {
    vdwktLJCUT_COMBGEOM, vdwktLJCUT_COMBLB, vdwktLJCUT_COMBNONE, vdwktLJFORCESWITCH, vdwktLJPOTSWITCH, vdwktLJEWALDCOMBGEOM, vdwktLJEWALDCOMBLB, vdwktNR = vdwktLJEWALDCOMBLB, vdwktNR_ref
};

/*! \brief Clears the force buffer.
 *
 * Either the whole buffer is cleared or only the parts used
 * by the current thread when nbat->bUseBufferFlags is set.
 * In the latter case output_index is the task/thread list/buffer index.
 */
void
clear_f(const nbnxn_atomdata_t *nbat, int output_index, real *f);

/*! \brief Clears the shift forces.
 */
void
clear_fshift(real *fshift);

/*! \brief Reduces the collected energy terms over the pair-lists/threads.
 */
void
reduce_energies_over_lists(const nbnxn_atomdata_t     *nbat,
                           int                         nlist,
                           real                       *Vvdw,
                           real                       *Vc);

#endif
