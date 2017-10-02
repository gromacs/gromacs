/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2017, by the GROMACS development team, led by
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

#include "gmxpre.h"

#include "nbnxn_kernel_prune.h"

#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdlib/nb_verlet.h"
#include "gromacs/mdlib/nbnxn_pairlist.h"
#include "gromacs/utility/gmxassert.h"

#include "nbnxn_kernel_ref_prune.h"
#include "simd_2xnn/nbnxn_kernel_simd_2xnn_prune.h"
#include "simd_4xn/nbnxn_kernel_simd_4xn_prune.h"


void nbnxn_kernel_cpu_prune(nonbonded_verlet_group_t *nbvg,
                            const nbnxn_atomdata_t   *nbat,
                            const rvec               *shift_vec,
                            real                      rlistInner)
{
    nbnxn_pairlist_set_t   *nbl_lists = &nbvg->nbl_lists;

    GMX_ASSERT(nbl_lists->nbl[0]->nciOuter >= 0, "nciOuter<0, which signals an invalid pair-list");

    // cppcheck-suppress unreadVariable
    int gmx_unused nthreads = gmx_omp_nthreads_get(emntNonbonded);
#pragma omp parallel for schedule(static) num_threads(nthreads)
    for (int i = 0; i < nbl_lists->nnbl; i++)
    {
        nbnxn_pairlist_t *nbl = nbl_lists->nbl[i];

        switch (nbvg->kernel_type)
        {
            case nbnxnk4xN_SIMD_4xN:
                nbnxn_kernel_prune_4xn(nbl, nbat, shift_vec, rlistInner);
                break;
            case nbnxnk4xN_SIMD_2xNN:
                nbnxn_kernel_prune_2xnn(nbl, nbat, shift_vec, rlistInner);
                break;
            case nbnxnk4x4_PlainC:
                nbnxn_kernel_prune_ref(nbl, nbat, shift_vec, rlistInner);
                break;
            default:
                GMX_RELEASE_ASSERT(false, "kernel type not handled (yet)");
        }
    }
}
