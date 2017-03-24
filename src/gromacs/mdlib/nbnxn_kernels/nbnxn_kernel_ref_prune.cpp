/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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

#include "nbnxn_kernel_ref_prune.h"

#include "gromacs/mdlib/nbnxn_consts.h"
#include "gromacs/mdlib/nbnxn_pairlist.h"
#include "gromacs/utility/gmxassert.h"


/* Prune a single nbnxn_pairlist_t entry with distance rlistInner */
void
nbnxn_kernel_prune_ref(nbnxn_pairlist_t *         nbl,
                       const nbnxn_atomdata_t *   nbat,
                       const rvec * gmx_restrict  shift_vec,
                       real                       rlistInner)
{
    const nbnxn_ci_t * gmx_restrict ciOuter  = nbl->ciOuter;
    nbnxn_ci_t       * gmx_restrict ciInner  = nbl->ci;

    const nbnxn_cj_t * gmx_restrict cjOuter   = nbl->cjOuter;
    nbnxn_cj_t       * gmx_restrict cjInner   = nbl->cj;

    const real       * gmx_restrict shiftvec = shift_vec[0];
    const real       * gmx_restrict x        = nbat->x;

    const real                      rlist2   = rlistInner*rlistInner;

    /* Use compile time constants to speed up the code */
    constexpr int c_xStride  = 3;
    GMX_ASSERT(c_xStride == nbat->xstride, "xStride should match nbat->xstride");
    constexpr int c_xiStride = 3;

    constexpr int c_iUnroll  = NBNXN_CPU_CLUSTER_I_SIZE;
    constexpr int c_jUnroll  = NBNXN_CPU_CLUSTER_I_SIZE;

    /* Initialize the new list as empty and add pairs that are in range */
    int nciInner = 0;
    int ncjInner = 0;
    for (int ciIndex = 0; ciIndex < nbl->nciOuter; ciIndex++)
    {
        const nbnxn_ci_t * gmx_restrict ciEntry = &ciOuter[ciIndex];

        /* Copy the original list entry to the pruned entry */
        ciInner[nciInner].ci           = ciEntry->ci;
        ciInner[nciInner].shift        = ciEntry->shift;
        ciInner[nciInner].cj_ind_start = ncjInner;

        /* Extract shift data */
        int ish = (ciEntry->shift & NBNXN_CI_SHIFT);
        int ci  = ciEntry->ci;

        /* Load i atom coordinates */
        real xi[c_iUnroll*c_xiStride];
        for (int i = 0; i < c_iUnroll; i++)
        {
            for (int d = 0; d < DIM; d++)
            {
                xi[i*c_xiStride + d] = x[(ci*c_iUnroll + i)*c_xStride + d] + shiftvec[ish*DIM + d];
            }
        }

        for (int cjind = ciEntry->cj_ind_start; cjind < ciEntry->cj_ind_end; cjind++)
        {
            /* j-cluster index */
            int  cj        = cjOuter[cjind].cj;

            bool isInRange = false;
            for (int i = 0; i < c_iUnroll && !isInRange; i++)
            {
                for (int j = 0; j < c_jUnroll; j++)
                {
                    int  aj  = cj*c_jUnroll + j;

                    real dx  = xi[i*c_xiStride + XX] - x[aj*c_xStride + XX];
                    real dy  = xi[i*c_xiStride + YY] - x[aj*c_xStride + YY];
                    real dz  = xi[i*c_xiStride + ZZ] - x[aj*c_xStride + ZZ];

                    real rsq = dx*dx + dy*dy + dz*dz;

                    if (rsq < rlist2)
                    {
                        isInRange = true;
                    }
                }
            }

            if (isInRange)
            {
                /* This cluster is in range, put it in the pruned list */
                cjInner[ncjInner++] = cjOuter[cjind];
            }
        }

        /* Check if there are any j's in the list, if so, add the i-entry */
        if (ncjInner > ciInner[nciInner].cj_ind_start)
        {
            ciInner[nciInner].cj_ind_end = ncjInner;
            nciInner++;
        }
    }

    nbl->nci = nciInner;
}
