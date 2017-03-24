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

#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdlib/nbnxn_consts.h"
#include "gromacs/mdlib/nbnxn_pairlist.h"
#include "gromacs/mdlib/nbnxn_simd.h"
#include "gromacs/utility/gmxassert.h"


#define UNROLLI    NBNXN_CPU_CLUSTER_I_SIZE
#define UNROLLJ    NBNXN_CPU_CLUSTER_I_SIZE

/* We could use nbat->xstride, but macros are faster */
#define X_STRIDE   3
/* Local i-atom buffer stride */
#define XI_STRIDE  3


/* Prune a single nbnxn_pairlist_t entry with distance rlistPrune */
void
nbnxn_kernel_prune_ref(nbnxn_pairlist_t        *nbl,
                       const nbnxn_atomdata_t  *nbat,
                       const rvec              *shift_vec,
                       real                     rlistPrune)
{
    const nbnxn_ci_t * gmx_restrict ciOrig   = nbl->ci0;
    nbnxn_ci_t       * gmx_restrict ciPruned = nbl->ci;

    const nbnxn_cj_t * gmx_restrict cjOrig   = nbl->cj0;
    nbnxn_cj_t       * gmx_restrict cjPruned = nbl->cj;

    const real       * gmx_restrict shiftvec = shift_vec[0];
    const real       * gmx_restrict x        = nbat->x;

    const real                      rlist2   = rlistPrune*rlistPrune;

    /* Initialize the new list as empty and add pairs that are in range */
    int nciPruned = 0;
    int ncjPruned = 0;
    for (int ciIndex = 0; ciIndex < nbl->nci0; ciIndex++)
    {
        const nbnxn_ci_t * gmx_restrict ciEntry = &ciOrig[ciIndex];

        /* Copy the original list entry to the pruned entry */
        ciPruned[nciPruned].ci           = ciEntry->ci;
        ciPruned[nciPruned].shift        = ciEntry->shift;
        ciPruned[nciPruned].cj_ind_start = ncjPruned;

        /* Extract shift data */
        int ish = (ciEntry->shift & NBNXN_CI_SHIFT);
        int ci  = ciEntry->ci;

        /* Load i atom coordinates */
        real xi[UNROLLI*XI_STRIDE];
        for (int i = 0; i < UNROLLI; i++)
        {
            for (int d = 0; d < DIM; d++)
            {
                xi[i*XI_STRIDE + d] = x[(ci*UNROLLI + i)*X_STRIDE + d] + shiftvec[ish*DIM + d];
            }
        }

        for (int cjind = ciEntry->cj_ind_start; cjind < ciEntry->cj_ind_end; cjind++)
        {
            /* j-cluster index */
            int  cj        = cjOrig[cjind].cj;

            bool isInRange = false;
            for (int i = 0; i < UNROLLI && !isInRange; i++)
            {
                for (int j = 0; j < UNROLLJ; j++)
                {
                    int  aj  = cj*UNROLLJ + j;

                    real dx  = xi[i*XI_STRIDE + XX] - x[aj*X_STRIDE + XX];
                    real dy  = xi[i*XI_STRIDE + YY] - x[aj*X_STRIDE + YY];
                    real dz  = xi[i*XI_STRIDE + ZZ] - x[aj*X_STRIDE + ZZ];

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
                cjPruned[ncjPruned++] = cjOrig[cjind];
            }
        }

        /* Check if there are any j's in the list, if so, add the i-entry */
        if (ncjPruned > ciPruned[nciPruned].cj_ind_start)
        {
            ciPruned[nciPruned].cj_ind_end = ncjPruned;
            nciPruned++;
        }
    }

    nbl->nci = nciPruned;
}
