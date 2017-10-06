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

#include "nbnxn_kernel_simd_4xn_prune.h"

#include "gromacs/mdlib/nbnxn_pairlist.h"
#include "gromacs/mdlib/nbnxn_simd.h"
#include "gromacs/utility/gmxassert.h"

#ifdef GMX_NBNXN_SIMD_4XN
#define GMX_SIMD_J_UNROLL_SIZE 1
#include "gromacs/mdlib/nbnxn_kernels/simd_4xn/nbnxn_kernel_simd_4xn_common.h"
#endif

/* Prune a single nbnxn_pairtlist_t entry with distance rlistInner */
void
nbnxn_kernel_prune_4xn(nbnxn_pairlist_t *         nbl,
                       const nbnxn_atomdata_t *   nbat,
                       const rvec * gmx_restrict  shift_vec,
                       real                       rlistInner)
{
#ifdef GMX_NBNXN_SIMD_4XN
    const nbnxn_ci_t * gmx_restrict ciOuter  = nbl->ciOuter;
    nbnxn_ci_t       * gmx_restrict ciInner  = nbl->ci;

    const nbnxn_cj_t * gmx_restrict cjOuter  = nbl->cjOuter;
    nbnxn_cj_t       * gmx_restrict cjInner  = nbl->cj;

    const real       * gmx_restrict shiftvec = shift_vec[0];
    const real       * gmx_restrict x        = nbat->x;

    const SimdReal                  rlist2_S(rlistInner*rlistInner);

    /* Initialize the new list as empty and add pairs that are in range */
    int nciInner = 0;
    int ncjInner = 0;
    for (int i = 0; i < nbl->nciOuter; i++)
    {
        const nbnxn_ci_t * gmx_restrict ciEntry = &ciOuter[i];

        /* Copy the original list entry to the pruned entry */
        ciInner[nciInner].ci           = ciEntry->ci;
        ciInner[nciInner].shift        = ciEntry->shift;
        ciInner[nciInner].cj_ind_start = ncjInner;

        /* Extract shift data */
        int      ish     = (ciEntry->shift & NBNXN_CI_SHIFT);
        int      ish3    = ish*3;
        int      ci      = ciEntry->ci;

        SimdReal shX_S   = SimdReal(shiftvec[ish3    ]);
        SimdReal shY_S   = SimdReal(shiftvec[ish3 + 1]);
        SimdReal shZ_S   = SimdReal(shiftvec[ish3 + 2]);

#if UNROLLJ <= 4
        int      sci     = ci*STRIDE;
        int      scix    = sci*DIM;
#else
        int      sci     = (ci >> 1)*STRIDE;
        int      scix    = sci*DIM + (ci & 1)*(STRIDE >> 1);
        sci             += (ci & 1)*(STRIDE >> 1);
#endif

        /* Load i atom data */
        int      sciy    = scix + STRIDE;
        int      sciz    = sciy + STRIDE;
        SimdReal ix_S0   = SimdReal(x[scix    ]) + shX_S;
        SimdReal ix_S1   = SimdReal(x[scix + 1]) + shX_S;
        SimdReal ix_S2   = SimdReal(x[scix + 2]) + shX_S;
        SimdReal ix_S3   = SimdReal(x[scix + 3]) + shX_S;
        SimdReal iy_S0   = SimdReal(x[sciy    ]) + shY_S;
        SimdReal iy_S1   = SimdReal(x[sciy + 1]) + shY_S;
        SimdReal iy_S2   = SimdReal(x[sciy + 2]) + shY_S;
        SimdReal iy_S3   = SimdReal(x[sciy + 3]) + shY_S;
        SimdReal iz_S0   = SimdReal(x[sciz    ]) + shZ_S;
        SimdReal iz_S1   = SimdReal(x[sciz + 1]) + shZ_S;
        SimdReal iz_S2   = SimdReal(x[sciz + 2]) + shZ_S;
        SimdReal iz_S3   = SimdReal(x[sciz + 3]) + shZ_S;

        for (int cjind = ciEntry->cj_ind_start; cjind < ciEntry->cj_ind_end; cjind++)
        {
            /* j-cluster index */
            int cj      = cjOuter[cjind].cj;

            /* Atom indices (of the first atom in the cluster) */
            int aj      = cj*UNROLLJ;
#if UNROLLJ == STRIDE
            int ajx     = aj*DIM;
#else
            int ajx     = (cj >> 1)*DIM*STRIDE + (cj & 1)*UNROLLJ;
#endif
            int ajy     = ajx + STRIDE;
            int ajz     = ajy + STRIDE;

            /* load j atom coordinates */
            SimdReal jx_S   = load<SimdReal>(x + ajx);
            SimdReal jy_S   = load<SimdReal>(x + ajy);
            SimdReal jz_S   = load<SimdReal>(x + ajz);

            /* Calculate distance */
            SimdReal dx_S0  = ix_S0 - jx_S;
            SimdReal dy_S0  = iy_S0 - jy_S;
            SimdReal dz_S0  = iz_S0 - jz_S;
            SimdReal dx_S1  = ix_S1 - jx_S;
            SimdReal dy_S1  = iy_S1 - jy_S;
            SimdReal dz_S1  = iz_S1 - jz_S;
            SimdReal dx_S2  = ix_S2 - jx_S;
            SimdReal dy_S2  = iy_S2 - jy_S;
            SimdReal dz_S2  = iz_S2 - jz_S;
            SimdReal dx_S3  = ix_S3 - jx_S;
            SimdReal dy_S3  = iy_S3 - jy_S;
            SimdReal dz_S3  = iz_S3 - jz_S;

            /* rsq = dx*dx+dy*dy+dz*dz */
            SimdReal rsq_S0 = norm2(dx_S0, dy_S0, dz_S0);
            SimdReal rsq_S1 = norm2(dx_S1, dy_S1, dz_S1);
            SimdReal rsq_S2 = norm2(dx_S2, dy_S2, dz_S2);
            SimdReal rsq_S3 = norm2(dx_S3, dy_S3, dz_S3);

            /* Do the cut-off check */
            SimdBool wco_S0 = (rsq_S0 < rlist2_S);
            SimdBool wco_S1 = (rsq_S1 < rlist2_S);
            SimdBool wco_S2 = (rsq_S2 < rlist2_S);
            SimdBool wco_S3 = (rsq_S3 < rlist2_S);

            wco_S0          = wco_S0 || wco_S1;
            wco_S2          = wco_S2 || wco_S3;
            wco_S0          = wco_S0 || wco_S2;

            /* Putting the assignment inside the conditional is slower */
            cjInner[ncjInner] = cjOuter[cjind];
            if (anyTrue(wco_S0))
            {
                ncjInner++;
            }
        }

        if (ncjInner > ciInner[nciInner].cj_ind_start)
        {
            ciInner[nciInner].cj_ind_end = ncjInner;
            nciInner++;
        }
    }

    nbl->nci = nciInner;

#else  /* GMX_NBNXN_SIMD_4XN */

    GMX_RELEASE_ASSERT(false, "4xN kernel called without 4xN support");

    GMX_UNUSED_VALUE(nbl);
    GMX_UNUSED_VALUE(nbat);
    GMX_UNUSED_VALUE(shift_vec);
    GMX_UNUSED_VALUE(rlistInner);

#endif /* GMX_NBNXN_SIMD_4XN */
}
