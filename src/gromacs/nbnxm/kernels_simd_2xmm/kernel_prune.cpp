/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2012- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */

#include "gmxpre.h"

#include "kernel_prune.h"

#include "gromacs/nbnxm/atomdata.h"
#include "gromacs/nbnxm/nbnxm_simd.h"
#include "gromacs/nbnxm/pairlist.h"
#include "gromacs/simd/simd.h"
#include "gromacs/simd/vector_operations.h"
#include "gromacs/utility/gmxassert.h"

/* Prune a single nbnxn_pairtlist_t entry with distance rlistInner */
void nbnxn_kernel_prune_2xnn(NbnxnPairlistCpu*              nbl,
                             const nbnxn_atomdata_t*        nbat,
                             gmx::ArrayRef<const gmx::RVec> shiftvec,
                             real                           rlistInner)
{
#if GMX_SIMD && GMX_SIMD_HAVE_HSIMD_UTIL_REAL
    using namespace gmx;

    constexpr KernelLayout kernelLayout = KernelLayout::r2xMM;

    // The number of j-clusters stored in a SIMD register
    constexpr int c_numJClustersPerSimdRegister = (kernelLayout == KernelLayout::r2xMM ? 2 : 1);

    // The i-cluster size
    constexpr int c_iClusterSize = 4;
    // The j-cluster size
    constexpr int c_jClusterSize(GMX_SIMD_REAL_WIDTH / c_numJClustersPerSimdRegister);

    // The stride of all atom data arrays
    constexpr int c_stride = std::max(c_iClusterSize, c_jClusterSize);

    /* We avoid push_back() for efficiency reasons and resize after filling */
    nbl->ci.resize(nbl->ciOuter.size());
    nbl->cj.resize(nbl->cjOuter.size());

    const nbnxn_ci_t* gmx_restrict ciOuter = nbl->ciOuter.data();
    nbnxn_ci_t* gmx_restrict       ciInner = nbl->ci.data();

    const nbnxn_cj_t* gmx_restrict cjOuter = nbl->cjOuter.data();
    nbnxn_cj_t* gmx_restrict       cjInner = nbl->cj.list_.data();

    const real* gmx_restrict x = nbat->x().data();

    const SimdReal rlist2_S(rlistInner * rlistInner);

    /* Initialize the new list count as empty and add pairs that are in range */
    int       nciInner = 0;
    int       ncjInner = 0;
    const int nciOuter = nbl->ciOuter.size();
    for (int i = 0; i < nciOuter; i++)
    {
        const nbnxn_ci_t* gmx_restrict ciEntry = &ciOuter[i];

        /* Copy the original list entry to the pruned entry */
        ciInner[nciInner].ci           = ciEntry->ci;
        ciInner[nciInner].shift        = ciEntry->shift;
        ciInner[nciInner].cj_ind_start = ncjInner;

        /* Extract shift data */
        int ish = (ciEntry->shift & NBNXN_CI_SHIFT);
        int ci  = ciEntry->ci;

        SimdReal shX_S = SimdReal(shiftvec[ish][XX]);
        SimdReal shY_S = SimdReal(shiftvec[ish][YY]);
        SimdReal shZ_S = SimdReal(shiftvec[ish][ZZ]);

        int scix;
        if constexpr (c_jClusterSize <= 4)
        {
            scix = ci * c_stride * DIM;
        }
        else
        {
            scix = (ci >> 1) * c_stride * DIM + (ci & 1) * (c_stride >> 1);
        }

        /* Load i atom data */
        int      sciy  = scix + c_stride;
        int      sciz  = sciy + c_stride;
        SimdReal ix_S0 = loadU1DualHsimd(x + scix) + shX_S;
        SimdReal ix_S2 = loadU1DualHsimd(x + scix + 2) + shX_S;
        SimdReal iy_S0 = loadU1DualHsimd(x + sciy) + shY_S;
        SimdReal iy_S2 = loadU1DualHsimd(x + sciy + 2) + shY_S;
        SimdReal iz_S0 = loadU1DualHsimd(x + sciz) + shZ_S;
        SimdReal iz_S2 = loadU1DualHsimd(x + sciz + 2) + shZ_S;

        for (int cjind = ciEntry->cj_ind_start; cjind < ciEntry->cj_ind_end; cjind++)
        {
            /* j-cluster index */
            int cj = cjOuter[cjind].cj;

            /* Atom indices (of the first atom in the cluster) */
            int ajx;
            if constexpr (c_jClusterSize == c_stride)
            {
                int aj = cj * c_jClusterSize;
                ajx    = aj * DIM;
            }
            else
            {
                ajx = (cj >> 1) * DIM * c_stride + (cj & 1) * c_jClusterSize;
            }
            int ajy = ajx + c_stride;
            int ajz = ajy + c_stride;

            /* load j atom coordinates */
            SimdReal jx_S = loadDuplicateHsimd(x + ajx);
            SimdReal jy_S = loadDuplicateHsimd(x + ajy);
            SimdReal jz_S = loadDuplicateHsimd(x + ajz);

            /* Calculate distance */
            SimdReal dx_S0 = ix_S0 - jx_S;
            SimdReal dy_S0 = iy_S0 - jy_S;
            SimdReal dz_S0 = iz_S0 - jz_S;
            SimdReal dx_S2 = ix_S2 - jx_S;
            SimdReal dy_S2 = iy_S2 - jy_S;
            SimdReal dz_S2 = iz_S2 - jz_S;

            /* rsq = dx*dx+dy*dy+dz*dz */
            SimdReal rsq_S0 = norm2(dx_S0, dy_S0, dz_S0);
            SimdReal rsq_S2 = norm2(dx_S2, dy_S2, dz_S2);

            /* Do the cut-off check */
            SimdBool wco_S0 = (rsq_S0 < rlist2_S);
            SimdBool wco_S2 = (rsq_S2 < rlist2_S);

            wco_S0 = wco_S0 || wco_S2;

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

    nbl->ci.resize(nciInner);
    nbl->cj.resize(ncjInner);

#else /* GMX_SIMD && GMX_SIMD_HAVE_HSIMD_UTIL_REAL */

    GMX_RELEASE_ASSERT(false, "2xNN kernel called without 2xNN support");

    GMX_UNUSED_VALUE(nbl);
    GMX_UNUSED_VALUE(nbat);
    GMX_UNUSED_VALUE(shiftvec);
    GMX_UNUSED_VALUE(rlistInner);

#endif /* GMX_SIMD && GMX_SIMD_HAVE_HSIMD_UTIL_REAL */
}
