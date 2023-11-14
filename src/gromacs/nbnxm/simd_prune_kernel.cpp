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

#include "simd_prune_kernel.h"

#include "gromacs/nbnxm/atomdata.h"
#include "gromacs/nbnxm/nbnxm_simd.h"
#include "gromacs/nbnxm/pairlist.h"
#include "gromacs/simd/simd.h"
#include "gromacs/simd/vector_operations.h"
#include "gromacs/utility/gmxassert.h"

#if GMX_SIMD
#    include "simd_load_store_functions.h"
#endif

//! Returns the base 2 log of \p x
template<int x>
gmx_unused static constexpr int sc_log2()
{
    static_assert(x == 2 || x == 4);

    switch (x)
    {
        case 2: return 1; break;
        case 4: return 2; break;
    }

    return 0;
}

template<KernelLayout kernelLayout>
void nbnxmSimdPruneKernel(NbnxnPairlistCpu*              nbl,
                          const nbnxn_atomdata_t&        nbat,
                          gmx::ArrayRef<const gmx::RVec> shiftvec,
                          real                           rlistInner)
{
#if GMX_SIMD
    using namespace gmx;

    // The number of j-clusters stored in a SIMD register
    constexpr int c_numJClustersPerSimdRegister = (kernelLayout == KernelLayout::r2xMM ? 2 : 1);

    // The i-cluster size
    constexpr int c_iClusterSize = 4;
    // The j-cluster size
    constexpr int c_jClusterSize(GMX_SIMD_REAL_WIDTH / c_numJClustersPerSimdRegister);

    // The stride of all atom data arrays
    constexpr int c_stride = std::max(c_iClusterSize, c_jClusterSize);

    // The number of 'i' SIMD registers
    constexpr int nR = c_iClusterSize / c_numJClustersPerSimdRegister;

    /* We avoid push_back() for efficiency reasons and resize after filling */
    nbl->ci.resize(nbl->ciOuter.size());
    nbl->cj.resize(nbl->cjOuter.size());

    const nbnxn_ci_t* gmx_restrict ciOuter = nbl->ciOuter.data();
    nbnxn_ci_t* gmx_restrict       ciInner = nbl->ci.data();

    const nbnxn_cj_t* gmx_restrict cjOuter = nbl->cjOuter.data();
    nbnxn_cj_t* gmx_restrict       cjInner = nbl->cj.list_.data();

    const real* gmx_restrict x = nbat.x().data();

    const SimdReal rlist2_S(rlistInner * rlistInner);

    /* Initialize the new list count as empty and add pairs that are in range */
    int       nciInner = 0;
    int       ncjInner = 0;
    const int nciOuter = nbl->ciOuter.size();
    for (int ciIndex = 0; ciIndex < nciOuter; ciIndex++)
    {
        const nbnxn_ci_t* gmx_restrict ciEntry = &ciOuter[ciIndex];

        /* Copy the original list entry to the pruned entry */
        ciInner[nciInner].ci           = ciEntry->ci;
        ciInner[nciInner].shift        = ciEntry->shift;
        ciInner[nciInner].cj_ind_start = ncjInner;

        /* Extract shift data */
        const int ish = (ciEntry->shift & NBNXN_CI_SHIFT);
        const int ci  = ciEntry->ci;

        const SimdReal iShiftX = SimdReal(shiftvec[ish][XX]);
        const SimdReal iShiftY = SimdReal(shiftvec[ish][YY]);
        const SimdReal iShiftZ = SimdReal(shiftvec[ish][ZZ]);

        int scix;
        if constexpr (c_jClusterSize <= c_iClusterSize)
        {
            scix = ci * c_stride * DIM;
        }
        else
        {
            scix = (ci >> 1) * c_stride * DIM + (ci & 1) * (c_stride >> 1);
        }

        /* Load i atom data */
        int                                       sciy = scix + c_stride;
        int                                       sciz = sciy + c_stride;
        std::array<std::array<SimdReal, DIM>, nR> xi;
        for (int i = 0; i < nR; i++)
        {
            xi[i][0] = loadIAtomData<kernelLayout>(x, scix, i) + iShiftX;
            xi[i][1] = loadIAtomData<kernelLayout>(x, sciy, i) + iShiftY;
            xi[i][2] = loadIAtomData<kernelLayout>(x, sciz, i) + iShiftZ;
        }

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
            const SimdReal jx_S = loadJAtomData<kernelLayout>(x, ajx);
            const SimdReal jy_S = loadJAtomData<kernelLayout>(x, ajy);
            const SimdReal jz_S = loadJAtomData<kernelLayout>(x, ajz);

            /* Calculate distance */
            std::array<std::array<SimdReal, DIM>, nR> d;
            for (int i = 0; i < nR; i++)
            {
                d[i][0] = xi[i][0] - jx_S;
                d[i][1] = xi[i][1] - jy_S;
                d[i][2] = xi[i][2] - jz_S;
            }

            /* rsq = dx*dx+dy*dy+dz*dz */
            std::array<SimdReal, nR> rsq;
            for (int i = 0; i < nR; i++)
            {
                rsq[i] = norm2(d[i][0], d[i][1], d[i][2]);
            }

            /* Do the cut-off check */
            std::array<SimdBool, nR> wco;
            for (int i = 0; i < nR; i++)
            {
                wco[i] = (rsq[i] < rlist2_S);
            }

            constexpr int numIterations = sc_log2<nR>();
            for (int iter = 0; iter < numIterations; iter++)
            {
                const int offset = (1 << iter);

                for (int i = 0; i < nR; i += 2 * offset)
                {
                    wco[i] = wco[i] || wco[i + offset];
                }
            }

            /* Putting the assignment inside the conditional is slower */
            cjInner[ncjInner] = cjOuter[cjind];
            if (anyTrue(wco[0]))
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

#else /* GMX_SIMD */

    GMX_RELEASE_ASSERT(false, "SIMD prune kernel called without SIMD support");

    GMX_UNUSED_VALUE(nbl);
    GMX_UNUSED_VALUE(nbat);
    GMX_UNUSED_VALUE(shiftvec);
    GMX_UNUSED_VALUE(rlistInner);

#endif /* GMX_SIMD */
}

#if GMX_HAVE_NBNXM_SIMD_2XMM
template void nbnxmSimdPruneKernel<KernelLayout::r2xMM>(NbnxnPairlistCpu*              nbl,
                                                        const nbnxn_atomdata_t&        nbat,
                                                        gmx::ArrayRef<const gmx::RVec> shiftvec,
                                                        real                           rlistInner);
#endif

#if GMX_HAVE_NBNXM_SIMD_4XM
template void nbnxmSimdPruneKernel<KernelLayout::r4xM>(NbnxnPairlistCpu*              nbl,
                                                       const nbnxn_atomdata_t&        nbat,
                                                       gmx::ArrayRef<const gmx::RVec> shiftvec,
                                                       real                           rlistInner);
#endif
