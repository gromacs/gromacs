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

#include "kernel_ref_prune.h"

#include "../nbnxm_geometry.h"

#include "gromacs/nbnxm/atomdata.h"
#include "gromacs/nbnxm/pairlist.h"
#include "gromacs/utility/gmxassert.h"

namespace gmx
{

/* Prune a single NbnxnPairlistCpu entry with distance rlistInner */
void nbnxn_kernel_prune_ref(NbnxnPairlistCpu*       nbl,
                            const nbnxn_atomdata_t* nbat,
                            ArrayRef<const RVec>    shiftvec,
                            real                    rlistInner)
{
    /* We avoid push_back() for efficiency reasons and resize after filling */
    nbl->ci.resize(nbl->ciOuter.size());
    nbl->cj.resize(nbl->cjOuter.size());

    const nbnxn_ci_t* gmx_restrict ciOuter = nbl->ciOuter.data();
    nbnxn_ci_t* gmx_restrict       ciInner = nbl->ci.data();

    const nbnxn_cj_t* gmx_restrict cjOuter = nbl->cjOuter.data();
    nbnxn_cj_t* gmx_restrict       cjInner = nbl->cj.list_.data();

    const real* gmx_restrict x = nbat->x().data();

    const real rlist2 = rlistInner * rlistInner;

    /* Use compile time constants to speed up the code */
    constexpr int c_xStride = 3;
    GMX_ASSERT(c_xStride == nbat->xstride, "xStride should match nbat->xstride");
    constexpr int c_xiStride = 3;

    constexpr int c_iUnroll = sc_iClusterSize(NbnxmKernelType::Cpu4x4_PlainC);
    constexpr int c_jUnroll = sc_jClusterSize(NbnxmKernelType::Cpu4x4_PlainC);

    /* Initialize the new list as empty and add pairs that are in range */
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
        int ish = (ciEntry->shift & NBNXN_CI_SHIFT);
        int ci  = ciEntry->ci;

        /* Load i atom coordinates */
        real xi[c_iUnroll * c_xiStride];
        for (int i = 0; i < c_iUnroll; i++)
        {
            for (int d = 0; d < DIM; d++)
            {
                xi[i * c_xiStride + d] = x[(ci * c_iUnroll + i) * c_xStride + d] + shiftvec[ish][d];
            }
        }

        for (int cjind = ciEntry->cj_ind_start; cjind < ciEntry->cj_ind_end; cjind++)
        {
            /* j-cluster index */
            int cj = cjOuter[cjind].cj;

            bool isInRange = false;
            for (int i = 0; i < c_iUnroll && !isInRange; i++)
            {
                for (int j = 0; j < c_jUnroll; j++)
                {
                    int aj = cj * c_jUnroll + j;

                    real dx = xi[i * c_xiStride + XX] - x[aj * c_xStride + XX];
                    real dy = xi[i * c_xiStride + YY] - x[aj * c_xStride + YY];
                    real dz = xi[i * c_xiStride + ZZ] - x[aj * c_xStride + ZZ];

                    real rsq = dx * dx + dy * dy + dz * dz;

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

    nbl->ci.resize(nciInner);
    nbl->cj.resize(ncjInner);
}

} // namespace gmx
