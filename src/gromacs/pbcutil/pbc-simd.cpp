/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015, by the GROMACS development team, led by
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
 * \brief This file defines a low-level function for SIMD PBC calculation.
 *
 * \author Berk Hess <hess@kth.se>
 *
 * \ingroup module_pbcutil
 */
#include "gmxpre.h"

#include "pbc-simd.h"

#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/simd/simd.h"

void set_pbc_simd(const t_pbc gmx_unused *pbc,
                  pbc_simd_t gmx_unused  *pbc_simd)
{
#ifdef GMX_SIMD_HAVE_REAL
    rvec inv_box_diag;
    int  d;

    /* Setting inv_bdiag to 0 effectively turns off PBC */
    clear_rvec(inv_box_diag);
    if (pbc != NULL)
    {
        for (d = 0; d < pbc->ndim_ePBC; d++)
        {
            inv_box_diag[d] = 1.0/pbc->box[d][d];
        }
    }

    pbc_simd->inv_bzz = gmx_simd_set1_r(inv_box_diag[ZZ]);
    pbc_simd->inv_byy = gmx_simd_set1_r(inv_box_diag[YY]);
    pbc_simd->inv_bxx = gmx_simd_set1_r(inv_box_diag[XX]);

    if (pbc != NULL)
    {
        pbc_simd->bzx = gmx_simd_set1_r(pbc->box[ZZ][XX]);
        pbc_simd->bzy = gmx_simd_set1_r(pbc->box[ZZ][YY]);
        pbc_simd->bzz = gmx_simd_set1_r(pbc->box[ZZ][ZZ]);
        pbc_simd->byx = gmx_simd_set1_r(pbc->box[YY][XX]);
        pbc_simd->byy = gmx_simd_set1_r(pbc->box[YY][YY]);
        pbc_simd->bxx = gmx_simd_set1_r(pbc->box[XX][XX]);
    }
    else
    {
        pbc_simd->bzx = gmx_simd_setzero_r();
        pbc_simd->bzy = gmx_simd_setzero_r();
        pbc_simd->bzz = gmx_simd_setzero_r();
        pbc_simd->byx = gmx_simd_setzero_r();
        pbc_simd->byy = gmx_simd_setzero_r();
        pbc_simd->bxx = gmx_simd_setzero_r();
    }
#endif
}
