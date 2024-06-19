/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2015- The GROMACS Authors
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
/*! \internal \file
 *
 * \brief This file defines a low-level function for SIMD PBC calculation.
 *
 * \author Berk Hess <hess@kth.se>
 *
 * \ingroup module_pbcutil
 */
#include "gmxpre.h"

#include "gromacs/pbcutil/pbc_simd.h"

#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/simd/simd.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

using namespace gmx; // TODO: Remove when this file is moved into gmx namespace

void set_pbc_simd(const t_pbc gmx_unused* pbc, real gmx_unused* pbc_simd)
{
#if GMX_SIMD_HAVE_REAL
    if (pbc != nullptr && pbc->pbcType != PbcType::No)
    {
        rvec inv_box_diag = { 0, 0, 0 };

        for (int d = 0; d < pbc->ndim_ePBC; d++)
        {
            inv_box_diag[d] = 1.0 / pbc->box[d][d];
        }

        store(pbc_simd + 0 * GMX_SIMD_REAL_WIDTH, SimdReal(inv_box_diag[ZZ]));
        store(pbc_simd + 1 * GMX_SIMD_REAL_WIDTH, SimdReal(pbc->box[ZZ][XX]));
        store(pbc_simd + 2 * GMX_SIMD_REAL_WIDTH, SimdReal(pbc->box[ZZ][YY]));
        store(pbc_simd + 3 * GMX_SIMD_REAL_WIDTH, SimdReal(pbc->box[ZZ][ZZ]));
        store(pbc_simd + 4 * GMX_SIMD_REAL_WIDTH, SimdReal(inv_box_diag[YY]));
        store(pbc_simd + 5 * GMX_SIMD_REAL_WIDTH, SimdReal(pbc->box[YY][XX]));
        store(pbc_simd + 6 * GMX_SIMD_REAL_WIDTH, SimdReal(pbc->box[YY][YY]));
        store(pbc_simd + 7 * GMX_SIMD_REAL_WIDTH, SimdReal(inv_box_diag[XX]));
        store(pbc_simd + 8 * GMX_SIMD_REAL_WIDTH, SimdReal(pbc->box[XX][XX]));
    }
    else
    {
        /* Setting inv_box_diag to zero leads to no PBC being applied */
        for (int i = 0; i < (DIM + DIM * (DIM + 1) / 2); i++)
        {
            store(pbc_simd + i * GMX_SIMD_REAL_WIDTH, SimdReal(0));
        }
    }
#endif
}
