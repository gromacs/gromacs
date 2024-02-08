/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
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

#include "pme_spline_work.h"

#include "gromacs/simd/simd.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

#include "pme_simd.h"

using namespace gmx; // TODO: Remove when this file is moved into gmx namespace

pme_spline_work::pme_spline_work(int gmx_unused order)
{
#ifdef PME_SIMD4_SPREAD_GATHER
    alignas(GMX_SIMD_ALIGNMENT) real tmp[GMX_SIMD4_WIDTH * 2];
    Simd4Real                        zero_S;
    Simd4Real                        real_mask_S0, real_mask_S1;
    int                              of, i;

    zero_S = setZero();

    /* Generate bit masks to mask out the unused grid entries,
     * as we only operate on order of the 8 grid entries that are
     * load into 2 SIMD registers.
     */
    for (of = 0; of < 2 * GMX_SIMD4_WIDTH - (order - 1); of++)
    {
        for (i = 0; i < 2 * GMX_SIMD4_WIDTH; i++)
        {
            tmp[i] = (i >= of && i < of + order ? -1.0 : 1.0);
        }
        real_mask_S0 = load4(tmp);
        real_mask_S1 = load4(tmp + GMX_SIMD4_WIDTH);
        mask_S0[of]  = (real_mask_S0 < zero_S);
        mask_S1[of]  = (real_mask_S1 < zero_S);
    }
#endif
}
