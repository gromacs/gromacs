/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015,2016,2017, by the GROMACS development team, led by
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
/*! \libinternal \file
 *
 * \brief This file contains a definition, declaration and inline function
 * for SIMD accelerated PBC calculations.
 *
 * \author Berk Hess <hess@kth.se>
 * \inlibraryapi
 * \ingroup module_pbcutil
 */

#ifndef GMX_PBCUTIL_PBC_SIMD_H
#define GMX_PBCUTIL_PBC_SIMD_H

#include "config.h"

#include "gromacs/pbcutil/pbc.h"
#include "gromacs/simd/simd.h"

using namespace gmx; // TODO: Remove when this file is moved into gmx namespace

struct gmx_domdec_t;

/*! \brief Set the SIMD PBC data from a normal t_pbc struct.
 *
 * \param pbc        Type of periodic boundary,
 *                   NULL can be passed for then no PBC will be used.
 * \param pbc_simd   Pointer to aligned memory with (DIM + DIM*(DIM+1)/2)
 *                   GMX_SIMD_REAL_WIDTH reals describing the box vectors
 *                   unrolled by GMX_SIMD_REAL_WIDTH.
 *                   These are sorted in a slightly non-standard
 *                   order so that we always issue the memory loads in order
 *                   (to improve prefetching) in pbc_correct_dx_simd().
 *                   The order is inv_bzz, bzx, bzy, bzz, inv_byy, byx, byy,
 *                   inv_bxx, and bxx.
 */
void set_pbc_simd(const t_pbc *pbc,
                  real        *pbc_simd);

#if GMX_SIMD_HAVE_REAL

/*! \brief Correct SIMD distance vector *dx,*dy,*dz for PBC using SIMD.
 *
 * For rectangular boxes all returned distance vectors are the shortest.
 * For triclinic boxes only distances up to half the smallest box diagonal
 * element are guaranteed to be the shortest. This means that distances from
 * 0.5/sqrt(2) times a box vector length (e.g. for a rhombic dodecahedron)
 * can use a more distant periodic image.
 * Note that this routine always does PBC arithmetic, even for dimensions
 * without PBC. But on modern processors the overhead of this, often called,
 * routine should be low. On e.g. Intel Haswell/Broadwell it takes 8 cycles.
 */
static gmx_inline void gmx_simdcall
pbc_correct_dx_simd(SimdReal         *dx,
                    SimdReal         *dy,
                    SimdReal         *dz,
                    const real       *pbc_simd)
{
    SimdReal shz, shy, shx;

    shz = round(*dz * load<SimdReal>(pbc_simd+0*GMX_SIMD_REAL_WIDTH)); // load inv_bzz
    *dx = *dx - shz * load<SimdReal>(pbc_simd+1*GMX_SIMD_REAL_WIDTH);  // load bzx
    *dy = *dy - shz * load<SimdReal>(pbc_simd+2*GMX_SIMD_REAL_WIDTH);  // load bzy
    *dz = *dz - shz * load<SimdReal>(pbc_simd+3*GMX_SIMD_REAL_WIDTH);  // load bzz

    shy = round(*dy * load<SimdReal>(pbc_simd+4*GMX_SIMD_REAL_WIDTH)); // load inv_byy
    *dx = *dx - shy * load<SimdReal>(pbc_simd+5*GMX_SIMD_REAL_WIDTH);  // load byx
    *dy = *dy - shy * load<SimdReal>(pbc_simd+6*GMX_SIMD_REAL_WIDTH);  // load byy

    shx = round(*dx * load<SimdReal>(pbc_simd+7*GMX_SIMD_REAL_WIDTH)); // load inv_bxx
    *dx = *dx - shx * load<SimdReal>(pbc_simd+8*GMX_SIMD_REAL_WIDTH);  // load bxx

}

/*! \brief Calculates the PBC corrected distance between SIMD coordinates.
 *
 * \param pbc_simd  SIMD formatted PBC information
 * \param x1        Packed coordinates of atom1, size 3*GMX_SIMD_REAL_WIDTH
 * \param x2        Packed coordinates of atom2, size 3*GMX_SIMD_REAL_WIDTH
 * \param dx        The PBC corrected distance x1 - x2
 *
 * This routine only returns the shortest distance correctd for PBC
 * when all atoms are in the unit-cell (aiuc).
 * This is the SIMD equivalent of the scalar version declared in pbc.h.
 */
static gmx_inline void gmx_simdcall
pbc_dx_aiuc(const real       *pbc_simd,
            const SimdReal   *x1,
            const SimdReal   *x2,
            SimdReal         *dx)
{
    for (int d = 0; d < DIM; d++)
    {
        dx[d] = x1[d] - x2[d];
    }
    pbc_correct_dx_simd(&dx[XX], &dx[YY], &dx[ZZ], pbc_simd);
}

#endif /* GMX_SIMD_HAVE_REAL */

#endif
