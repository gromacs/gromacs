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
#include "gromacs/utility/fatalerror.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! \cond INTERNAL */

/*! \brief Structure containing the PBC setup for SIMD PBC calculations.
 *
 * Without SIMD this is a dummy struct, so it can be declared and passed.
 * This can avoid some ifdef'ing.
 */
typedef struct {
#ifdef GMX_SIMD_HAVE_REAL
    gmx_simd_real_t inv_bzz; /**< 1/box[ZZ][ZZ] */
    gmx_simd_real_t inv_byy; /**< 1/box[YY][YY] */
    gmx_simd_real_t inv_bxx; /**< 1/box[XX][XX] */
    gmx_simd_real_t bzx;     /**< box[ZZ][XX] */
    gmx_simd_real_t bzy;     /**< box[ZZ][YY] */
    gmx_simd_real_t bzz;     /**< box[ZZ][ZZ] */
    gmx_simd_real_t byx;     /**< box[YY][XX] */
    gmx_simd_real_t byy;     /**< box[YY][YY] */
    gmx_simd_real_t bxx;     /**< bo[XX][XX] */
#else
    int             dum;     /**< Dummy variable to avoid empty struct */
#endif
} pbc_simd_t;

/*! \endcond */

/*! \brief Set the SIMD PBC data from a normal t_pbc struct.
 *
 * NULL can be passed for \p pbc, then no PBC will be used.
 */
void set_pbc_simd(const t_pbc *pbc,
                  pbc_simd_t  *pbc_simd);

#if defined GMX_SIMD_HAVE_REAL

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
pbc_correct_dx_simd(gmx_simd_real_t  *dx,
                    gmx_simd_real_t  *dy,
                    gmx_simd_real_t  *dz,
                    const pbc_simd_t *pbc)
{
    gmx_simd_real_t shz, shy, shx;

#if defined _MSC_VER && _MSC_VER < 1700 && !defined(__ICL)
    /* The caller side should make sure we never end up here.
     * TODO Black-list _MSC_VER < 1700 when it's old enough, so we can rid
     * of this code complication.
     */
    gmx_incons("pbc_correct_dx_simd was called for code compiled with MSVC 2010 or older, which produces incorrect code (probably corrupts memory) and therefore this function should not have been called");
#endif

    shz = gmx_simd_round_r(gmx_simd_mul_r(*dz, pbc->inv_bzz));
    *dx = gmx_simd_fnmadd_r(shz, pbc->bzx, *dx);
    *dy = gmx_simd_fnmadd_r(shz, pbc->bzy, *dy);
    *dz = gmx_simd_fnmadd_r(shz, pbc->bzz, *dz);

    shy = gmx_simd_round_r(gmx_simd_mul_r(*dy, pbc->inv_byy));
    *dx = gmx_simd_fnmadd_r(shy, pbc->byx, *dx);
    *dy = gmx_simd_fnmadd_r(shy, pbc->byy, *dy);

    shx = gmx_simd_round_r(gmx_simd_mul_r(*dx, pbc->inv_bxx));
    *dx = gmx_simd_fnmadd_r(shx, pbc->bxx, *dx);
}

#endif /* GMX_SIMD_HAVE_REAL */

#ifdef __cplusplus
}
#endif

#endif
