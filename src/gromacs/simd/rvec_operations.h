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
 * \brief Load and store of Gromacs rvec (distances) to and from SIMD registers
 *
 * \author Berk Hess <hess@kth.se>
 *
 * \inlibraryapi
 * \ingroup module_simd
 */

#ifndef GMX_SIMD_RVEC_OPERATIONS_H
#define GMX_SIMD_RVEC_OPERATIONS_H

#include "config.h"

#include "gromacs/math/vec.h"
#include "gromacs/simd/simd.h"

/*! \cond libapi */
/*! \addtogroup module_simd */
/*! \{ */

#ifdef GMX_SIMD_HAVE_REAL

/*! \brief Store differences between indexed rvecs in SIMD registers.
 *
 * Returns SIMD register with the difference vectors:
 *     v[pair_index[i*2]] - v[pair_index[i*2 + 1]]
 *
 * \param[in]     v           Array of rvecs
 * \param[in]     pair_index  Index pairs for GMX_SIMD_REAL_WIDTH vector pairs
 * \param[in,out] buf_aligned Aligned tmp buffer of size 3*GMX_SIMD_REAL_WIDTH
 * \param[out]    dx          SIMD register with x difference
 * \param[out]    dy          SIMD register with y difference
 * \param[out]    dz          SIMD register with z difference
 */
static gmx_inline void gmx_simdcall
gmx_simd_gather_rvec_dist_pair_index(const rvec      *v,
                                     const int       *pair_index,
                                     real gmx_unused *buf_aligned,
                                     gmx_simd_real_t *dx,
                                     gmx_simd_real_t *dy,
                                     gmx_simd_real_t *dz)
{
#if defined GMX_SIMD4_HAVE_REAL && defined GMX_SIMD4_HAVE_MASKLOAD3 && defined GMX_SIMD4_HAVE_SIMD_TRANSPOSE
    int              i;
    gmx_simd4_real_t d[GMX_SIMD_REAL_WIDTH];
    gmx_simd_real_t  tmp;

    for (i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
    {
        d[i] = gmx_simd4_sub_r(gmx_simd4_maskload3_r(&(v[pair_index[i*2 + 0]][0])),
                               gmx_simd4_maskload3_r(&(v[pair_index[i*2 + 1]][0])));
    }

    gmx_simd4_transpose_to_simd_r(d, dx, dy, dz, &tmp);
#else
    int i, m;

    for (i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
    {
        /* Store the distances packed and aligned */
        for (m = 0; m < DIM; m++)
        {
            buf_aligned[m*GMX_SIMD_REAL_WIDTH + i] =
                v[pair_index[i*2]][m] - v[pair_index[i*2 + 1]][m];
        }
    }
    *dx = gmx_simd_load_r(buf_aligned + 0*GMX_SIMD_REAL_WIDTH);
    *dy = gmx_simd_load_r(buf_aligned + 1*GMX_SIMD_REAL_WIDTH);
    *dz = gmx_simd_load_r(buf_aligned + 2*GMX_SIMD_REAL_WIDTH);
#endif
}

/*! \brief Store differences between indexed rvecs in SIMD registers.
 *
 * Returns SIMD register with the difference vectors:
 *     v[index0[i]] - v[index1[i]]
 *
 * \param[in]     v           Array of rvecs
 * \param[in]     index0      Index into the vector array
 * \param[in]     index1      Index into the vector array
 * \param[in,out] buf_aligned Aligned tmp buffer of size 3*GMX_SIMD_REAL_WIDTH
 * \param[out]    dx          SIMD register with x difference
 * \param[out]    dy          SIMD register with y difference
 * \param[out]    dz          SIMD register with z difference
 */
static gmx_inline void gmx_simdcall
gmx_simd_gather_rvec_dist_two_index(const rvec      *v,
                                    const int       *index0,
                                    const int       *index1,
                                    real gmx_unused *buf_aligned,
                                    gmx_simd_real_t *dx,
                                    gmx_simd_real_t *dy,
                                    gmx_simd_real_t *dz)
{
#if defined GMX_SIMD4_HAVE_REAL && defined GMX_SIMD4_HAVE_MASKLOAD3 && defined GMX_SIMD4_HAVE_SIMD_TRANSPOSE
    int              i;
    gmx_simd4_real_t d[GMX_SIMD_REAL_WIDTH];
    gmx_simd_real_t  tmp;

    for (i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
    {
        d[i] = gmx_simd4_sub_r(gmx_simd4_maskload3_r(&(v[index0[i]][0])),
                               gmx_simd4_maskload3_r(&(v[index1[i]][0])));
    }

    gmx_simd4_transpose_to_simd_r(d, dx, dy, dz, &tmp);
#else
    int i, m;

    for (i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
    {
        /* Store the distances packed and aligned */
        for (m = 0; m < DIM; m++)
        {
            buf_aligned[m*GMX_SIMD_REAL_WIDTH + i] =
                v[index0[i]][m] - v[index1[i]][m];
        }
    }
    *dx = gmx_simd_load_r(buf_aligned + 0*GMX_SIMD_REAL_WIDTH);
    *dy = gmx_simd_load_r(buf_aligned + 1*GMX_SIMD_REAL_WIDTH);
    *dz = gmx_simd_load_r(buf_aligned + 2*GMX_SIMD_REAL_WIDTH);
#endif
}

/*! \brief Stores SIMD vector into multiple rvecs.
 *
 * \param[in]     x           SIMD register with x-components of the vectors
 * \param[in]     y           SIMD register with y-components of the vectors
 * \param[in]     z           SIMD register with z-components of the vectors
 * \param[in,out] buf_aligned Aligned tmp buffer of size 3*GMX_SIMD_REAL_WIDTH
 * \param[out]    v           Array of GMX_SIMD_REAL_WIDTH rvecs
 */
static gmx_inline void gmx_simdcall
gmx_simd_store_vec_to_rvec(gmx_simd_real_t  x,
                           gmx_simd_real_t  y,
                           gmx_simd_real_t  z,
                           real gmx_unused *buf_aligned,
                           rvec            *v)
{
#if defined GMX_SIMD4_HAVE_REAL && defined GMX_SIMD4_HAVE_MASKSTORE3 && defined GMX_SIMD4_HAVE_SIMD_TRANSPOSE
    int              i;
    gmx_simd4_real_t s4[GMX_SIMD_REAL_WIDTH];
    gmx_simd_real_t  zero = gmx_simd_setzero_r();

    gmx_simd_transpose_to_simd4_r(x, y, z, zero, s4);

    for (i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
    {
        gmx_simd4_maskstore3_r(v[i], s4[i]);
    }
#else
    int i, m;

    gmx_simd_store_r(buf_aligned + 0*GMX_SIMD_REAL_WIDTH, x);
    gmx_simd_store_r(buf_aligned + 1*GMX_SIMD_REAL_WIDTH, y);
    gmx_simd_store_r(buf_aligned + 2*GMX_SIMD_REAL_WIDTH, z);

    for (i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
    {
        for (m = 0; m < DIM; m++)
        {
            v[i][m] = buf_aligned[m*GMX_SIMD_REAL_WIDTH + i];
        }
    }
#endif
}

#endif /* GMX_SIMD_HAVE_REAL */

/*! \} */
/*! \endcond */

#endif /* GMX_SIMD_RVEC_OPERATIONS_H */
