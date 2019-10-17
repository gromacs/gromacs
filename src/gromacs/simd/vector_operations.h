/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013,2014,2015,2017,2019, by the GROMACS development team, led by
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
 * \brief SIMD operations corresponding to Gromacs rvec datatypes.
 *
 * \author Erik Lindahl <erik.lindahl@scilifelab.se>
 *
 * \inlibraryapi
 * \ingroup module_simd
 */

#ifndef GMX_SIMD_VECTOR_OPERATIONS_H
#define GMX_SIMD_VECTOR_OPERATIONS_H

#include "config.h"

#include "gromacs/simd/simd.h"

namespace gmx
{

#if GMX_SIMD

/*! \cond libapi */
/*! \addtogroup module_simd */
/*! \{ */

/* This check is not actually required, but it must be true if the
 * code below actualy declares anything, and it makes it easy for
 * check-source to know that this file depends on simd.h (though
 * symbols like GMX_SIMD_HAVE_FLOAT are actually defined in its
 * implementation headers). */
#    if GMX_SIMD_HAVE_REAL || defined DOXYGEN

#        if GMX_SIMD_HAVE_FLOAT || defined DOXYGEN
/*! \brief SIMD float inner product of multiple float vectors.
 *
 * \param ax X components of first vectors
 * \param ay Y components of first vectors
 * \param az Z components of first vectors
 * \param bx X components of second vectors
 * \param by Y components of second vectors
 * \param bz Z components of second vectors
 *
 * \return Element i will be res[i] = ax[i]*bx[i]+ay[i]*by[i]+az[i]*bz[i].
 *
 * \note The SIMD part is that we calculate many scalar products in one call.
 */
static inline SimdFloat gmx_simdcall
                        iprod(SimdFloat ax, SimdFloat ay, SimdFloat az, SimdFloat bx, SimdFloat by, SimdFloat bz)
{
    SimdFloat ret;

    ret = ax * bx;
    ret = ay * by + ret;
    ret = az * bz + ret;

    return ret;
}

/*! \brief SIMD float norm squared of multiple vectors.
 *
 * \param ax X components of vectors
 * \param ay Y components of vectors
 * \param az Z components of vectors
 *
 * \return Element i will be res[i] = ax[i]*ax[i]+ay[i]*ay[i]+az[i]*az[i].
 *
 * \note This corresponds to the scalar product of the vector with itself, but
 * the compiler might be able to optimize it better with identical vectors.
 */
static inline SimdFloat gmx_simdcall norm2(SimdFloat ax, SimdFloat ay, SimdFloat az)
{
    SimdFloat ret;

    ret = ax * ax;
    ret = ay * ay + ret;
    ret = az * az + ret;

    return ret;
}

/*! \brief SIMD float cross-product of multiple vectors.
 *
 * \param ax X components of first vectors
 * \param ay Y components of first vectors
 * \param az Z components of first vectors
 * \param bx X components of second vectors
 * \param by Y components of second vectors
 * \param bz Z components of second vectors
 * \param[out] cx X components of cross product vectors
 * \param[out] cy Y components of cross product vectors
 * \param[out] cz Z components of cross product vectors
 *
 * \returns void
 *
 * This calculates C = A x B, where the cross denotes the cross product.
 * The arguments x/y/z denotes the different components, and each element
 * corresponds to a separate vector.
 */
static inline void gmx_simdcall cprod(SimdFloat  ax,
                                      SimdFloat  ay,
                                      SimdFloat  az,
                                      SimdFloat  bx,
                                      SimdFloat  by,
                                      SimdFloat  bz,
                                      SimdFloat* cx,
                                      SimdFloat* cy,
                                      SimdFloat* cz)
{
    *cx = ay * bz;
    *cx = fnma(az, by, *cx);

    *cy = az * bx;
    *cy = fnma(ax, bz, *cy);

    *cz = ax * by;
    *cz = fnma(ay, bx, *cz);
}
#        endif // GMX_SIMD_HAVE_FLOAT

#        if GMX_SIMD_HAVE_DOUBLE || defined DOXYGEN
/*! \brief SIMD double inner product of multiple double vectors.
 *
 * \param ax X components of first vectors
 * \param ay Y components of first vectors
 * \param az Z components of first vectors
 * \param bx X components of second vectors
 * \param by Y components of second vectors
 * \param bz Z components of second vectors
 *
 * \return Element i will be res[i] = ax[i]*bx[i]+ay[i]*by[i]+az[i]*bz[i].
 *
 * \note The SIMD part is that we calculate many scalar products in one call.
 */
static inline SimdDouble gmx_simdcall
                         iprod(SimdDouble ax, SimdDouble ay, SimdDouble az, SimdDouble bx, SimdDouble by, SimdDouble bz)
{
    SimdDouble ret;

    ret = ax * bx;
    ret = ay * by + ret;
    ret = az * bz + ret;

    return ret;
}

/*! \brief SIMD double norm squared of multiple vectors.
 *
 * \param ax X components of vectors
 * \param ay Y components of vectors
 * \param az Z components of vectors
 *
 * \return Element i will be res[i] = ax[i]*ax[i]+ay[i]*ay[i]+az[i]*az[i].
 *
 * \note This corresponds to the scalar product of the vector with itself, but
 * the compiler might be able to optimize it better with identical vectors.
 */
static inline SimdDouble gmx_simdcall norm2(SimdDouble ax, SimdDouble ay, SimdDouble az)
{
    SimdDouble ret;

    ret = ax * ax;
    ret = ay * ay + ret;
    ret = az * az + ret;

    return ret;
}

/*! \brief SIMD double cross-product of multiple vectors.
 *
 * \param ax X components of first vectors
 * \param ay Y components of first vectors
 * \param az Z components of first vectors
 * \param bx X components of second vectors
 * \param by Y components of second vectors
 * \param bz Z components of second vectors
 * \param[out] cx X components of cross product vectors
 * \param[out] cy Y components of cross product vectors
 * \param[out] cz Z components of cross product vectors
 *
 * \returns void
 *
 * This calculates C = A x B, where the cross denotes the cross product.
 * The arguments x/y/z denotes the different components, and each element
 * corresponds to a separate vector.
 */
static inline void gmx_simdcall cprod(SimdDouble  ax,
                                      SimdDouble  ay,
                                      SimdDouble  az,
                                      SimdDouble  bx,
                                      SimdDouble  by,
                                      SimdDouble  bz,
                                      SimdDouble* cx,
                                      SimdDouble* cy,
                                      SimdDouble* cz)
{
    *cx = ay * bz;
    *cx = *cx - az * by;

    *cy = az * bx;
    *cy = *cy - ax * bz;

    *cz = ax * by;
    *cz = *cz - ay * bx;
}
#        endif // GMX_SIMD_HAVE_DOUBLE


#        if GMX_SIMD4_HAVE_FLOAT || defined DOXYGEN
/*! \brief SIMD4 float norm squared of multiple vectors.
 *
 * \param ax X components of vectors
 * \param ay Y components of vectors
 * \param az Z components of vectors
 *
 * \return Element i will be res[i] = ax[i]*ax[i]+ay[i]*ay[i]+az[i]*az[i].
 *
 * \note This corresponds to the scalar product of the vector with itself, but
 * the compiler might be able to optimize it better with identical vectors.
 */
static inline Simd4Float gmx_simdcall norm2(Simd4Float ax, Simd4Float ay, Simd4Float az)
{
    Simd4Float ret;

    ret = ax * ax;
    ret = ay * ay + ret;
    ret = az * az + ret;

    return ret;
}

#        endif // GMX_SIMD4_HAVE_FLOAT

#        if GMX_SIMD4_HAVE_DOUBLE || defined DOXYGEN
/*! \brief SIMD4 double norm squared of multiple vectors.
 *
 * \param ax X components of vectors
 * \param ay Y components of vectors
 * \param az Z components of vectors
 *
 * \return Element i will be res[i] = ax[i]*ax[i]+ay[i]*ay[i]+az[i]*az[i].
 *
 * \note This corresponds to the scalar product of the vector with itself, but
 * the compiler might be able to optimize it better with identical vectors.
 */
static inline Simd4Double gmx_simdcall norm2(Simd4Double ax, Simd4Double ay, Simd4Double az)
{
    Simd4Double ret;

    ret = ax * ax;
    ret = ay * ay + ret;
    ret = az * az + ret;

    return ret;
}

#        endif // GMX_SIMD4_HAVE_DOUBLE

#    endif // GMX_SIMD_HAVE REAL || defined DOXYGEN

/*! \} */
/*! \endcond */

#endif // GMX_SIMD

} // namespace gmx

#endif // GMX_SIMD_VECTOR_OPERATIONS_H
