/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
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
#if GMX_SIMD_HAVE_REAL || defined DOXYGEN

#if GMX_SIMD_HAVE_FLOAT || defined DOXYGEN
/*! \brief SIMD float inner product of multiple float vectors.
 *
 * For normal usage you should always call the real-precision \ref gmx::simdIprod.
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
simdIprodF(SimdFloat ax, SimdFloat ay, SimdFloat az,
           SimdFloat bx, SimdFloat by, SimdFloat bz)
{
    SimdFloat ret;

    ret = simdMulF(ax, bx);
    ret = simdFmaddF(ay, by, ret);
    ret = simdFmaddF(az, bz, ret);

    return ret;
}

/*! \brief SIMD float norm squared of multiple vectors.
 *
 * For normal usage you should always call the real-precision \ref gmx::simdNorm2.
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
static inline SimdFloat gmx_simdcall
simdNorm2F(SimdFloat ax, SimdFloat ay, SimdFloat az)
{
    SimdFloat ret;

    ret = simdMulF(ax, ax);
    ret = simdFmaddF(ay, ay, ret);
    ret = simdFmaddF(az, az, ret);

    return ret;
}

/*! \brief Calculating r^2 is the same as evaluating the norm of dx*dx.
 *
 * For details, see \ref gmx::simdNorm2F.
 */
#define simdCalcRsqF simdNorm2F

/*! \brief SIMD float cross-product of multiple vectors.
 *
 * For normal usage you should always call the real-precision \ref gmx::simdCprod.
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
static inline void gmx_simdcall
simdCprodF(SimdFloat ax, SimdFloat ay, SimdFloat az,
           SimdFloat bx, SimdFloat by, SimdFloat bz,
           SimdFloat *cx, SimdFloat *cy, SimdFloat *cz)
{
    *cx = simdMulF(ay, bz);
    *cx = simdFnmaddF(az, by, *cx);

    *cy = simdMulF(az, bx);
    *cy = simdFnmaddF(ax, bz, *cy);

    *cz = simdMulF(ax, by);
    *cz = simdFnmaddF(ay, bx, *cz);
}
#endif /* GMX_SIMD_HAVE_FLOAT */

#if GMX_SIMD_HAVE_DOUBLE || defined DOXYGEN
/*! \brief SIMD double inner product of multiple double vectors.
 *
 * \copydetails simdIprodF
 */
static inline SimdDouble gmx_simdcall
simdIprodD(SimdDouble ax, SimdDouble ay, SimdDouble az,
           SimdDouble bx, SimdDouble by, SimdDouble bz)
{
    SimdDouble ret;

    ret = simdMulD(ax, bx);
    ret = simdFmaddD(ay, by, ret);
    ret = simdFmaddD(az, bz, ret);

    return ret;
}

/*! \brief SIMD double norm squared of multiple vectors.
 *
 * \copydetails simdNorm2F
 */
static inline SimdDouble gmx_simdcall
simdNorm2D(SimdDouble ax, SimdDouble ay, SimdDouble az)
{
    SimdDouble ret;

    ret = simdMulD(ax, ax);
    ret = simdFmaddD(ay, ay, ret);
    ret = simdFmaddD(az, az, ret);

    return ret;
}

/*! \brief Calculating r^2 is the same as evaluating the norm of dx*dx.
 *
 * For details, see \ref gmx::simdNorm2D.
 */
#define simdCalcRsqD simdNorm2D

/*! \brief SIMD double cross-product of multiple vectors.
 *
 * \copydetails simdCprodF
 */
static inline void gmx_simdcall
simdCprodD(SimdDouble ax, SimdDouble ay, SimdDouble az,
           SimdDouble bx, SimdDouble by, SimdDouble bz,
           SimdDouble *cx, SimdDouble *cy, SimdDouble *cz)
{
    *cx = simdMulD(ay, bz);
    *cx = simdFnmaddD(az, by, *cx);

    *cy = simdMulD(az, bx);
    *cy = simdFnmaddD(ax, bz, *cy);

    *cz = simdMulD(ax, by);
    *cz = simdFnmaddD(ay, bx, *cz);
}
#endif /* GMX_SIMD_HAVE_DOUBLE */


#if GMX_SIMD4_HAVE_FLOAT || defined DOXYGEN
/*! \brief SIMD4 float inner product of four float vectors.
 *
 * \copydetails simdNorm2F
 */
static inline Simd4Float gmx_simdcall
simd4Norm2F(Simd4Float ax, Simd4Float ay, Simd4Float az)
{
    Simd4Float ret;

    ret = simd4MulF(ax, ax);
    ret = simd4FmaddF(ay, ay, ret);
    ret = simd4FmaddF(az, az, ret);

    return ret;
}

/*! \brief Calculating r^2 is the same as evaluating the norm of dx*dx.
 *
 * For details, see \ref gmx::simd4Norm2F
 */
#define simd4CalcRsqF simd4Norm2F

#endif /* GMX_SIMD4_HAVE_FLOAT */

#if GMX_SIMD4_HAVE_DOUBLE || defined DOXYGEN
/*! \brief SIMD4 double norm squared of multiple vectors.
 *
 * \copydetails simdNorm2F
 */
static inline Simd4Double gmx_simdcall
simd4Norm2D(Simd4Double ax, Simd4Double ay, Simd4Double az)
{
    Simd4Double ret;

    ret = simd4MulD(ax, ax);
    ret = simd4FmaddD(ay, ay, ret);
    ret = simd4FmaddD(az, az, ret);

    return ret;
}

/*! \brief Calculating r^2 is the same as evaluating the norm of dx*dx.
 *
 * For details, see \ref gmx::simd4Norm2D.
 */
#define simd4CalcRsqD simd4Norm2D

#endif /* GMX_SIMD4_HAVE_DOUBLE */


#ifdef GMX_DOUBLE
/* Documented for the single branch below */
#    define simdIprod      simdIprodD
#    define simdNorm2      simdNorm2D
#    define simdCalcRsq   simdCalcRsqD
#    define simdCprod      simdCprodD
#    define simd4Norm2     simd4Norm2D
#    define simd4CalcRsq  simd4CalcRsqD
#else /* GMX_DOUBLE */

/*! \brief SIMD real inner product of multiple real vectors.
 *
 * This will call \ref gmx::simdIprodD if GMX_DOUBLE is defined, otherwise
 * \ref gmx::simdIprodF.
 *
 * \copydetails simdIprodF
 */
#    define simdIprod      simdIprodF

/*! \brief SIMD real norm squared of multiple real vectors.
 *
 * This will call \ref gmx::simdNorm2D if GMX_DOUBLE is defined, otherwise
 * \ref gmx::simdNorm2F.
 *
 * \copydetails simdNorm2F
 */
#    define simdNorm2      simdNorm2F

/*! \brief Calculating r^2 is the same as evaluating the norm of dx*dx.
 *
 * This will call \ref gmx::simdCalcRsqD if GMX_DOUBLE is defined, otherwise
 * \ref gmx::simdCalcRsqF.
 *
 * \copydetails simdCalcRsqF
 */
#    define simdCalcRsq   simdCalcRsqF

/*! \brief SIMD real cross-product of multiple real vectors.
 *
 * This will call \ref gmx::simdCprodD if GMX_DOUBLE is defined, otherwise
 * \ref gmx::simdCprodF.
 *
 * \copydetails simdCprodF
 */
#    define simdCprod      simdCprodF

/*! \brief SIMD4 real norm squared of multiple vectors.
 *
 * This will call \ref gmx::simd4Norm2D if GMX_DOUBLE is defined, otherwise
 * \ref gmx::simd4Norm2F.
 *
 * \copydetails simd4Norm2F
 */
#    define simd4Norm2     simd4Norm2F

/*! \brief Calculating r^2 is the same as evaluating the norm of dx*dx.
 *
 * This will call \ref gmx::simd4CalcRsqD if GMX_DOUBLE is defined, otherwise
 * \ref gmx::simd4CalcRsqF.
 *
 * \copydetails simd4CalcRsqF
 */
#    define simd4CalcRsq  simd4CalcRsqF

#endif /* GMX_DOUBLE */

#endif /* GMX_SIMD_HAVE REAL || defined DOXYGEN */

/*! \} */
/*! \endcond */

#endif /* GMX_SIMD */

}      // namespace gmx

#endif /* GMX_SIMD_VECTOR_OPERATIONS_H */
