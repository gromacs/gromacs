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

/*! \cond libapi */
/*! \addtogroup module_simd */
/*! \{ */

/* This check is not actually required, but it must be true if the
 * code below actualy declares anything, and it makes it easy for
 * check-source to know that this file depends on simd.h (though
 * symbols like GMX_SIMD_HAVE_FLOAT are actually defined in its
 * implementation headers). */
#if (defined GMX_SIMD_HAVE_REAL) || (defined DOXYGEN)

#if (defined GMX_SIMD_HAVE_FLOAT) || (defined DOXYGEN)
/*! \brief SIMD float inner product of multiple float vectors.
 *
 * For normal usage you should always call the real-precision \ref gmx_simd_iprod_r.
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
static gmx_inline gmx_simd_float_t gmx_simdcall
gmx_simd_iprod_f(gmx_simd_float_t ax, gmx_simd_float_t ay, gmx_simd_float_t az,
                 gmx_simd_float_t bx, gmx_simd_float_t by, gmx_simd_float_t bz)
{
    gmx_simd_float_t ret;

    ret = gmx_simd_mul_f(ax, bx);
    ret = gmx_simd_fmadd_f(ay, by, ret);
    ret = gmx_simd_fmadd_f(az, bz, ret);

    return ret;
}

/*! \brief SIMD float norm squared of multiple vectors.
 *
 * For normal usage you should always call the real-precision \ref gmx_simd_norm2_r.
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
static gmx_inline gmx_simd_float_t gmx_simdcall
gmx_simd_norm2_f(gmx_simd_float_t ax, gmx_simd_float_t ay, gmx_simd_float_t az)
{
    gmx_simd_float_t ret;

    ret = gmx_simd_mul_f(ax, ax);
    ret = gmx_simd_fmadd_f(ay, ay, ret);
    ret = gmx_simd_fmadd_f(az, az, ret);

    return ret;
}

/*! \brief Calculating r^2 is the same as evaluating the norm of dx*dx.
 *
 * For details, see \ref gmx_simd_norm2_f.
 */
#define gmx_simd_calc_rsq_f gmx_simd_norm2_f

/*! \brief SIMD float cross-product of multiple vectors.
 *
 * For normal usage you should always call the real-precision \ref gmx_simd_cprod_r.
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
static gmx_inline void gmx_simdcall
gmx_simd_cprod_f(gmx_simd_float_t ax, gmx_simd_float_t ay, gmx_simd_float_t az,
                 gmx_simd_float_t bx, gmx_simd_float_t by, gmx_simd_float_t bz,
                 gmx_simd_float_t *cx, gmx_simd_float_t *cy, gmx_simd_float_t *cz)
{
    *cx = gmx_simd_mul_f(ay, bz);
    *cx = gmx_simd_fnmadd_f(az, by, *cx);

    *cy = gmx_simd_mul_f(az, bx);
    *cy = gmx_simd_fnmadd_f(ax, bz, *cy);

    *cz = gmx_simd_mul_f(ax, by);
    *cz = gmx_simd_fnmadd_f(ay, bx, *cz);
}
#endif /* GMX_SIMD_HAVE_FLOAT */

#if (defined GMX_SIMD_HAVE_DOUBLE) || (defined DOXYGEN)
/*! \brief SIMD double inner product of multiple double vectors.
 *
 * \copydetails gmx_simd_iprod_f
 */
static gmx_inline gmx_simd_double_t gmx_simdcall
gmx_simd_iprod_d(gmx_simd_double_t ax, gmx_simd_double_t ay, gmx_simd_double_t az,
                 gmx_simd_double_t bx, gmx_simd_double_t by, gmx_simd_double_t bz)
{
    gmx_simd_double_t ret;

    ret = gmx_simd_mul_d(ax, bx);
    ret = gmx_simd_fmadd_d(ay, by, ret);
    ret = gmx_simd_fmadd_d(az, bz, ret);

    return ret;
}

/*! \brief SIMD double norm squared of multiple vectors.
 *
 * \copydetails gmx_simd_norm2_f
 */
static gmx_inline gmx_simd_double_t gmx_simdcall
gmx_simd_norm2_d(gmx_simd_double_t ax, gmx_simd_double_t ay, gmx_simd_double_t az)
{
    gmx_simd_double_t ret;

    ret = gmx_simd_mul_d(ax, ax);
    ret = gmx_simd_fmadd_d(ay, ay, ret);
    ret = gmx_simd_fmadd_d(az, az, ret);

    return ret;
}

/*! \brief Calculating r^2 is the same as evaluating the norm of dx*dx.
 *
 * For details, see \ref gmx_simd_norm2_d.
 */
#define gmx_simd_calc_rsq_d gmx_simd_norm2_d

/*! \brief SIMD double cross-product of multiple vectors.
 *
 * \copydetails gmx_simd_cprod_f
 */
static gmx_inline void gmx_simdcall
gmx_simd_cprod_d(gmx_simd_double_t ax, gmx_simd_double_t ay, gmx_simd_double_t az,
                 gmx_simd_double_t bx, gmx_simd_double_t by, gmx_simd_double_t bz,
                 gmx_simd_double_t *cx, gmx_simd_double_t *cy, gmx_simd_double_t *cz)
{
    *cx = gmx_simd_mul_d(ay, bz);
    *cx = gmx_simd_fnmadd_d(az, by, *cx);

    *cy = gmx_simd_mul_d(az, bx);
    *cy = gmx_simd_fnmadd_d(ax, bz, *cy);

    *cz = gmx_simd_mul_d(ax, by);
    *cz = gmx_simd_fnmadd_d(ay, bx, *cz);
}
#endif /* GMX_SIMD_HAVE_DOUBLE */


#if (defined GMX_SIMD4_HAVE_FLOAT) || (defined DOXYGEN)
/*! \brief SIMD4 float inner product of four float vectors.
 *
 * \copydetails gmx_simd_norm2_f
 */
static gmx_inline gmx_simd4_float_t gmx_simdcall
gmx_simd4_norm2_f(gmx_simd4_float_t ax, gmx_simd4_float_t ay, gmx_simd4_float_t az)
{
    gmx_simd4_float_t ret;

    ret = gmx_simd4_mul_f(ax, ax);
    ret = gmx_simd4_fmadd_f(ay, ay, ret);
    ret = gmx_simd4_fmadd_f(az, az, ret);

    return ret;
}

/*! \brief Calculating r^2 is the same as evaluating the norm of dx*dx.
 *
 * For details, see \ref gmx_simd4_norm2_f
 */
#define gmx_simd4_calc_rsq_f gmx_simd4_norm2_f

#endif /* GMX_SIMD4_HAVE_FLOAT */

#if (defined GMX_SIMD4_HAVE_DOUBLE)  || (defined DOXYGEN)
/*! \brief SIMD4 double norm squared of multiple vectors.
 *
 * \copydetails gmx_simd_norm2_f
 */
static gmx_inline gmx_simd4_double_t gmx_simdcall
gmx_simd4_norm2_d(gmx_simd4_double_t ax, gmx_simd4_double_t ay, gmx_simd4_double_t az)
{
    gmx_simd4_double_t ret;

    ret = gmx_simd4_mul_d(ax, ax);
    ret = gmx_simd4_fmadd_d(ay, ay, ret);
    ret = gmx_simd4_fmadd_d(az, az, ret);

    return ret;
}

/*! \brief Calculating r^2 is the same as evaluating the norm of dx*dx.
 *
 * For details, see \ref gmx_simd4_norm2_d.
 */
#define gmx_simd4_calc_rsq_d gmx_simd4_norm2_d

#endif /* GMX_SIMD4_HAVE_DOUBLE */


#ifdef GMX_DOUBLE
/* Documented for the single branch below */
#    define gmx_simd_iprod_r      gmx_simd_iprod_d
#    define gmx_simd_norm2_r      gmx_simd_norm2_d
#    define gmx_simd_calc_rsq_r   gmx_simd_calc_rsq_d
#    define gmx_simd_cprod_r      gmx_simd_cprod_d
#    define gmx_simd4_norm2_r     gmx_simd4_norm2_d
#    define gmx_simd4_calc_rsq_r  gmx_simd4_calc_rsq_d
#else /* GMX_DOUBLE */

/*! \brief SIMD real inner product of multiple real vectors.
 *
 * This will call \ref gmx_simd_iprod_d if GMX_DOUBLE is defined, otherwise
 * \ref gmx_simd_iprod_f.
 *
 * \copydetails gmx_simd_iprod_f
 */
#    define gmx_simd_iprod_r      gmx_simd_iprod_f

/*! \brief SIMD real norm squared of multiple real vectors.
 *
 * This will call \ref gmx_simd_norm2_d if GMX_DOUBLE is defined, otherwise
 * \ref gmx_simd_norm2_f.
 *
 * \copydetails gmx_simd_norm2_f
 */
#    define gmx_simd_norm2_r      gmx_simd_norm2_f

/*! \brief Calculating r^2 is the same as evaluating the norm of dx*dx.
 *
 * This will call \ref gmx_simd_calc_rsq_d if GMX_DOUBLE is defined, otherwise
 * \ref gmx_simd_calc_rsq_f.
 *
 * \copydetails gmx_simd_calc_rsq_f
 */
#    define gmx_simd_calc_rsq_r   gmx_simd_calc_rsq_f

/*! \brief SIMD real cross-product of multiple real vectors.
 *
 * This will call \ref gmx_simd_cprod_d if GMX_DOUBLE is defined, otherwise
 * \ref gmx_simd_cprod_f.
 *
 * \copydetails gmx_simd_cprod_f
 */
#    define gmx_simd_cprod_r      gmx_simd_cprod_f

/*! \brief SIMD4 real norm squared of multiple vectors.
 *
 * This will call \ref gmx_simd4_norm2_d if GMX_DOUBLE is defined, otherwise
 * \ref gmx_simd4_norm2_f.
 *
 * \copydetails gmx_simd4_norm2_f
 */
#    define gmx_simd4_norm2_r     gmx_simd4_norm2_f

/*! \brief Calculating r^2 is the same as evaluating the norm of dx*dx.
 *
 * This will call \ref gmx_simd4_calc_rsq_d if GMX_DOUBLE is defined, otherwise
 * \ref gmx_simd4_calc_rsq_f.
 *
 * \copydetails gmx_simd4_calc_rsq_f
 */
#    define gmx_simd4_calc_rsq_r  gmx_simd4_calc_rsq_f

#endif /* GMX_DOUBLE */

#endif /* (defined GMX_SIMD_HAVE REAL) || (defined DOXYGEN) */

/*! \} */
/*! \endcond */

#endif /* GMX_SIMD_VECTOR_OPERATIONS_H */
