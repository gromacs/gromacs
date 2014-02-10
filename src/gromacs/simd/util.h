/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014, by the GROMACS development team, led by
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
#ifndef GMX_SIMD_SIMD_UTIL_H_
#define GMX_SIMD_SIMD_UTIL_H_

#ifdef GMX_SIMD_HAVE_FLOAT

/*! \name Single precision SIMD utility functions
 *
 *  \note In most cases you should use the real-precision functions instead.
 *  \{
 */

/*******************************************
 * SINGLE PRECISION SIMD UTILITY FUNCTIONS *
 *******************************************/

#ifndef GMX_SIMD_HAVE_REDUCE
/*! \brief Computes the sum of the first 2 elements of vector x
 *
 * You should normally call the real-precision routine \ref gmx_simd_reduce2_r.
 *
 * \param x Argument.
 * \param b Aligned temporary memory buffer
 * \return sum of first 2 elements of x
 */
static gmx_inline float
gmx_simd_reduce2_f(gmx_simd_float_t x, float* b)
{
    gmx_simd_store_f(b, x);
    return b[0]+b[1];
}

#if GMX_SIMD_FLOAT_WIDTH == 2 || defined DOXYGEN
/*! \brief Computes the sum of the elements of vector x
 *
 * You should normally call the real-precision routine \ref gmx_simd_reduce_r.
 *
 * \param x Argument.
 * \param b Aligned temporary memory buffer
 * \return sum of elements of x
 */
static gmx_inline float gmx_simd_reduce_f(gmx_simd_float_t x, float* b)
{
    gmx_simd_store_f(b, x);
    return b[0]+b[1];
}
#elif GMX_SIMD_FLOAT_WIDTH == 4
static gmx_inline float gmx_simd_reduce_f(gmx_simd_float_t x, float* b)
{
    gmx_simd_store_f(b, x);
    return b[0]+b[1]+b[2]+b[3];
}
#elif GMX_SIMD_FLOAT_WIDTH == 8
static gmx_inline float gmx_simd_reduce_f(gmx_simd_float_t x, float* b)
{
    gmx_simd_store_f(b, x);
    return b[0]+b[1]+b[2]+b[3]+b[4]+b[5]+b[6]+b[7];
}
#elif GMX_SIMD_FLOAT_WIDTH == 16
/* This is getting ridiculous, SIMD horizontal adds would help,
 * but this is not performance critical (only used to reduce energies)
 */
static gmx_inline float gmx_simd_reduce_f(gmx_simd_float_t x, float* b)
{
    gmx_simd_store_f(b, x);
    return b[0]+b[1]+b[2]+b[3]+b[4]+b[5]+b[6]+b[7]+b[8]+b[9]+b[10]+b[11]+b[12]+b[13]+b[14]+b[15];
}
#else
#error "unsupported simd width configuration"
#endif
#endif /* GMX_SIMD_HAVE_REDUCE */
#endif /* GMX_SIMD_HAVE_FLOAT */

/*! \} */

#ifdef GMX_SIMD_HAVE_DOUBLE

/*! \name Double precision SIMD utility functions
 *
 *  \note In most cases you should use the real-precision functions instead.
 *  \{
 */

/*******************************************
 * DOUBLE PRECISION SIMD UTILITY FUNCTIONS *
 *******************************************/

#ifndef GMX_SIMD_HAVE_REDUCE
/*! \brief Computes the sum of the first 2 elements of vector x
 *
 * \copydetails gmx_simd_reduce2_f
 */
static gmx_inline double
gmx_simd_reduce2_d(gmx_simd_double_t x, double* b)
{
    gmx_simd_store_d(b, x);
    return b[0]+b[1];
}

#if GMX_SIMD_DOUBLE_WIDTH == 2 || defined DOXYGEN
/*! \brief Computes the sum of the elements of vector x
 *
 * \copydetails gmx_simd_reduce_f
 */
static gmx_inline double gmx_simd_reduce_d(gmx_simd_double_t x, double* b)
{
    gmx_simd_store_d(b, x);
    return b[0]+b[1];
}
#elif GMX_SIMD_DOUBLE_WIDTH == 4
static gmx_inline double gmx_simd_reduce_d(gmx_simd_double_t x, double* b)
{
    gmx_simd_store_d(b, x);
    return b[0]+b[1]+b[2]+b[3];
}
#elif GMX_SIMD_DOUBLE_WIDTH == 8
static gmx_inline double gmx_simd_reduce_d(gmx_simd_double_t x, double* b)
{
    gmx_simd_store_d(b, x);
    return b[0]+b[1]+b[2]+b[3]+b[4]+b[5]+b[6]+b[7];
}
#elif GMX_SIMD_DOUBLE_WIDTH == 16
/* This is getting ridiculous, SIMD horizontal adds would help,
 * but this is not performance critical (only used to reduce energies)
 */
static gmx_inline double gmx_simd_reduce_d(gmx_simd_double_t x, double* b)
{
    gmx_simd_store_d(b, x);
    return b[0]+b[1]+b[2]+b[3]+b[4]+b[5]+b[6]+b[7]+b[8]+b[9]+b[10]+b[11]+b[12]+b[13]+b[14]+b[15];
}
#else
#error "unsupported simd width configuration"
#endif
#endif /* GMX_SIMD_HAVE_REDUCE */

/*! \} */

#endif /* GMX_SIMD_HAVE_DOUBLE */

/*! \name SIMD4 utility functions
 *
 *  \{
 */

#ifdef GMX_SIMD4_HAVE_FLOAT

/*************************************************************************
 * SINGLE PRECISION SIMD4 UTILITY FUNCTIONS *
 *************************************************************************/

#ifndef GMX_SIMD_HAVE_REDUCE
/*! \brief Computes the sum of the elements of SIMD4 x
 *
 * You should normally call the real-precision routine \ref gmx_simd4_reduce_r.
 *
 * \param x Argument.
 * \param b Aligned temporary memory buffer
 * \return sum of elements of x
 */
static gmx_inline float
gmx_simd4_reduce_f(gmx_simd4_float_t x, float* b)
{
    gmx_simd4_store_f(b, x);
    return b[0]+b[1]+b[2]+b[3];
}
#endif /* GMX_SIMD_HAVE_REDUCE */

/*! \} */

#endif /* GMX_SIMD4_HAVE_FLOAT */

/*! \name SIMD4 utility functions
 *
 *  \{
 */

#ifdef GMX_SIMD4_HAVE_DOUBLE
/*************************************************************************
 * DOULBE PRECISION SIMD4 UTILITY FUNCTIONS *
 *************************************************************************/

#ifndef GMX_SIMD_HAVE_REDUCE
/*! \brief Computes the sum of the elements of SIMD4 x
 *
 * \copydetails gmx_simd4_reduce_f
 */
static gmx_inline double
gmx_simd4_reduce_d(gmx_simd4_double_t x, double* b)
{
    gmx_simd4_store_d(b, x);
    return b[0]+b[1]+b[2]+b[3];
}
#endif /* GMX_SIMD_HAVE_REDUCE */

/*! \} */

#endif /* GMX_SIMD4_HAVE_DOUBLE */

/* Set defines based on default Gromacs precision */
#ifdef GMX_DOUBLE
/* Documentation in single branch below */
#    define gmx_simd_reduce2_r        gmx_simd_reduce2_d
#    define gmx_simd_reduce_r         gmx_simd_reduce_d
#    define gmx_simd4_reduce_r        gmx_simd4_reduce_d

#else /* GMX_DOUBLE */

/*! \name Real-precision SIMD utility functions
 *
 *  These are the ones you should typically call in Gromacs.
 * \{
 */

/*! \brief Computes the sum of the first 2 elements of vector x
 *
 * \copydetails gmx_simd_reduce2_f
 */
#    define gmx_simd_reduce2_r           gmx_simd_reduce2_f

/*! \brief Computes the sum of the elements of vector x
 *
 * \copydetails gmx_simd_reduce_f
 */
#    define gmx_simd_reduce_r           gmx_simd_reduce_f

/*! \brief Computes the sum of the elements of SIMD4 x
 *
 * \copydetails gmx_simd4_reduce_f
 */
#    define gmx_simd4_reduce_r           gmx_simd4_reduce_f

/*! \} */

#endif /* GMX_DOUBLE */

#endif /* GMX_SIMD_SIMD_UTIL_H_ */
