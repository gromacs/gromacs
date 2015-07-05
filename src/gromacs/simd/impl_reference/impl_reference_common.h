/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015, by the GROMACS development team, led by
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

#ifndef GMX_SIMD_IMPL_REFERENCE_COMMON_H
#define GMX_SIMD_IMPL_REFERENCE_COMMON_H

/*! \libinternal \file
 *
 * \brief Reference SIMD implementation, including SIMD documentation.
 *
 * \author Erik Lindahl <erik.lindahl@scilifelab.se>
 *
 * \ingroup module_simd
 */

/*! \cond libapi */
/*! \addtogroup module_simd */
/*! \{ */

/*! \name SIMD implementation capability definitions
 *  \{
 */

/* We list all the capability definitions in the main/wrapper SIMD header for
 * each implementation, so we don't forget to set values for unsupported
 * features to 0.
 */

/*! \brief 1 if any SIMD support is present, otherwise 0 */
#define GMX_SIMD                             0

/*! \brief
 * 1 when SIMD float support is present, otherwise 0
 *
 * You should only use this to specifically check for single precision SIMD,
 * support, even when the rest of Gromacs uses double precision.
 * \sa GMX_SIMD_HAVE_REAL, GMX_SIMD_HAVE_DOUBLE
 */
#define GMX_SIMD_HAVE_FLOAT                  1

/*! \brief 1 if SIMD double support is present, otherwise 0 */
#define GMX_SIMD_HAVE_DOUBLE                 1

/*! \brief 1 if the SIMD implementation supports unaligned loads, otherwise 0 */
#define GMX_SIMD_HAVE_LOADU                  1

/*! \brief 1 if the SIMD implementation supports unaligned stores, otherwise 0 */
#define GMX_SIMD_HAVE_STOREU                 1

/*! \brief 1 if SIMD impl has logical operations on floating-point data, otherwise 0 */
#define GMX_SIMD_HAVE_LOGICAL                1

/*! \brief 1 if SIMD fused multiply-add uses hardware instructions, otherwise 0 */
#define GMX_SIMD_HAVE_FMA                    0

/*! \brief 1 if the SIMD fraction has a direct hardware instruction, otherwise 0 */
#define GMX_SIMD_HAVE_FRACTION               0

/*! \brief 1 if the SIMD implementation has \ref gmx_simd_fint32_t, otherwise 0 */
#define GMX_SIMD_HAVE_FINT32                 1

/*! \brief Support for extracting integers from \ref gmx_simd_fint32_t (1/0 for present/absent) */
#define GMX_SIMD_HAVE_FINT32_EXTRACT         1

/*! \brief 1 if SIMD logical ops are supported for \ref gmx_simd_fint32_t, otherwise 0 */
#define GMX_SIMD_HAVE_FINT32_LOGICAL         1

/*! \brief 1 if SIMD arithmetic ops are supported for \ref gmx_simd_fint32_t, otherwise 0 */
#define GMX_SIMD_HAVE_FINT32_ARITHMETICS     1

/*! \brief 1 if the SIMD implementation has \ref gmx_simd_dint32_t, otherwise 0.
 *
 * \note The Gromacs SIMD module works entirely with 32 bit integers, both
 * in single and double precision, since some platforms do not support 64 bit
 * SIMD integers at all. In particular, this means it is up to each
 * implementation to get this working even if the architectures internal
 * representation uses 64 bit integers when converting to/from double SIMD
 * variables. For now we will try HARD to use conversions, packing or shuffling
 * so the integer datatype has the same width as the floating-point type, i.e.
 * if you use double precision SIMD with a width of 8, we want the integers
 * we work with to also use a SIMD width of 8 to make it easy to load/store
 * indices from arrays. This refers entirely to the function calls
 * and how many integers we load/store in one call; the actual SIMD registers
 * might be wider for integers internally (e.g. on x86 gmx_simd_dint32_t will
 * only fill half the register), but this is none of the user's business.
 * While this works for all current architectures, and we think it will work
 * for future ones, we might have to alter this decision in the future. To
 * avoid rewriting every single instance that refers to the SIMD width we still
 * provide separate defines for the width of SIMD integer variables that you
 * should use.
 */
#define GMX_SIMD_HAVE_DINT32                 1

/*! \brief Support for extracting integer from \ref gmx_simd_dint32_t (1/0 for present/absent) */
#define GMX_SIMD_HAVE_DINT32_EXTRACT         1

/*! \brief 1 if logical operations are supported for \ref gmx_simd_dint32_t, otherwise 0 */
#define GMX_SIMD_HAVE_DINT32_LOGICAL         1

/*! \brief 1 if SIMD arithmetic ops are supported for \ref gmx_simd_dint32_t, otherwise 0 */
#define GMX_SIMD_HAVE_DINT32_ARITHMETICS     1

/*! \brief 1 if implementation provides \ref gmx_simd4_float_t, otherwise 0 */
#define GMX_SIMD4_HAVE_FLOAT                 1

/*! \brief 1 if the implementation provides \ref gmx_simd4_double_t, otherwise 0 */
#define GMX_SIMD4_HAVE_DOUBLE                1

#ifdef GMX_SIMD_REF_FLOAT_WIDTH
#    define GMX_SIMD_FLOAT_WIDTH             GMX_SIMD_REF_FLOAT_WIDTH
#else
/*! \brief Width of the \ref gmx_simd_float_t datatype. */
#    define GMX_SIMD_FLOAT_WIDTH             4
#endif

#ifdef GMX_SIMD_REF_DOUBLE_WIDTH
#    define GMX_SIMD_DOUBLE_WIDTH            GMX_SIMD_REF_DOUBLE_WIDTH
#else
/*! \brief Width of the \ref gmx_simd_double_t datatype. */
#    define GMX_SIMD_DOUBLE_WIDTH            4
#endif

/*! \brief Width of the \ref gmx_simd_fint32_t datatype. */
#define GMX_SIMD_FINT32_WIDTH            GMX_SIMD_FLOAT_WIDTH

/*! \brief Width of the \ref gmx_simd_dint32_t datatype. */
#define GMX_SIMD_DINT32_WIDTH            GMX_SIMD_DOUBLE_WIDTH

/*! \brief The SIMD4 type is always four units wide, but this makes code more explicit */
#define GMX_SIMD4_WIDTH                      4

/*! \brief Accuracy of SIMD 1/sqrt(x) lookup. Used to determine number of iterations. */
#define GMX_SIMD_RSQRT_BITS                 23

/*! \brief Accuracy of SIMD 1/x lookup. Used to determine number of iterations. */
#define GMX_SIMD_RCP_BITS                   23

/*! \} */

/*! \} */
/*! \endcond */

#endif /* GMX_SIMD_IMPL_REFERENCE_COMMON_H */
