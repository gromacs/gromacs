/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015,2017, by the GROMACS development team, led by
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

#ifndef GMX_SIMD_IMPL_REFERENCE_DEFINITIONS_H
#define GMX_SIMD_IMPL_REFERENCE_DEFINITIONS_H

/*! \libinternal \file
 *
 * \brief Reference SIMD implementation, including SIMD documentation.
 *
 * \author Erik Lindahl <erik.lindahl@scilifelab.se>
 *
 * \ingroup module_simd
 */
namespace gmx
{

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

//! \brief 1 if any SIMD support is present, otherwise 0.
#define GMX_SIMD                                                 1

/*! \brief 1 when SIMD float support is present, otherwise 0
 *
 * You should only use this to specifically check for single precision SIMD,
 * support, even when the rest of Gromacs uses double precision.
 */
#define GMX_SIMD_HAVE_FLOAT                                      1

//! \brief 1 if SIMD double support is present, otherwise 0
#define GMX_SIMD_HAVE_DOUBLE                                     1

//! \brief 1 if the SIMD implementation supports unaligned loads, otherwise 0
#define GMX_SIMD_HAVE_LOADU                                      1

//! \brief 1 if the SIMD implementation supports unaligned stores, otherwise 0
#define GMX_SIMD_HAVE_STOREU                                     1

/*! \brief 1 if the SIMD implementation has fused-multiply add hardware
 *
 * \note All the fused multiply-add functions are always available and can be
 *       used in any code (by executing separate multiply and add ops), but in
 *       a few very tight loops you might be able to save a few instructions
 *       with a separate non-FMA code path.
 */
#define GMX_SIMD_HAVE_FMA                                        0

//! \brief 1 if SIMD impl has logical operations on floating-point data, otherwise 0
#define GMX_SIMD_HAVE_LOGICAL                                    1

//! \brief Support for extracting integers from \ref gmx::SimdFInt32 (1/0 for present/absent)
#define GMX_SIMD_HAVE_FINT32_EXTRACT                             1

//! \brief 1 if SIMD logical ops are supported for \ref gmx::SimdFInt32, otherwise 0
#define GMX_SIMD_HAVE_FINT32_LOGICAL                             1

//! \brief 1 if SIMD arithmetic ops are supported for \ref gmx::SimdFInt32, otherwise 0
#define GMX_SIMD_HAVE_FINT32_ARITHMETICS                         1

//! \brief Support for extracting integer from \ref gmx::SimdDInt32 (1/0 for present/absent)
#define GMX_SIMD_HAVE_DINT32_EXTRACT                             1

//! \brief 1 if logical operations are supported for \ref gmx::SimdDInt32, otherwise 0
#define GMX_SIMD_HAVE_DINT32_LOGICAL                             1

//! \brief 1 if SIMD arithmetic ops are supported for \ref gmx::SimdDInt32, otherwise 0
#define GMX_SIMD_HAVE_DINT32_ARITHMETICS                         1

/*! \brief 1 if implementation provides single precision copysign()
 *
 *  Only used in simd_math.h to selectively override the generic implementation.
 */
#define GMX_SIMD_HAVE_NATIVE_COPYSIGN_FLOAT                      0

/*! \brief 1 if implementation provides single precision 1/sqrt(x) N-R iterations faster than simd_math.h
 *
 *  Only used in simd_math.h to selectively override the generic implementation.
 */
#define GMX_SIMD_HAVE_NATIVE_RSQRT_ITER_FLOAT                    0

/*! \brief 1 if implementation provides single precision 1/x N-R iterations faster than simd_math.h
 *
 *  Only used in simd_math.h to selectively override the generic implementation.
 */
#define GMX_SIMD_HAVE_NATIVE_RCP_ITER_FLOAT                      0

/*! \brief 1 if implementation provides single precision log() faster than simd_math.h
 *
 *  Only used in simd_math.h to selectively override the generic implementation.
 */
#define GMX_SIMD_HAVE_NATIVE_LOG_FLOAT                           0

/*! \brief 1 if implementation provides single precision exp2() faster than simd_math.h
 *
 *  Only used in simd_math.h to selectively override the generic implementation.
 */
#define GMX_SIMD_HAVE_NATIVE_EXP2_FLOAT                          0

/*! \brief 1 if implementation provides single precision exp() faster than simd_math.h
 *
 *  Only used in simd_math.h to selectively override the generic implementation.
 */
#define GMX_SIMD_HAVE_NATIVE_EXP_FLOAT                           0

/*! \brief 1 if implementation provides double precision copysign()
 *
 *  Only used in simd_math.h to selectively override the generic implementation.
 */
#define GMX_SIMD_HAVE_NATIVE_COPYSIGN_DOUBLE                     0

/*! \brief 1 if implementation provides double precision 1/sqrt(x) N-R iterations faster than simd_math.h
 *
 *  Only used in simd_math.h to selectively override the generic implementation.
 */
#define GMX_SIMD_HAVE_NATIVE_RSQRT_ITER_DOUBLE                   0

/*! \brief 1 if implementation provides double precision 1/x N-R iterations faster than simd_math.h
 *
 *  Only used in simd_math.h to selectively override the generic implementation.
 */
#define GMX_SIMD_HAVE_NATIVE_RCP_ITER_DOUBLE                     0

/*! \brief 1 if implementation provides double precision log() faster than simd_math.h
 *
 *  Only used in simd_math.h to selectively override the generic implementation.
 */
#define GMX_SIMD_HAVE_NATIVE_LOG_DOUBLE                          0

/*! \brief 1 if implementation provides double precision exp2() faster than simd_math.h
 *
 *  Only used in simd_math.h to selectively override the generic implementation.
 */
#define GMX_SIMD_HAVE_NATIVE_EXP2_DOUBLE                         0

/*! \brief 1 if implementation provides double precision exp() faster than simd_math.h
 *
 *  Only used in simd_math.h to selectively override the generic implementation.
 */
#define GMX_SIMD_HAVE_NATIVE_EXP_DOUBLE                          0

//! \brief 1 if \ref gmx::gatherLoadUBySimdIntTranspose is present, otherwise 0
#define GMX_SIMD_HAVE_GATHER_LOADU_BYSIMDINT_TRANSPOSE_FLOAT     1

//! \brief 1 if \ref gmx::gatherLoadUBySimdIntTranspose is present, otherwise 0
#define GMX_SIMD_HAVE_GATHER_LOADU_BYSIMDINT_TRANSPOSE_DOUBLE    1

//! \brief 1 if float half-register load/store/reduce utils present, otherwise 0
#define GMX_SIMD_HAVE_HSIMD_UTIL_FLOAT                           1

//! \brief 1 if double half-register load/store/reduce utils present, otherwise 0
#define GMX_SIMD_HAVE_HSIMD_UTIL_DOUBLE                          1

#ifdef GMX_SIMD_REF_FLOAT_WIDTH
#    define GMX_SIMD_FLOAT_WIDTH                                 GMX_SIMD_REF_FLOAT_WIDTH
#else
//! \brief Width of the \ref gmx::SimdFloat datatype
#    define GMX_SIMD_FLOAT_WIDTH                                 4
#endif

#ifdef GMX_SIMD_REF_DOUBLE_WIDTH
#    define GMX_SIMD_DOUBLE_WIDTH                                GMX_SIMD_REF_DOUBLE_WIDTH
#else
//! \brief Width of the \ref gmx::SimdDouble datatype
#    define GMX_SIMD_DOUBLE_WIDTH                                4
#endif

#if GMX_SIMD_FLOAT_WIDTH >= 8 || defined DOXYGEN //set in simd.h for GMX_SIMD_FLOAT_WIDTH<=4
//! \brief 1 if float 4xN load utils present, otherwise 0
#define GMX_SIMD_HAVE_4NSIMD_UTIL_FLOAT                          1
#endif

#if GMX_SIMD_DOUBLE_WIDTH >= 8 || defined DOXYGEN //set in simd.h for GMX_SIMD_DOUBLE_WIDTH<=4
//! \brief 1 if double 4xN load utils present, otherwise 0
#define GMX_SIMD_HAVE_4NSIMD_UTIL_DOUBLE                         1
#endif

//! \brief 1 if implementation provides \ref gmx::Simd4Float, otherwise 0.
#define GMX_SIMD4_HAVE_FLOAT                                     1

//! \brief 1 if the implementation provides \ref gmx::Simd4Double, otherwise 0.
#define GMX_SIMD4_HAVE_DOUBLE                                    1

//! \brief Width of the \ref gmx::SimdFInt32 datatype.
#define GMX_SIMD_FINT32_WIDTH                                    GMX_SIMD_FLOAT_WIDTH

//! \brief Width of the \ref gmx::SimdDInt32 datatype.
#define GMX_SIMD_DINT32_WIDTH                                    GMX_SIMD_DOUBLE_WIDTH

//! \brief The SIMD4 type is always four units wide, but this makes code more explicit
#define GMX_SIMD4_WIDTH                                          4

//! \brief Accuracy of SIMD 1/sqrt(x) lookup. Used to determine number of iterations.
#define GMX_SIMD_RSQRT_BITS                                      23

//! \brief Accuracy of SIMD 1/x lookup. Used to determine number of iterations.
#define GMX_SIMD_RCP_BITS                                        23

//! \}

//! \}

//! \endcond

}
#endif // GMX_SIMD_IMPL_REFERENCE_DEFINITIONS_H
