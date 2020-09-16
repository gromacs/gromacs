/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2020 Research Organization for Information Science and Technology (RIST).
 * Copyright (c) 2020, by the GROMACS development team, led by
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

/*
 * armv8+sve support to GROMACS was contributed by the Research Organization for
 * Information Science and Technology (RIST).
 */

#ifndef GMX_SIMD_IMPL_ARM_SVE_H
#define GMX_SIMD_IMPL_ARM_SVE_H

#include "impl_arm_sve_definitions.h"
#include "impl_arm_sve_general.h"
#include "impl_arm_sve_simd4_double.h"
#include "impl_arm_sve_simd4_float.h"
#include "impl_arm_sve_simd_double.h"
#include "impl_arm_sve_simd_float.h"
#include "impl_arm_sve_util_double.h"
#include "impl_arm_sve_util_float.h"

/*
 * armv8+sve support is implemented via the ARM C Language Extensions (ACLE)
 * that introduces scalable/sizeless types such as svfloat32_t (vector of FP32),
 * svint32_t (vector of int) or svbool_t (bit mask).
 * These sizeless types are not a fit for GROMACS that needs to know the vector
 * length at cmake time. GCC 10 and later have the -msve-vector-bits=<len> option
 * in order to fix the vector length at compile time. This feature is currently
 * planned in LLVM 12.
 * Even with this option, svfloat32_t is considered as sizeless
 * by the compiler. However, a float __attribute((vector_size(len/8))) variable is
 * automagically converted as a svfloat32_t and can hence be used to invoke the SIMD
 * intrinsics specified by ACLE.
 * Unfortunately, there is no such thing for svbool_t, so a bit mask (SimdBool) is implemented
 * as a vector of integers.
 */
#endif // GMX_SIMD_IMPL_ARM_SVE_H
