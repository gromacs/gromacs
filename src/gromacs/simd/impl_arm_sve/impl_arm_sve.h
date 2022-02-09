/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */

/*
 * armv8+sve support to GROMACS was contributed by the Research Organization for
 * Information Science and Technology (RIST).
 * Copyright (c) 2020 Research Organization for Information Science and Technology (RIST).
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
 * We use fixed-size SVE types with __atributte__((arm_sve_vector_bits(...)))
 * to set the size. (See Arm C Language Extensions for SVE - version 6, chapter 3.7.3.)
 */

#endif // GMX_SIMD_IMPL_ARM_SVE_H
