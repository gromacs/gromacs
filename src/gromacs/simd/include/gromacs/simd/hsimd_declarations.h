/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2022- The GROMACS Authors
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

#ifndef GMX_SIMD_HSIMD_DECLARATIONS_H
#define GMX_SIMD_HSIMD_DECLARATIONS_H

/*! \libinternal \file
 * \brief Declares all Hsimd functions that are not supported
 *
 * Including the file simd.h present in the same module conditionally gives access
 * to SIMD functions that load, store and operate on full SIMD width data. On some
 * architectures functions are available for loading and storing half SIMD width
 * data (and duplicating it); this functionality is termed "HSIMD".
 * Including this file is useful for avoiding C-preprocessing conditionals
 * of the form: if GMX_SIMD_HAVE_HSIMD_UTIL_REAL
 * in templated functions that use Hsimd functions within constexpr conditionals.
 *
 * This file can be included unconditionally, as it only declares the HSIMD functions
 * when they are not defined through including simd.h.
 *
 * \author Berk Hess <hess@kth.se>
 *
 * \ingroup module_simd
 */

#include "gromacs/simd/simd.h"

namespace gmx
{

#if GMX_SIMD_HAVE_FLOAT
#    if !GMX_SIMD_HAVE_HSIMD_UTIL_FLOAT

template<int dummy = 0>
SimdFloat gmx_simdcall loadDualHsimd(const float gmx_unused* m0, const float gmx_unused* m1)
{
    static_assert(((void)dummy, GMX_SIMD_HAVE_HSIMD_UTIL_FLOAT),
                  "Hsimd function used in unconditionally compiled code");
    return {};
}

template<int dummy = 0>
SimdFloat gmx_simdcall loadDuplicateHsimd(const float gmx_unused* m)
{
    static_assert(((void)dummy, GMX_SIMD_HAVE_HSIMD_UTIL_FLOAT),
                  "Hsimd function used in unconditionally compiled code");
    return {};
}

template<int dummy = 0>
SimdFloat gmx_simdcall loadU1DualHsimd(const float gmx_unused* m)
{
    static_assert(((void)dummy, 0), "Hsimd function used without Hsimd support");
    return {};
}

template<int dummy = 0>
void gmx_simdcall storeDualHsimd(float gmx_unused* m0, float gmx_unused* m1, SimdFloat gmx_unused a)
{
    static_assert(((void)dummy, 0), "Hsimd function used without Hsimd support");
}

template<int dummy = 0>
void gmx_simdcall incrDualHsimd(float gmx_unused* m0, float gmx_unused* m1, SimdFloat gmx_unused a)
{
    static_assert(((void)dummy, 0), "Hsimd function used without Hsimd support");
}

template<int dummy = 0>
void gmx_simdcall decr3Hsimd(float gmx_unused*    m,
                             SimdFloat gmx_unused a0,
                             SimdFloat gmx_unused a1,
                             SimdFloat gmx_unused a2)
{
    static_assert(((void)dummy, 0), "Hsimd function used without Hsimd support");
}

template<int align>
void gmx_simdcall gatherLoadTransposeHsimd(const float gmx_unused* base0,
                                           const float gmx_unused*       base1,
                                           const std::int32_t gmx_unused offset[],
                                           SimdFloat gmx_unused* v0,
                                           SimdFloat gmx_unused* v1)
{
    static_assert(((void)align, 0), "Hsimd function used without Hsimd support");
}

template<int dummy = 0>
float gmx_simdcall reduceIncr4ReturnSumHsimd(float gmx_unused* m, SimdFloat gmx_unused v0, SimdFloat gmx_unused v1)
{
    static_assert(((void)dummy, 0), "Hsimd function used without Hsimd support");
    return 0;
}

#    endif // !GMX_SIMD_HAVE_HSIMD_UTIL_FLOAT
#endif     // GMX_SIMD_HAVE_FLOAT

#if GMX_SIMD_HAVE_DOUBLE
#    if !GMX_SIMD_HAVE_HSIMD_UTIL_DOUBLE

template<int dummy = 0>
SimdDouble gmx_simdcall loadDualHsimd(const double gmx_unused* m0, const double gmx_unused* m1)
{
    static_assert(((void)dummy, 0), "Hsimd function used without Hsimd support");
    return {};
}


template<int dummy = 0>
SimdDouble gmx_simdcall loadDuplicateHsimd(const double gmx_unused* m)
{
    static_assert(((void)dummy, 0), "Hsimd function used without Hsimd support");
    return {};
}

template<int dummy = 0>
SimdDouble gmx_simdcall loadU1DualHsimd(const double gmx_unused* m)
{
    static_assert(((void)dummy, 0), "Hsimd function used without Hsimd support");
    return {};
}

template<int dummy = 0>
void gmx_simdcall storeDualHsimd(double gmx_unused* m0, double gmx_unused* m1, SimdDouble gmx_unused a)
{
    static_assert(((void)dummy, 0), "Hsimd function used without Hsimd support");
}

template<int dummy = 0>
void gmx_simdcall incrDualHsimd(double gmx_unused* m0, double gmx_unused* m1, SimdDouble gmx_unused a)
{
    static_assert(((void)dummy, 0), "Hsimd function used without Hsimd support");
}

template<int dummy = 0>
void gmx_simdcall decr3Hsimd(double gmx_unused*    m,
                             SimdDouble gmx_unused a0,
                             SimdDouble gmx_unused a1,
                             SimdDouble gmx_unused a2)
{
    static_assert(((void)dummy, 0), "Hsimd function used without Hsimd support");
}

template<int align>
void gmx_simdcall gatherLoadTransposeHsimd(const double gmx_unused* base0,
                                           const double gmx_unused*      base1,
                                           const std::int32_t gmx_unused offset[],
                                           SimdDouble gmx_unused* v0,
                                           SimdDouble gmx_unused* v1)
{
    static_assert(((void)align, 0), "Hsimd function used without Hsimd support");
}

template<int dummy = 0>
double gmx_simdcall reduceIncr4ReturnSumHsimd(double gmx_unused*    m,
                                              SimdDouble gmx_unused v0,
                                              SimdDouble gmx_unused v1)
{
    static_assert(((void)dummy, 0), "Hsimd function used without Hsimd support");
    return 0;
}

#    endif // !GMX_SIMD_HAVE_HSIMD_UTIL_DOUBLE
#endif     // GMX_SIMD_HAVE_DOUBLE

} // namespace gmx

#endif // GMX_SIMD_HSIMD_DECLARATIONS_H
