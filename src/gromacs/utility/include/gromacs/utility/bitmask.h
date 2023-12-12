/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2014- The GROMACS Authors
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
/*! \libinternal \file
 * \brief
 * Declares gmx_bitmask_t and associated functions
 *
 * \author Roland Schulz <roland@utk.edu>
 * \inlibraryapi
 * \ingroup module_utility
 */

#ifndef GMX_MDLIB_BITMASK_H
#define GMX_MDLIB_BITMASK_H

#include "config.h" /* for GMX_MAX_OPENMP_THREADS */

#include <cstring>

#include <algorithm>
#include <array>

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/stringutil.h"

/*! \brief Size of bitmask. Has to be 32 or multiple of 64. */
#ifndef BITMASK_SIZE
#    define BITMASK_SIZE GMX_OPENMP_MAX_THREADS
#endif

#if BITMASK_SIZE != 32 && BITMASK_SIZE % 64 != 0
#    error BITMASK_SIZE has to be 32 or a multiple of 64.
#endif

#if BITMASK_SIZE <= 64 || defined DOXYGEN
#    if BITMASK_SIZE == 32
typedef uint32_t gmx_bitmask_t;
#    else
typedef uint64_t gmx_bitmask_t; /**< bitmask type */
#    endif

/*! \brief Initialize all bits to 0 */
inline static void bitmask_clear(gmx_bitmask_t* m)
{
    *m = 0;
}

/*! \brief Set bit at position b to 1. */
inline static void bitmask_set_bit(gmx_bitmask_t* m, int b)
{
    *m |= gmx_bitmask_t(1) << b;
}

/*! \brief Initialize all bits: bit b to 1, others to 0 */
inline static void bitmask_init_bit(gmx_bitmask_t* m, int b)
{
    *m = gmx_bitmask_t(1) << b;
}

/*! \brief Initialize all bits: all bits below b to 1, others to 0 */
inline static void bitmask_init_low_bits(gmx_bitmask_t* m, int b)
{
    *m = (gmx_bitmask_t(1) << b) - 1;
}

/*! \brief Test if bit b is set */
inline static bool bitmask_is_set(gmx_bitmask_t m, int b)
{
    return (m & (gmx_bitmask_t(1) << b)) != 0;
}

/*! \brief Test if both bitmasks have no common bits enabled */
inline static bool bitmask_is_disjoint(gmx_bitmask_t a, gmx_bitmask_t b)
{
    return (a & b) == 0U;
}

/*! \brief Test if both bitmasks are equal */
inline static bool bitmask_is_equal(gmx_bitmask_t a, gmx_bitmask_t b)
{
    return a == b;
}

/*! \brief Test if bitmask has no enabled bits */
inline static bool bitmask_is_zero(gmx_bitmask_t m)
{
    return m == 0U;
}

/*! \brief Set all bits enabled in either mask and write into a */
inline static void bitmask_union(gmx_bitmask_t* a, gmx_bitmask_t b)
{
    *a |= b;
}
#else
#    define BITMASK_ALEN (BITMASK_SIZE / 64)
using gmx_bitmask_t = std::array<uint64_t, BITMASK_ALEN>;

inline static void bitmask_clear(gmx_bitmask_t* m)
{
    m->fill(0);
}

inline static void bitmask_set_bit(gmx_bitmask_t* m, int b)
{
    (*m)[b / 64] |= (static_cast<uint64_t>(1) << (b % 64));
}

inline static void bitmask_init_bit(gmx_bitmask_t* m, int b)
{
    bitmask_clear(m);
    (*m)[b / 64] = (static_cast<uint64_t>(1) << (b % 64));
}

inline static void bitmask_init_low_bits(gmx_bitmask_t* m, int b)
{
    std::memset(m->data(), 255, b / 64 * 8);
    (*m)[b / 64] = (static_cast<uint64_t>(1) << (b % 64)) - 1;
    std::memset(m->data() + (b / 64 + 1), 0, (BITMASK_ALEN - b / 64 - 1) * 8);
}

inline static bool bitmask_is_set(gmx_bitmask_t m, int b)
{
    return (m[b / 64] & (static_cast<uint64_t>(1) << (b % 64))) != 0;
}

inline static bool bitmask_is_disjoint(gmx_bitmask_t a, gmx_bitmask_t b)
{
    bool r = true;
    for (int i = 0; i < BITMASK_ALEN; i++)
    {
        r = r && ((a[i] & b[i]) == 0U);
    }
    return r;
}

inline static bool bitmask_is_equal(gmx_bitmask_t a, gmx_bitmask_t b)
{
    bool r = true;
    for (int i = 0; i < BITMASK_ALEN; i++)
    {
        r = r && (a[i] == b[i]);
    }
    return r;
}

inline static bool bitmask_is_zero(gmx_bitmask_t m)
{
    bool r = true;
    for (int i = 0; i < BITMASK_ALEN; i++)
    {
        r = r && (m[i] == 0U);
    }
    return r;
}

inline static void bitmask_union(gmx_bitmask_t* a, gmx_bitmask_t b)
{
    for (int i = 0; i < BITMASK_ALEN; i++)
    {
        (*a)[i] |= b[i];
    }
}
#endif
// In bitmask.h because only current use is for bitmask.

//! Convert uint32_t to hex string
inline static std::string to_hex_string(uint32_t m)
{
    return gmx::formatString("%08" PRIx32, m);
}
//! Convert uint64_t to hex string
inline static std::string to_hex_string(uint64_t m)
{
    return gmx::formatString("%016" PRIx64, m);
}
//! Convert container of integers to hex string
template<typename C>
inline static std::string to_hex_string(C m)
{
    std::string ret;
    for (auto it = m.rbegin(); it < m.rend(); ++it)
    {
        ret += to_hex_string(*it);
    }
    return ret;
}

#endif
