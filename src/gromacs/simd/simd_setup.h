/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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

    \brief   helper class that determines parameters and types for setting up SIMD calculations

    \author R. Thomas Ullmann <tullman@gwdg.de>

    \copyright GROMACS license

    \date May 2016

    \inlibraryapi
    \ingroup module_simd
 */
#ifndef GMX_SIMD_SIMD_SETUP_H
#define GMX_SIMD_SIMD_SETUP_H

#include <memory>
#include <type_traits>

#include "gromacs/simd/simd.h"
#include "gromacs/utility/alignedallocator.h"

namespace gmx
{

/*! \brief  setup of default memory allocation and data types for SIMD calculations with the GROMACS SIMD module

    \tparam   T   scalar data type

    \libinternal
    \ingroup module_simd
 */
template<typename T>
class SimdSetup
{
    public:
        //! scalar float type for which SIMD data types and parameters shall be determined
        typedef   T    float_type;

#if GMX_SIMD

        //! true if smaller blocks of data can be processed separately with SimdWidth 4
        static constexpr bool doSimd_       = ((static_cast<bool>(GMX_SIMD != 0) &&
                                                ((std::is_same<float_type, float>::value  && static_cast<bool>(GMX_SIMD_HAVE_FLOAT != 0) ) ||
                                                 (std::is_same<float_type, double>::value && static_cast<bool>(GMX_SIMD_HAVE_DOUBLE != 0)))) ? true : false);
        //! the SIMD width, i.e., the number of SIMD vector elements processed in a SIMD operation
        static constexpr size_t simdWidth_  = ((doSimd_ && std::is_same<float_type, float>::value  && static_cast<bool>(GMX_SIMD_HAVE_FLOAT != 0)) ? GMX_SIMD_FLOAT_WIDTH :
                                               (doSimd_ && std::is_same<float_type, double>::value && static_cast<bool>(GMX_SIMD_HAVE_DOUBLE != 0) ? GMX_SIMD_DOUBLE_WIDTH : 1u));
        //! the SIMD4 width, should be 4
        static constexpr size_t simd4Width_ = (doSimd_ &&
                                               ((std::is_same<float_type, float>::value  && static_cast<bool>(GMX_SIMD4_HAVE_FLOAT != 0) ) ||
                                                (std::is_same<float_type, double>::value && static_cast<bool>(GMX_SIMD4_HAVE_DOUBLE != 0))) ? GMX_SIMD4_WIDTH : 1u);
        //! true if smaller blocks of data can be processed separately with SimdWidth 4
        static constexpr bool doSimd4_      = ((simd4Width_ < simdWidth_ && (simdWidth_ % simd4Width_) == 0) &&
                                               ((std::is_same<float_type, float>::value  && static_cast<bool>(GMX_SIMD4_HAVE_FLOAT != 0) ) ||
                                                (std::is_same<float_type, double>::value && static_cast<bool>(GMX_SIMD4_HAVE_DOUBLE != 0))) ? true : false);
        //! default alignment in byte
        static constexpr size_t alignment_  = alignof(float_type) * simdWidth_;
        //! default padding in byte
        static constexpr size_t padding_    = 0u;

        //! choose the SimdFloat or SimdDouble to match float_type, SIMD vectors with implementation dependent internal storage
        typedef typename std::conditional<std::is_same<float_type, float>::value, gmx::SimdFloat,  gmx::SimdDouble>::type   simd_type;
#if GMX_SIMD4_HAVE_FLOAT && GMX_SIMD4_HAVE_DOUBLE
        //! choose Simd4Float or Simd4Double if available and if the SimdWidth is not 4 or smaller anyway (availability already checked while assigning doSimd4_
        typedef typename std::conditional<doSimd4_,
                                          typename std::conditional<std::is_same<float_type, float>::value, gmx::Simd4Float, gmx::Simd4Double>::type,
                                          float>::type simd4_type;
#elif GMX_SIMD4_HAVE_FLOAT
        //! choose Simd4Float or Simd4Double if available and if the SimdWidth is not 4 or smaller anyway (availability already checked while assigning doSimd4_
        typedef typename std::conditional<doSimd4_,
                                          typename std::conditional<std::is_same<float_type, float>::value, gmx::Simd4Float, float>::type,
                                          float>::type simd4_type;
#else
        typedef                           float        simd4_type;
#endif
#else   // !GMX_SIMD
        //! true if SIMD computation is available
        static constexpr bool   doSimd_       = false;
        //! the SIMD width, i.e., the number of SIMD vector elements processed in a SIMD operation
        static constexpr size_t simdWidth_  = 1u;
        //! the SIMD4 width, should be 4
        static constexpr size_t simd4Width_ = 1u;
        //! true if smaller blocks of data can be processed separately with SimdWidth 4
        static constexpr bool   doSimd4_      = false;
        //! default alignment in byte
        static constexpr size_t alignment_  = alignof(float_type);
        //! default padding in byte
        static constexpr size_t padding_    = 0u;
#endif
        //! default allocator type
        typedef typename std::conditional<doSimd_, gmx::AlignedAllocator<float_type>, std::allocator<float_type> >::type allocator_type;
};


} // end namespace gmx

#endif
