/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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
 * \brief Declares SimdArrayRef
 *
 * \author Roland Schulz <roland.schulz@intel.com>
 * \inlibraryapi
 * \ingroup module_simd
 */
#ifndef GMX_SIMD_SIMD_MEMORY_H
#define GMX_SIMD_SIMD_MEMORY_H

#include "gromacs/simd/simd.h"
#include "gromacs/utility/gmxassert.h"

namespace gmx
{


//Possible extension:
// - Make this more complete like ArrayRef (possible by sharing code with ArrayRef)
// - Add a unalinged version

template<typename T>
class SimdArrayRef
{
    static constexpr int simdWidth = SimdTraits<T>::width;

    public:
        /* TODO: for child patch:
           SimdArrayRef(std::vector<T, gmx::Allocator<T, gmx::AlignedAllocationPolicy> > &v)
            : begin_((!v.empty()) ? &v[0] : nullptr),
              end_((!v.empty()) ? &v[0] + v.size() : nullptr) {}
         */

        /*! \brief
         * Constructs a reference to a particular range.
         *
         * \param[in] begin  Pointer to the beginning of a range.
         * \param[in] end    Pointer to the end of a range.
         *
         * Passed pointers must remain valid for the lifetime of this object.
         * Begin pointer as to be alinged according to SIMD requirement.
         * The size (end-begin) has to be either multiple of the SIMD width,
         * or sufficient padding after the end has to be guranteed so that
         * load/stores with full SIMD width is legal for the last element.
         *
         */
        SimdArrayRef(T* begin, T* end)
            : begin_(begin), end_(end)
        {
            GMX_ASSERT(end >= begin, "Invalid range");
        }

        //! Returns the size of the container.
        size_t size() const { return (end_-begin_+simdWidth-1)/simdWidth; }

        //! Access container element.
        SimdReference<T> operator[](size_t n)
        {
            return load(begin_+n*simdWidth);
        }
    private:
        T* const     begin_;
        T* const     end_;
};

//! Construct SimdArrayRef from a pointer and size
//! \related SimdArrayRef
template <typename T>
SimdArrayRef<T> simdArrayRefFromArray(T *begin, size_t size)
{
    return SimdArrayRef<T>(begin, begin+size);
}


} // namespace gmx

#endif
