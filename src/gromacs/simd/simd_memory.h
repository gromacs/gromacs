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

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/gmxassert.h"

namespace gmx
{
namespace internal
{

/*! \internal
 * \brief STL-like container for aligned SIMD type. Used as ArrayRef<SimdReal>.
 *
 * Should provide the same interface as ArrayRef. Any missing functions should be
 * added as needed. The pointer type (used e.g. for initialization) is a real
 * pointer. The reference type (used e.g. for operator[] and iterator dereference)
 * is SimdReference which executes the aligned load/store as is appropriate. For
 * both iterator and element access, the access happens in blocks of Simd width.
 * Meaning that a[1] refers to the 2nd SIMD vector and thus reals 8-15 for 8-wide
 * SIMD. The starting address has to be aligned and the length has to be multiple
 * of the SIMD width.
 */
template<typename T>
class SimdArrayRef
{
    private:
        typedef typename std::remove_const<T>::type non_const_value_type;
        static constexpr int simdWidth = SimdTraits<T>::width;

    public:
        //! Type for representing size of the container.
        typedef size_t    size_type;
        //! Type for representing difference between two container indices.
        using difference_type   = std::ptrdiff_t;
        //! Type of values stored in the container.
        using value_type        = typename SimdTraits<T>::type;
        //! Pointer to a container element.
        using pointer           = T*;
        //! Reference to a container element.
        using reference         = SimdReference<T>;

        class iterator
        {
            public:
                using difference_type   = difference_type;
                using value_type        = value_type;
                using pointer           = pointer;
                using reference         = reference;
                using iterator_category = std::random_access_iterator_tag;

                explicit iterator(pointer p = 0) : p_(p)
                {
                    GMX_ASSERT(reinterpret_cast<size_type>(p)%(simdWidth*sizeof(T)) == 0,
                               "Trying to create aligned iterator for non aligned address.");
                }
                iterator &operator++()
                {
                    p_ += simdWidth;
                    return *this;
                }
                iterator operator++(int)
                {
                    iterator retval = *this;
                    ++(*this);
                    return retval;
                }
                iterator &operator--()
                {
                    p_ -= simdWidth;
                    return *this;
                }
                iterator operator--(int)
                {
                    iterator retval = *this;
                    --(*this);
                    return retval;
                }
                iterator &operator+=(difference_type d)
                {
                    p_ += simdWidth * d;
                    return *this;
                }
                iterator &operator-=(difference_type d)
                {
                    p_ -= simdWidth * d;
                    return *this;
                }
                iterator operator+(difference_type d) { return iterator(p_ + simdWidth*d); }
                iterator operator-(difference_type d) { return iterator(p_ - simdWidth*d); }

                bool operator==(iterator other) const { return p_ == other.p_; }
                bool operator!=(iterator other) const { return p_ != other.p_; }
                bool operator< (iterator other) const { return p_ <  other.p_; }
                bool operator> (iterator other) const { return p_ >  other.p_; }
                bool operator<=(iterator other) const { return p_ <= other.p_; }
                bool operator>=(iterator other) const { return p_ >= other.p_; }

                reference operator*() const { return load(p_); }

                operator pointer() const { return p_; }
            private:
                pointer p_;
        };

        //! Standard reverse iterator.
        typedef std::reverse_iterator<iterator>       reverse_iterator;

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
         * Begin pointer has to be aligned according to SIMD requirement.
         * The size (end-begin) has to be either multiple of the SIMD width,
         * or sufficient padding after the end has to be guaranteed so that
         * load/stores with full SIMD width is legal for the last element.
         *
         */
        SimdArrayRef(T* begin, T* end)
            : begin_(begin), end_(end)
        {
            GMX_ASSERT(end >= begin, "Invalid range");
            GMX_ASSERT(reinterpret_cast<size_type>(begin)%(simdWidth*sizeof(T)) == 0,
                       "Aligned ArrayRef requires aligned starting address");
            GMX_ASSERT(reinterpret_cast<size_type>(end)%(simdWidth*sizeof(T)) == 0,
                       "Size of ArrayRef needs to be divisible by type size");
        }

        /*! \brief
         * Constructs a reference to const data from a reference to non-const data.
         *
         * Constructs a ArrayRef<const T> from a ArrayRef<T>.
         */
        template<typename = T> //Otherwise useless template argument
                               //to avoid this being used as copy constructor
        SimdArrayRef(const SimdArrayRef<non_const_value_type> &o) :
            begin_(o.begin()), end_(o.end()) {}

        //! Returns the size of the container.
        size_type size() const { return (end_-begin_)/simdWidth; }
        //! Returns an iterator to the beginning of the container.
        iterator begin() const { return iterator(begin_); }
        //! Returns an iterator to the end of the container.
        iterator end() const { return iterator(end_); }

        //! Access container element.
        reference operator[](size_type n)
        {
            return load(begin_+n*simdWidth);
        }
    private:
        T* const     begin_;
        T* const     end_;
};

}   //namespace internal

/* Specialize ArraryRef<SimdReal>
 * So far only an aligned version is implemented. The constructor verifies that
 * a ArrayRef<SimdReal> is constructed for only a properly aligned data.
 */
#if GMX_SIMD_HAVE_FLOAT
template<>
class ArrayRef<SimdFloat>  : public internal::SimdArrayRef<float>
{
    using Base = internal::SimdArrayRef<float>;
    using Base::Base;
};
template<>
class ArrayRef<const SimdFloat>  : public internal::SimdArrayRef<const float>
{
    using Base = internal::SimdArrayRef<const float>;
    using Base::Base;
};
#endif
#if GMX_SIMD_HAVE_DOUBLE
template<>
class ArrayRef<SimdDouble> : public internal::SimdArrayRef<double>
{
    using Base = internal::SimdArrayRef<double>;
    using Base::Base;
};
template<>
class ArrayRef<const SimdDouble> : public internal::SimdArrayRef<const double>
{
    using Base = internal::SimdArrayRef<const double>;
    using Base::Base;
};
#endif

} // namespace gmx

#endif
