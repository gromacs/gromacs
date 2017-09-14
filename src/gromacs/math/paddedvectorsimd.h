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
/*! \file
 * \brief
 * Declares gmx::PaddedVector
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \inpublicapi
 * \ingroup module_math
 */
// TODO eventually this will replace paddedvector.h
#ifndef GMX_MATH_PADDEDVECTORSIMD_H
#define GMX_MATH_PADDEDVECTORSIMD_H

#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/simd/simd.h"
#include "gromacs/utility/alignedallocator.h"
#include "gromacs/utility/arrayref.h"

namespace gmx
{

namespace
{

}   // namespace

/*! \class vector
 *  \brief vector interface
 *
 * \tparam T the type of objects within the container
 * \tparam Alloc the allocator used. Can be a standard compliant allocator or an advanced allocator as pt::malloc_allocator
 *
 * This is the main vector interface. It is compliant with the C++11 standard of std::vector.
 * Thus, as a reference, the documentation available here
 * (http://en.cppreference.com/w/cpp/container/vector) works for this object.
 */
template <class BaseType, size_t simdWidthOfBaseType, size_t packSizeOfBaseType = 1, class Allocator = std::allocator<BaseType> >
class PaddedVector
{
    public:
        typedef BaseType value_type;
        typedef typename Allocator::template rebind<BaseType>::other allocator_type;
        typedef size_t size_type;
        typedef std::ptrdiff_t difference_type;
        typedef value_type &reference;
        typedef value_type const &const_reference;
        // TODO enable_if only if BaseType is POD (or scalar, if that is possible)
        static constexpr bool is_pod = std::is_pod<value_type>::value;

    private:
        typedef std::vector<BaseType, Allocator> storage_type;
    public:
        typedef typename storage_type::pointer       pointer;
        typedef typename storage_type::const_pointer const_pointer;

        typedef typename storage_type::iterator        iterator;
        typedef typename storage_type::const_iterator  const_iterator;
        typedef std::reverse_iterator<iterator>        reverse_iterator;
        typedef std::reverse_iterator<const_iterator>  const_reverse_iterator;

        typedef ArrayRef<BaseType> array_ref;
        typedef ConstArrayRef<BaseType> const_array_ref;

    private:
        /*! \brief Helper function for computing the padded size.
         *
         * Care is needed to compute in bytes and report in the
         * units of size. from \p size given \p the padding width. */
        std::size_t computePaddedSize(const std::size_t numElements)
        {
            constexpr std::size_t simdWidthOfBaseTypeLessOne = simdWidthOfBaseType - 1;
            return ((numElements * packSizeOfBaseType + simdWidthOfBaseTypeLessOne) / simdWidthOfBaseType) * simdWidthOfBaseType;
        }

        /*! \brief Helper function for updating unpadded_end_ after
         * inserting/emplacing/appending/prepending \c numNewElements elements.
         *
         * \param[in]  oldBegin        The value of storage_.begin() before the pointer was
         *                             potentially invalidated by adding the element.
         * \param[in]  numNewElements  The number of new elements in the container.
         */
        void updateUnpaddedEnd(const pointer oldBegin, const size_type numNewElements)
        {
            if (begin() != oldBegin)
            {
                auto old_unpadded_size = unpadded_end_ - oldBegin;
                unpadded_end_ = begin() + old_unpadded_size;
            }
            unpadded_end_ += numNewElements;
        }

    public:
        PaddedVector() :
            storage_(),
            unpadded_end_(data())
        {}

        // TODO This should also be specialized by allocator, but
        // storage_ doesn't have such a constructor before
        // C++14. Resolve.
        explicit PaddedVector(size_type count) :
            storage_(computePaddedSize(count)),
            unpadded_end_(data() + count)
        {}

        explicit PaddedVector(allocator_type const &alloc) :
            storage_(alloc),
            unpadded_end_(data())
        {}

        explicit PaddedVector(PaddedVector const &o) :
            storage_(o.storage_),
            unpadded_end_(begin() + o.unpadded_size())
        {}

        explicit PaddedVector(PaddedVector &&o) :
            storage_(std::move(o.storage_)),
            unpadded_end_(o.unpadded_end())
        {}

        explicit PaddedVector(std::initializer_list<value_type> const &il) :
            storage_(0),
            unpadded_end_(storage_.data())
        {
            storage_type temp(il);
            auto         temp_initial_size = temp.size();
            temp.resize(computePaddedSize(temp_initial_size));
            storage_.swap(temp);
            unpadded_end_ = data() + temp_initial_size;
        }

        void reserve(const size_type new_extent)
        {
            auto unpaddedSize = unpadded_end() - begin();
            /* v.reserve(13) should allocate enough memory so that
               v.resize(13) does not reallocate. This means that the
               new extent should be large enough for the padded
               storage for a vector whose size is new_extent. */
            auto new_padded_extent = computePaddedSize(new_extent);
            storage_.reserve(new_padded_extent);
            unpadded_end_ = data() + unpaddedSize;
        }

        void resize(const size_type newSize)
        {
            auto new_padded_size = computePaddedSize(newSize);
            storage_.resize(new_padded_size);
            unpadded_end_ = data() + newSize;
        }

        void resize(const size_type newSize, value_type const &v)
        {
            auto new_padded_size = computePaddedSize(newSize);
            storage_.resize(new_padded_size, v);
            unpadded_end_ = data() + newSize;
        }

        size_type size() const { return storage_.size(); }

        size_type capacity() const { return storage_.capacity(); }
        bool empty() const { return size() == 0; }

        size_type unpadded_size() const { return unpadded_end() - begin(); }

        // TODO
        /*
           template <class InputIterator>
           void assign(InputIterator first, InputIterator last)
           {
           // TODO
           }

           void assign(std::initializer_list<value_type> const& li)
           {
            assign(li.begin(), li.end());
           }

           void assign(size_type const n, value_type const& v)
           {
           }
         */
    public:
        void push_back(value_type const &v)
        {
            auto old_begin        = begin();
            auto emplaced_element = storage_.push_back(v);
            updateUnpaddedEnd(old_begin, 1);
            return emplaced_element;
        }

        void pop_back()
        {
            storage_.pop_back();
            unpadded_end_--;
        }

        void swap(PaddedVector &x)
        {
            std::swap(storage_, x.storage_);
            std::swap(unpadded_end_, x.unpadded_end_);
        }

        void clear()
        {
            storage_.clear();
            unpadded_end_ = data();
        }

        template <class ... Args>
        void emplace_back(Args && ... args)
        {
            auto old_begin        = begin();
            auto emplaced_element = storage_.emplace_back(std::forward<Args>(args) ...);
            updateUnpaddedEnd(old_begin, 1);
            return emplaced_element;
        }

        template< class InputIt>
        iterator insert(typename std::enable_if<std::is_base_of<std::input_iterator_tag, typename std::iterator_traits<InputIt>::iterator_category>::value, const_iterator>::type pos, InputIt first, InputIt last)
        {
            auto old_begin        = begin();
            auto emplaced_element = storage_.insert(pos, first, last);
            updateUnpaddedEnd(old_begin, last - first);
            return emplaced_element;
        }

        iterator insert(const_iterator pos, size_type const count, value_type const &v)
        {
            auto old_begin        = begin();
            auto emplaced_element = storage_.insert(pos, count, v);
            updateUnpaddedEnd(old_begin, count);
            return emplaced_element;
        }

        template<class ... Args>
        iterator emplace(const_iterator pos, Args && ... args)
        {
            auto old_begin        = begin();
            auto emplaced_element = storage_.emplace(pos, std::forward<Args>(args) ...);
            updateUnpaddedEnd(old_begin, 1);
            return emplaced_element;
        }

        iterator insert(const_iterator pos, value_type &&v)
        {
            return emplace(pos, std::move(v));
        }

        iterator insert(const_iterator pos, value_type const &v)
        {
            return insert(pos, 1, v);
        }

        iterator insert(const_iterator pos, std::initializer_list<value_type> const &il)
        {
            return insert(pos, il.begin(), il.end());
        }

        // TODO need to decide upon, and implement correct error handling for out-of-range access
        reference at(const size_type i)
        {
            return storage_.at(i);
        }

        const_reference at(const size_type i) const
        {
            return storage_.at(i);
        }

        reference       operator[](const size_type i) { return at(i); }
        const_reference operator[](const size_type i) const { return at(i); }

        reference       front()       { return *begin(); }
        const_reference front() const { return *begin(); }

        reference       back()       { return *end(); }
        const_reference back() const { return *end(); }

        pointer       data()       noexcept { return storage_.data(); }
        const_pointer data() const noexcept { return storage_.data(); }

        iterator       begin()        { return storage_.begin(); }
        iterator       end()          { return storage_.end(); }
        iterator       unpadded_end() { return iterator(unpadded_end_); }

        const_iterator cbegin()        { return const_iterator(begin()); }
        const_iterator cend()          { return const_iterator(end()); }
        const_iterator unpadded_cend() { return const_iterator(unpadded_end_); }

        const_iterator begin()        const { return storage_.begin(); }
        const_iterator end()          const { return storage_.end(); }
        const_iterator unpadded_end() const { return const_iterator(unpadded_end_); }

        const_iterator cbegin()        const { return const_iterator(begin()); }
        const_iterator cend()          const { return const_iterator(end()); }
        const_iterator unpadded_cend() const { return const_iterator(unpadded_end_); }

        reverse_iterator rbegin()        { return reverse_iterator(end()); }
        reverse_iterator rend()          { return reverse_iterator(begin()); }
        reverse_iterator unpadded_rend() { return reverse_iterator(unpadded_end_); }

        const_reverse_iterator crbegin()        { return const_reverse_iterator(end()); }
        const_reverse_iterator crend()          { return const_reverse_iterator(begin()); }
        const_reverse_iterator unpadded_crend() { return const_reverse_iterator(unpadded_end_); }

        const_reverse_iterator rbegin()        const { return const_reverse_iterator(end()); }
        const_reverse_iterator rend()          const { return const_reverse_iterator(begin()); }
        const_reverse_iterator unpadded_rend() const { return const_reverse_iterator(unpadded_end_); }

        const_reverse_iterator crbegin()        const { return const_reverse_iterator(end()); }
        const_reverse_iterator crend()          const { return const_reverse_iterator(begin()); }
        const_reverse_iterator unpadded_crend() const { return const_reverse_iterator(unpadded_end_); }

        array_ref unpaddedArrayRef()
        {
            return arrayRefFromPointers<BaseType>(data(), unpadded_end_);
        }

        const_array_ref unpaddedConstArrayRef() const
        {
            return constArrayRefFromPointers<BaseType>(data(), unpadded_end_);
        }
        /* TODO where suitable (enforced with std::enable_if?), make
         * other getter functions that return ArrayRef<RVec> and
         * ConstArray<RVec> views of the unpadded data. */
    public:
        inline bool operator==(PaddedVector const &o) const
        {
            return storage_type::equals(begin(), size(), o.begin(), o.size());
        }

        inline bool operator<(PaddedVector const &o) const
        {
            return std::lexicographical_compare(begin(), end(), o.begin(), o.end());
        }

    public:
        PaddedVector &operator=(PaddedVector const &o)
        {
            if (&o != this)
            {
                storage_ = o.storage_;
            }
            return *this;
        }

        PaddedVector &operator=(PaddedVector &&o)
        {
            if (&o != this)
            {
                storage_ = std::move(o.storage_);
            }
            return *this;
        }

    private:
        storage_type storage_;
        pointer      unpadded_end_;
};

template <class BaseType>
using PaddedSimdRVecVector = PaddedVector<real, GMX_SIMD_REAL_WIDTH, 3, AlignedAllocator<RVec> >;

//! \copydoc ArrayRef::fromVector()
//! \related ArrayRef
template <class BaseType, size_t simdWidthOfBaseType, size_t packSizeOfBaseType = 1, class Allocator = std::allocator<BaseType> >
ArrayRef<BaseType> arrayRefFromVector(typename PaddedVector<BaseType, simdWidthOfBaseType, packSizeOfBaseType, Allocator>::iterator begin,
                                      typename PaddedVector<BaseType, simdWidthOfBaseType, packSizeOfBaseType, Allocator>::iterator end)
{
    return ArrayRef<BaseType>::fromPointers(begin, end);
}

//! \copydoc ConstArrayRef::fromVector()
//! \related ConstArrayRef
template <class BaseType, size_t simdWidthOfBaseType, size_t packSizeOfBaseType = 1, class Allocator = std::allocator<BaseType> >
ConstArrayRef<BaseType> constArrayRefFromVector(typename PaddedVector<BaseType, simdWidthOfBaseType, packSizeOfBaseType, Allocator>::const_iterator begin,
                                                typename PaddedVector<BaseType, simdWidthOfBaseType, packSizeOfBaseType, Allocator>::const_iterator end)
{
    return ConstArrayRef<BaseType>::fromPointers(begin, end);
}

} // namespace gmx

#endif
