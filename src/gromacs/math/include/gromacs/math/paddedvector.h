/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2016- The GROMACS Authors
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
/*! \file
 * \brief
 * Declares gmx::PaddedRVecVector
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \inpublicapi
 * \ingroup module_math
 */
#ifndef GMX_MATH_PADDEDVECTOR_H
#define GMX_MATH_PADDEDVECTOR_H

#include <algorithm>
#include <utility>
#include <vector>

#include "gromacs/math/arrayrefwithpadding.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/alignedallocator.h"
#include "gromacs/utility/allocator.h"
#include "gromacs/utility/arrayref.h"

namespace gmx
{

namespace detail
{

/*! \brief Traits classes for handling padding for types used with PaddedVector
 *
 * Only the base types of the SIMD module are supported for
 * PaddedVector, because the purpose of the padding is to permit
 * SIMD-width operations from the SIMD module.
 *
 * \todo Consider explicitly tying these types to the SimdTrait
 * types. That would require depending on the SIMD module, or
 * extracting the traits from it. This would also permit
 * maxSimdWidthOfBaseType to be set more efficiently, e.g. as a
 * metaprogramming max over the maximum width from different
 * implementations.
 */
template<typename T>
struct PaddingTraits
{
};

template<>
struct PaddingTraits<int32_t>
{
    using SimdBaseType                          = int32_t;
    static constexpr int maxSimdWidthOfBaseType = 16;
};

template<>
struct PaddingTraits<float>
{
    using SimdBaseType                          = float;
    static constexpr int maxSimdWidthOfBaseType = GMX_FLOAT_MAX_SIMD_WIDTH;
};

template<>
struct PaddingTraits<double>
{
    using SimdBaseType                          = double;
    static constexpr int maxSimdWidthOfBaseType = GMX_DOUBLE_MAX_SIMD_WIDTH;
};

template<>
struct PaddingTraits<BasicVector<float>>
{
    using SimdBaseType                          = float;
    static constexpr int maxSimdWidthOfBaseType = GMX_FLOAT_MAX_SIMD_WIDTH;
};

template<>
struct PaddingTraits<BasicVector<double>>
{
    using SimdBaseType                          = double;
    static constexpr int maxSimdWidthOfBaseType = GMX_DOUBLE_MAX_SIMD_WIDTH;
};

/*! \brief Returns the allocation size for PaddedVector that contains
 * \c numElements elements plus padding for SIMD operations.
 *
 * \param[in] numElements  The number of T elements for which data will be stored.
 * \returns                The number of T elements that must be allocated
 *                         (ie >= numElements).
 */
template<typename T>
Index computePaddedSize(Index numElements)
{
    // We don't need padding if there is no access.
    if (numElements == 0)
    {
        return 0;
    }

    // We sometimes load a whole extra element when doing 4-wide SIMD
    // operations (which might e.g. be an RVec) so we need to pad for
    // that.
    Index simdScatterAccessSize = numElements + 1;

    // For SIMD updates based on RVec, we might load starting from
    // the last RVec element, so that sets the minimum extent of the
    // padding. That extent must take the initialized allocation up to
    // the SIMD width of the base type multiplied by the width of T in
    // that base type. But since storage_ contains RVec, we only have
    // to tell it the number of elements, which means to round up to
    // the next SIMD width.
    //
    // We don't want a dependence on the SIMD module for the actual
    // SIMD width of the base type, so we use maximum for the base
    // type via the traits. A little extra padding won't really hurt.
    constexpr int maxSimdWidth = PaddingTraits<T>::maxSimdWidthOfBaseType;
    Index simdFlatAccessSize   = (numElements + (maxSimdWidth - 1)) / maxSimdWidth * maxSimdWidth;

    return std::max(simdScatterAccessSize, simdFlatAccessSize);
}

//! Helper function to insert padding elements for most T.
template<typename T, typename AllocatorType>
inline void insertPaddingElements(std::vector<T, AllocatorType>* v, Index newPaddedSize)
{
    // Ensure the padding region is initialized to zero. There is no
    // way to insert a number of default-initialized elements. So we
    // have to provide a value for those elements, which anyway suits
    // this use case.
    v->insert(v->end(), newPaddedSize - v->size(), 0);
}

//! Specialization of helper function to insert padding elements, used for BasicVector<T>.
template<typename T, typename AllocatorType>
inline void insertPaddingElements(std::vector<BasicVector<T>, AllocatorType>* v, Index newPaddedSize)
{
    // Ensure the padding region is initialized to zero.
    v->insert(v->end(), newPaddedSize - v->size(), BasicVector<T>(0, 0, 0));
}

} // namespace detail

/*! \brief PaddedVector is a container of elements in contiguous
 * storage that allocates extra memory for safe SIMD-style loads for
 * operations used in GROMACS.
 *
 * \tparam T the type of objects within the container
 * \tparam Allocator the allocator used. Can be any standard-compliant
 * allocator, such gmx::Allocator used for alignment and/or pinning.
 *
 * The interface resembles std::vector. However, access
 * intended to include padded elements must be via ArrayRef objects
 * explicitly created to view those elements. Most other aspects of
 * this vector refer to the unpadded view, e.g. iterators, data(),
 * size().
 *
 * The underlying storage is allocated with extra elements, properly
 * initialized, that ensure that any operations accessing the any
 * non-additional element that operate on memory equivalent to a full
 * SIMD lane do so on allocated memory that has been initialized, so
 * that memory traps will not occur, and arithmetic operations will
 * not cause e.g. floating-point exceptions so long as the values in
 * the padded elements are properly managed.
 *
 * Proper initialization is tricker than it would first appear, since
 * we intend this container to be used with scalar and class types
 * (e.g. RVec). Resize and construction operations use "default
 * insertion" which leads to zero initialization for the former, and
 * calling the default constructor for the latter. BasicVector has a
 * default constructor that leaves the elements uninitialized, which
 * is particularly risky for elements only present as padding. Thus
 * the implementation specifically initializes the padded elements to
 * zero, which makes no difference to the scalar template
 * instantiations, and makes the BasicVector ones safer to use.
 *
 * Because the allocator can be configured, the memory allocation can
 * have other attributes such as SIMD alignment or being pinned to
 * physical memory for efficient transfers. The default allocator
 * ensures alignment, but std::allocator also works.
 */
template<typename T, typename Allocator = Allocator<T, AlignedAllocationPolicy>>
class PaddedVector
{
public:
    //! Standard helper types
    //! \{
    using value_type      = T;
    using allocator_type  = Allocator;
    using size_type       = Index;
    using reference       = value_type&;
    using const_reference = const value_type&;
    using storage_type    = std::vector<T, allocator_type>;
    using pointer         = typename storage_type::pointer;
    using const_pointer   = typename storage_type::const_pointer;
    using iterator        = typename storage_type::iterator;
    using const_iterator  = typename storage_type::const_iterator;
    using difference_type = typename storage_type::iterator::difference_type;
    //! \}

    PaddedVector() : storage_(), unpaddedEnd_(begin()) {}
    /*! \brief Constructor that specifies the initial size. */
    explicit PaddedVector(size_type count, const allocator_type& allocator = Allocator()) :
        storage_(count, allocator), unpaddedEnd_(begin() + count)
    {
        // The count elements have been default inserted, and now
        // the padding elements are added
        resizeWithPadding(count);
    }
    /*! \brief Constructor that specifies the initial size and an element to copy. */
    explicit PaddedVector(size_type count, value_type const& v, const allocator_type& allocator = Allocator()) :
        storage_(count, v, allocator), unpaddedEnd_(begin() + count)
    {
        // The count elements have been default inserted, and now
        // the padding elements are added
        resizeWithPadding(count);
    }
    //! Default constructor with allocator
    explicit PaddedVector(allocator_type const& allocator) :
        storage_(allocator), unpaddedEnd_(begin())
    {
    }
    //! Copy constructor
    PaddedVector(PaddedVector const& o) : storage_(o.storage_), unpaddedEnd_(begin() + o.size()) {}
    /*! \brief Move constructor
     *
     * Leaves \c o in a valid state (ie the destructor can be
     * called). */
    PaddedVector(PaddedVector&& o) noexcept :
        storage_(std::exchange(o.storage_, {})), unpaddedEnd_(o.unpaddedEnd_)
    {
    }
    /*! \brief Move constructor using \c alloc for the new vector.
     *
     * Note that \c alloc is another instance of the same allocator
     * type as used for \c PaddedVector. This makes sense e.g. for
     * stateful allocators such as HostAllocator used in
     * PaddedHostVector.
     *
     * Leaves \c o in a valid state (ie. the destructor can be
     * called). */
    PaddedVector(PaddedVector&& o, const Allocator& alloc) noexcept :
        storage_(alloc), unpaddedEnd_(begin())
    {
        if (alloc == o.storage_.get_allocator())
        {
            std::swap(storage_, o.storage_);
            unpaddedEnd_ = o.unpaddedEnd_;
        }
        else
        {
            // If the allocator compares differently, we must
            // reallocate and copy.
            auto unpaddedSize = o.size();
            resizeWithPadding(unpaddedSize);
            std::copy(o.begin(), o.end(), storage_.begin());
            unpaddedEnd_ = begin() + unpaddedSize;
        }
    }
    //! Construct from an initializer list
    PaddedVector(std::initializer_list<value_type> const& il) :
        storage_(il), unpaddedEnd_(storage_.end())
    {
        // We can't choose the padding until we know the size of
        // the normal vector, so we have to make the storage_ and
        // then resize it.
        resizeWithPadding(storage_.size());
    }
    //! Reserve storage for the container to contain newExtent elements, plus the required padding.
    void reserveWithPadding(const size_type newExtent)
    {
        auto unpaddedSize = end() - begin();
        /* v.reserve(13) should allocate enough memory so that
           v.resize(13) does not reallocate. This means that the
           new extent should be large enough for the padded
           storage for a vector whose size is newExtent. */
        auto newPaddedExtent = detail::computePaddedSize<T>(newExtent);
        storage_.reserve(newPaddedExtent);
        unpaddedEnd_ = begin() + unpaddedSize;
    }
    //! Resize the container to contain newSize elements, plus the required padding.
    void resizeWithPadding(const size_type newSize)
    {
        // When the contained type is e.g. a scalar, then the
        // default initialization behaviour is to zero all
        // elements, which is OK, but we have to make sure that it
        // happens for the elements in the padded region when the
        // vector is shrinking.
        auto newPaddedSize = detail::computePaddedSize<T>(newSize);
        // Make sure there is room for padding if we need to grow.
        storage_.reserve(newPaddedSize);
        // Make the unpadded size correct, with any additional
        // elements initialized by the default constructor. It is
        // particularly important to destruct former elements when
        // newSize is smaller than the old size.
        storage_.resize(newSize);
        // Ensure the padding region is zeroed if required.
        detail::insertPaddingElements(&storage_, newPaddedSize);
        unpaddedEnd_ = begin() + newSize;
    }
    //! Return the size of the view without the padding.
    size_type size() const { return end() - begin(); }
    //! Return the container size including the padding.
    size_type paddedSize() const { return storage_.size(); }
    //! Return whether the storage is empty.
    bool empty() const { return storage_.empty(); }
    //! Swap two PaddedVectors
    void swap(PaddedVector& x) noexcept
    {
        std::swap(storage_, x.storage_);
        std::swap(unpaddedEnd_, x.unpaddedEnd_);
    }
    //! Clear the vector, ie. set size to zero and remove padding.
    void clear()
    {
        storage_.clear();
        unpaddedEnd_ = begin();
    }
    //! Iterator getters refer to a view without padding.
    //! \{
    pointer       data() noexcept { return storage_.data(); }
    const_pointer data() const noexcept { return storage_.data(); }

    iterator begin() { return storage_.begin(); }
    iterator end() { return iterator(unpaddedEnd_); }

    const_iterator cbegin() { return const_iterator(begin()); }
    const_iterator cend() { return const_iterator(unpaddedEnd_); }

    const_iterator begin() const { return storage_.begin(); }
    const_iterator end() const { return const_iterator(unpaddedEnd_); }

    const_iterator cbegin() const { return const_iterator(begin()); }
    const_iterator cend() const { return const_iterator(unpaddedEnd_); }
    //! \}
    // TODO should these do bounds checking for the unpadded range? In debug mode?
    //! Indexing operator.
    reference operator[](int i) { return storage_[i]; }
    //! Indexing operator as const.
    const_reference operator[](int i) const { return storage_[i]; }
    //! Returns an ArrayRef of elements that includes the padding region, e.g. for use in SIMD code.
    ArrayRefWithPadding<T> arrayRefWithPadding()
    {
        return ArrayRefWithPadding<T>(data(), data() + size(), data() + paddedSize());
    }
    //! Returns an ArrayRef of const elements that includes the padding region, e.g. for use in SIMD code.
    ArrayRefWithPadding<const T> constArrayRefWithPadding() const
    {
        return ArrayRefWithPadding<const T>(data(), data() + size(), data() + paddedSize());
    }
    /*! \brief Returns an rvec * pointer for containers of RVec, for use with legacy code.
     *
     * \todo Use std::is_same_v when CUDA 11 is a requirement.
     */
    template<typename AlsoT = T, typename = typename std::enable_if<std::is_same<AlsoT, RVec>::value>>
    rvec* rvec_array()
    {
        return as_rvec_array(data());
    }
    /*! \brief Returns a const rvec * pointer for containers of RVec, for use with legacy code.
     *
     * \todo Use std::is_same_v when CUDA 11 is a requirement.
     */
    template<typename AlsoT = T, typename = typename std::enable_if<std::is_same<AlsoT, RVec>::value>>
    const rvec* rvec_array() const
    {
        return as_rvec_array(data());
    }
    //! Copy assignment operator
    PaddedVector& operator=(PaddedVector const& o)
    {
        if (&o != this)
        {
            storage_     = o.storage_;
            unpaddedEnd_ = begin() + o.size();
        }
        return *this;
    }
    //! Move assignment operator
    PaddedVector& operator=(PaddedVector&& o) noexcept
    {
        if (&o != this)
        {
            auto oSize     = o.size();
            storage_       = std::move(o.storage_);
            unpaddedEnd_   = begin() + oSize;
            o.unpaddedEnd_ = o.begin();
        }
        return *this;
    }
    //! Getter for the allocator
    allocator_type get_allocator() const { return storage_.get_allocator(); }

private:
    storage_type storage_;
    iterator     unpaddedEnd_;
};

} // namespace gmx

// TODO These are hacks to avoid littering gmx:: all over code that is
// almost all destined to move into the gmx namespace at some point.
// An alternative would be about 20 files with using statements.
using gmx::PaddedVector; //NOLINT(google-global-names-in-headers)

#endif
