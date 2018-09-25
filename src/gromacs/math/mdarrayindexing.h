/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
/*! \internal
 * \file
 * \brief
 * Declares and imlements methods for indexing of multidinemsional data.
 * Implements a subset of C++ propoal N3851 "Multidimensional bounds, index and array_view"
 * and mdspan.
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \ingroup module_griddata
 */
#ifndef GMX_MATH_MDARRAYINDEXING
#define GMX_MATH_MDARRAYINDEXING

#include <array>
#include <numeric>

namespace gmx
{

/*! \internal
 * \brief Multidimensional offset from zero to index an array.
 * Indexing is zero based. Refer to N3851.
 * \tparam Rank the rank of the data
 */
template<size_t Rank>
class offset
{
    public:
        //! \brief Expose the rank template parameter
        static constexpr size_t rank = Rank;
        //! \brief Reference data type.
        using reference              = std::ptrdiff_t&;
        //! \brief Const reference data type.
        using const_reference        = const std::ptrdiff_t&;
        //! \brief Size type for the size of bounds.
        using size_type              = size_t;
        //! \brief Value type.
        using value_type             = std::ptrdiff_t;

        /*! \brief Default constructor.
         */
        constexpr offset() noexcept = default;

        /*! \brief Construct offset from plain C array of size Rank.
         * \param[in] offset The offset as array.
         * NOTE Not part of N3851.
         * This allows conversion from gmx::IVec for three dimensional grids.
         */
        offset(const int (&offset)[Rank])
        {
            for (size_t dim = 0; dim < Rank; ++dim)
            {
                offset_[dim] = offset[dim];
            }
        }

        /*! \brief Construct offset from array.
         * \param[in] offset The offset as array.
         * NOTE Not part of N3851.
         */
        constexpr offset(const std::array<value_type, Rank> &offset) : offset_ {offset}
        {}

        /*! \brief Construct offset from initializer list.
         * \param[in] il Initializer list
         *
         * Avoids double brace construction syntax for array.
         */
        offset(std::initializer_list<value_type> il)
        {
            std::copy(il.begin(), il.end(), offset_.data());
        }

        /*! \brief Access the n-th offset.
         * \param[in] n The dimension of the offset.
         */
        reference operator[](size_type n) { return offset_[n]; }
        /*! \brief Access the constant reference to the n-th offset.
         * \param[in] n The dimension of the offset.
         */
        constexpr const_reference operator[](size_type n) const { return offset_[n]; }
    private:
        //!\brief The offset.
        std::array<value_type, Rank> offset_ = {};
};


template <size_t Rank> class bounds_iterator;

/*! \internal
 * \brief The bounds of a multidimensional array.
 * The total size of the array is the product of bounds in all dimensions.
 * Multidimensional iteration within bounds is achieved with bounds_iterator
 * that are provided by the begin() and end() method.
 * \tparam Rank The dimension of the array
 */
template<size_t Rank>
class bounds
{
    public:
        //! \brief Expose the rank template parameter
        static constexpr size_t rank = Rank;
        //! \brief Reference data type.
        using reference              = std::ptrdiff_t&;
        //! \brief Const reference data type.
        using const_reference        = const std::ptrdiff_t&;
        //! \brief Iterator type.
        using iterator               = bounds_iterator<Rank>;
        //! \brief Const iterator type.
        using const_iterator         = bounds_iterator<Rank>;
        //! \brief Size type for the size of bounds.
        using size_type              = size_t;
        //! \brief Value type.
        using value_type             = std::ptrdiff_t;


        /*!\brief The default constructor.
         */
        constexpr bounds() noexcept = default;

        /*!\brief Construct bounds from intializer list.
         * Avoids double brace construction syntax for array.
         */
        bounds(std::initializer_list<value_type> il)
        {
            std::copy(il.begin(), il.end(), bounds_.data());
        }

        /*! \brief Construct bounds from plain C array.
         * \param[in] bounds The bounds as array.
         * NOTE Not part of N3851.
         * This allows conversion from gmx::IVec for three dimensional grids.
         */
        bounds(const int (&bounds)[Rank])
        {
            for (size_t dim = 0; dim < Rank; ++dim)
            {
                bounds_[dim] = bounds[dim];
            }
        }

        /*! \brief Construct bounds from an array.
         * \param[in] bounds The bounds as array.
         * NOTE Not part of N3851.
         */
        bounds(const std::array<value_type, Rank> &bounds)
        {
            std::copy(std::begin(bounds), std::end(bounds), std::begin(bounds_));
        }

        /*!\brief The total size of the array.
         */
        constexpr size_type size() const noexcept
        {
            return std::accumulate(std::begin(bounds_), std::end(bounds_), 1, std::multiplies<int>());
        }

        /*!\brief Tests if an offset is contained within the bounds of an array.
         * \param[in] offset The offset to be checked if within bounds.
         * \returns offset larger zero and smaller bounds for all dimension smaller than Rank.
         */
        bool contains(const offset<Rank> &offset) const noexcept
        {
            for (size_t dim = 0; dim < Rank; ++dim)
            {
                if ( (offset[dim] < 0 ) || (offset[dim] >= bounds_[dim]))
                {
                    return false;
                }
            }
            return true;
        }

        /*!\brief Iterator to begin of the multidimensional array.
         * \returns iterator at begin
         */
        const_iterator begin() const noexcept
        {
            return *this;
        }
        /*!\brief Iterator to end of the multidimensional array.
         * \returns iterator at end
         */
        const_iterator end() const noexcept
        {
            return {*this, bounds_}; // define end as offset at bounds
        }
        /*!\brief Accesses array extent at dimension.
         * \param[in] n Dimension to access.
         * \returns Array extent
         */
        reference operator[](size_type n) { return bounds_[n]; }
        /*!\brief Accesses bound at dimension.
         * \param[in] n Dimension to access.
         * \returns Array extent
         */
        constexpr const_reference operator[](size_type n) const { return bounds_[n]; }

    private:
        //! The bounds of an N-dimensional array.
        std::array<value_type, Rank> bounds_ = {};
};

/*! \brief Check if two offsets are equal.
 * \tparam Rank the dimension of the offsets.
 * \param[in] lhs first offset to compare
 * \param[in] rhs second offset to compare
 * \returns true if two offsets are equal.
 */
template <size_t Rank>
bool operator==(const offset<Rank> &lhs, const offset<Rank> &rhs) noexcept
{
    for (size_t i = 0; i < Rank; i++)
    {
        if (lhs[i] != rhs[i])
        {
            return false;
        }
    }
    return true;
};

/*! \brief Check if two offsets are not equal.
 * \tparam Rank the dimension of the offsets.
 * \param[in] lhs first offset to compare
 * \param[in] rhs second offset to compare
 * \returns true if two offsets are not equal.
 */
template <size_t Rank>
constexpr bool operator!=(const offset<Rank> &lhs, const offset<Rank> &rhs) noexcept
{
    return !(lhs == rhs);
};

/*! \brief Check if two bounds are equal.
 * \tparam Rank the dimension of the bounds.
 * \param[in] lhs first bounds to compare
 * \param[in] rhs second bounds to compare
 * \returns true if two bounds are equal.
 */
template <size_t Rank>
bool operator==(const bounds<Rank> &lhs, const bounds<Rank> &rhs) noexcept
{
    for (size_t i = 0; i < Rank; i++)
    {
        if (lhs[i] != rhs[i])
        {
            return false;
        }
    }
    return true;
};

/*! \brief Check if two bounds are different.
 * \tparam Rank the dimension of the bounds.
 * \param[in] lhs first bounds to compare
 * \param[in] rhs second bounds to compare
 * \returns true if two bounds are not equal.
 */
template <size_t Rank>
constexpr bool operator!=(const bounds<Rank> &lhs, const bounds<Rank> &rhs) noexcept
{
    return !(lhs == rhs);
};

/*! \internal
 * \brief Proxy iterator over multidimensional array elements within bounds.
 *
 * Steps through elements in a multidimensional array.
 * The highest rank index is incremented first. If any index hits the array
 * bounds, it is reset to zero and the next lower rank index is incremented.
 * E.g.: Within bounds (3,5,2), the iterator moves
 * (2,3,0) -> (2,3,1) -> (2,4,0) -> (2,4,1) -> end
 *
 * NOTE in conformance with N3851, this is a proxy iterator that does not return
 * access to an underlying object, but rather a multidimensional offset.
 */
template <size_t Rank>
class bounds_iterator
{
    public:
        /*! \brief Type of iterator declared here.
         * NOTE: in contrast to N3851, only forward iterator routines are
         * implemented, because this all needed in GROMACS right now.
         * Therefore there is no difference_type=std::ptrdiff_t declared here.
         */
        using iterator_category = std::forward_iterator_tag;
        //! \brief Dereferencing iterator to comply std::iterator_traits.
        using value_type        = offset<Rank>;
        //! \brief Pointer to iterator value type to comply std::iterator_traits.
        using pointer           = offset<Rank>*;
        /*! \brief Refernce to iterator value type to comply std::iterator_traits.
         * NOTE: Dereferencing this iterator returns an offset rather than access to
         * an underlying object in compliance with N3851.
         */
        using reference         = const offset<Rank>;

        /*!\brief Construct an iterator within array bounds at given offset.
         * \param[in] bounds The bounds of the array to iterate over.
         * \param[in] offset The offset from where to start iteration.
         */
        bounds_iterator(const bounds<Rank> &bounds, const offset<Rank> &offset) noexcept
            : bounds_(bounds), offset_(offset)
        {
        }
        /*!\brief Construct an iterator within array bounds starting at zero.
         * \param[in] bounds The bounds of the array to iterate over.
         */
        bounds_iterator(const bounds<Rank> &bounds) noexcept
            : bounds_(bounds), offset_() {}

        /*!\brief Test iterators for equality.
         * \param[in] rhs The other iterator.
         * \returns true if this and compared iterator point to the same array offset.
         */
        bool operator==(const bounds_iterator &rhs) const
        {
            return offset_ == rhs.offset_;
        }
        /*!\brief Test iterators for inequality.
         * \param[in] rhs The other iterator.
         * \returns true if this and compared iterator do not point to the same array offset.
         */
        bool operator!=(const bounds_iterator &rhs) const
        {
            return !(*this == rhs);
        }

        /*!\brief Go to next offset in multidimensional array.
         * \returns iterator to offset position in multidimensional array.
         */
        bounds_iterator &operator++()
        {
            for (int dim = Rank-1; dim >= 0; --dim)
            {
                ++offset_[dim];
                if (offset_[dim] < bounds_[dim])
                {
                    return *this;
                }
                else
                {
                    offset_[dim] = 0;
                }
            }

            // iterator hit bounds in all dimensions and will
            // dereference to (0,...,0) now. Setting to end instead.
            for (size_t dim = 0; dim < Rank; ++dim)
            {
                offset_[dim] = bounds_[dim];
            }

            return *this;

        }
        /*!\brief Post-fix increment to reach next position in array.
         */
        bounds_iterator operator++(int)
        {
            bounds_iterator tmp(*this);
            ++(*this);
            return tmp;
        }
        /*!\brief Dereference iterator to offset in grid.
         * \returns the multidimensional offset the iterator points to.
         */
        reference operator*() const { return offset_; }

    private:
        //! The bounds within which the iterator moves.
        bounds<Rank> bounds_;
        //! The current iterator position
        offset<Rank> offset_;
};


} // namespace gmx

#endif
