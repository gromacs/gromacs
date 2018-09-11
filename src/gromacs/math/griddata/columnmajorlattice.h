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
/*! \file
 * \brief
 * Declares gmx::ColumnMajorLattice template class.
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \inpublicapi
 * \ingroup module_griddata
 */
#ifndef GMX_MATH_COLUMNMAJORLATTICE
#define GMX_MATH_COLUMNMAJORLATTICE

#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/iserializer.h"

#include <algorithm>
#include <array>
#include <numeric>
#include <functional>

namespace gmx
{

//! N-dimensional integer index.
template<size_t Rank>
class offset
{
    public:
        static constexpr size_t rank = Rank;
        using reference              = std::ptrdiff_t&;
        using const_reference        = const std::ptrdiff_t&;
        using size_type              = size_t;
        using value_type             = std::ptrdiff_t;
        constexpr offset() noexcept
        {}

        offset(const std::array<value_type, Rank> &offset)
        {
            offset_ = offset;
        }

        offset(std::initializer_list<value_type> il)
        {
            std::copy(il.begin(), il.end(), offset_.data());
        };

        reference operator[](size_type n) { return offset_[n]; }
        constexpr const_reference operator[](size_type n) const { return offset_[n]; }
    private:
        std::array<value_type, Rank> offset_ = {};
};

template <size_t Rank> class bounds_iterator;

template<size_t Rank>
class bounds
{
    public:
        static constexpr size_t rank = Rank;
        using reference              = std::ptrdiff_t&;
        using const_reference        = const std::ptrdiff_t&;
        using iterator               = bounds_iterator<Rank>;
        using const_iterator         = bounds_iterator<Rank>;
        using size_type              = size_t;
        using value_type             = std::ptrdiff_t;

        constexpr bounds() noexcept {};

        bounds(std::initializer_list<value_type> il)
        {
            std::copy(il.begin(), il.end(), bounds_.data());
        };

        constexpr size_type size() const noexcept
        {
            return std::accumulate(std::begin(bounds_), std::end(bounds_), 1, std::multiplies<int>());
        };

        bool contains(const offset<Rank> &latticeIndex) const noexcept
        {
            for (size_t dim = 0; dim < Rank-1; ++dim)
            {
                if ( (latticeIndex[dim] < 0) || (bounds_[dim] <= latticeIndex[dim]))
                {
                    return false;
                }
            }
            return true;
        };

        const_iterator begin() const noexcept
        {
            return *this;
        };

        const_iterator end() const noexcept
        {
            return {*this, bounds_}; // if iterator offset is at bounds, the end is reached
        }

        reference       operator[](size_type n) { return bounds_[n]; }
        constexpr const_reference operator[](size_type n) const { return bounds_[n]; }

    private:
        std::array<value_type, Rank> bounds_ = {};
};

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

template <size_t Rank>
constexpr bool operator!=(const offset<Rank> &lhs, const offset<Rank> &rhs) noexcept
{
    return !(lhs == rhs);
};

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

template <size_t Rank>
constexpr bool operator!=(const bounds<Rank> &lhs, const bounds<Rank> &rhs) noexcept
{
    return !(lhs == rhs);
};

/*!\brief
 * Note: The iteration order deviates from the cpp core proposal:
 * first dimension moves fastest.
 */
template <size_t Rank>
class bounds_iterator
{
    public:
        using iterator_category = std::forward_iterator_tag;
        using value_type        = offset<Rank>;
        using difference_type   = std::ptrdiff_t;
        using pointer           = offset<Rank>*;
        using reference         = const offset<Rank>;

        bounds_iterator(const bounds<Rank> &bounds, offset<Rank> offset) noexcept
            : bounds_(bounds), offset_(offset)
        {
        }

        bounds_iterator(const bounds<Rank> &bounds) noexcept
            : bounds_(bounds), offset_() {}

        bool operator==(const bounds_iterator &rhs) const
        {
            return offset_ == rhs.offset_;
        }

        bool operator!=(const bounds_iterator &rhs) const
        {
            return !(*this == rhs);
        }

        bounds_iterator &operator++()
        {
            for (size_t dim = 0; dim < Rank; ++dim)
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

            // this is the end
            for (size_t dim = 0; dim < Rank; ++dim)
            {
                offset_[dim] = bounds_[dim];
            }

            return *this;

        };

        bounds_iterator  operator++(int)
        {
            bounds_iterator tmp(*this);
            ++(*this);
            return tmp;
        }

        reference operator*() const { return offset_; }
        pointer   operator->() const { return &offset_; }

    private:
        bounds<Rank> bounds_;
        offset<Rank> offset_;
        ptrdiff_t    position_;
};


/*! \brief
 * An N-dimensional lattice with column major indexing.
 *
 * The first dimension varies quickest and is most contigiuous.
 * Vector-like multi indices are linearised via
 * i_linear= i_0 + extend_0 * (i_1 + extend_1 * (i_2 + extend_2 * ( ... ))).
 * Lattice extend may not be altered, construct a new lattice instead.
 *
 * \tparam int Number of lattice dimensions.
 * \throws gmx::RangeError if indices are out of lattice bounds.
 */
// template <size_t Rank> class ColumnMajorLattice
// {
//     public:
//
//         /*! \brief
//          * Constructs Lattice by setting its extend and checking extend validity.
//          *
//          * Indices span (0,...,0) to (extend[0]-1,...,extend[N]-1)
//          * \TODO check for overflow of number of lattice points
//          * \throws gmx::RangeError if any extend is not larger than zero.
//          */
//         /*! \brief
//          * The lattice extend as N-dimensional index.
//          */
//         /*! \brief
//          * Compare lattices.
//          */
//         /*! \brief
//          * Lattice are unequal if not equal.
//          */
//         //  auto  stride = bounds_.size();
//         //  // Build a multi index from a linear index,
//         //  // starting from the slowest moving dimension (from the back)
//         //  // Integer division by the stride yields the vector index,
//         //  // the remainder is used for calculating the remaining vector indices
//         //  auto currentLatticeIndex = offset_.rbegin();
//         //  for (auto currentDimensionExtend = bounds_.rbegin();
//         //       currentDimensionExtend != extend_.rend(); ++currentDimensionExtend)
//         //  {
//         //      stride              /= *currentDimensionExtend;
//         //      *currentLatticeIndex = linearIndex / stride;
//         //      linearIndex         -= *currentLatticeIndex * stride;
//         //      ++currentLatticeIndex;
//         //  }
//         //
//         /*! \brief
//          * Returns one-dimensional lattice index from N-dimensional multi index.
//          *
//          * i_0 + extend_0 * (i_1 + extend_1 * (i_2 + extend_2 * ( ... ))).
//          * \throws gmx::RangeError if latticeIndex is not within the lattice.
//          * \returns non-negative integer lattice index
//          */
//         int lineariseVectorIndex(const MultiIndex &latticeIndex) const
//         {
//             auto linearIndex            = *latticeIndex.rbegin();
//             auto currentDimensionExtend = extend_.rbegin();
//             for (auto currentDimensionIndex = ++latticeIndex.rbegin();
//                  currentDimensionIndex != latticeIndex.rend();
//                  ++currentDimensionIndex)
//             {
//                 ++currentDimensionExtend;
//                 linearIndex = *currentDimensionIndex + *currentDimensionExtend * linearIndex;
//             }
//
//             return linearIndex;
//         }
//
//         /*! \brief
//          * Generate the multi index from a linear lattice index.
//          *
//          * The inverse of lineariseVectorIndex.
//          *
//          * \param[in] linearIndex non-negative and smaller than number of lattice points.
//          * \returns the vectorised linear index
//          * \throws gmx::RangeError if linearIndex exceeds the number of lattice points.
//          */
//         MultiIndex vectoriseLinearIndex(int linearIndex) const
//         {
//             MultiIndex result;
//             auto       stride = getNumLatticePoints();
//             // Build a multi index from a linear index,
//             // starting from the slowest moving dimension (from the back)
//             // Integer division by the stride yields the vector index,
//             // the remainder is used for calculating the remaining vector indices
//             auto currentLatticeIndex = result.rbegin();
//             for (auto currentDimensionExtend = extend_.rbegin();
//                  currentDimensionExtend != extend_.rend(); ++currentDimensionExtend)
//             {
//                 stride              /= *currentDimensionExtend;
//                 *currentLatticeIndex = linearIndex / stride;
//                 linearIndex         -= *currentLatticeIndex * stride;
//                 ++currentLatticeIndex;
//             }
//             return result;
//         }
//
//         /*! \brief
//          * Check if latticeIndex is a valid index for this lattice.
//          *
//          * \param[in] latticeIndex Lattice index to check
//          * \returns True if latticeIndex is in lattice
//          */
//         bool inLattice(const MultiIndex &latticeIndex) const
//         {
//         }
//
//         /*\brief Serialise from master to all other nodes.
//          * TODO: this should live somewhere else, since non-parallel uses of
//          * the lattice indices will never need to see the broadcast functionality
//          */
//         void readwrite(ISerializer * serializer)
//         {
//             for (auto &e : extend_)
//             {
//                 serializer->doInt(&e);
//             }
//         }
//
//     private:
//         //! \brief the extend of the lattice
//         MultiIndex extend_;
// };

} // gmx

#endif
