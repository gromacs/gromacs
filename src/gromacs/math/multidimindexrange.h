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
/*! \libinternal
 * \file
 * \brief GROMACS extensions to mdspan.
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \inlibraryapi
 * \ingroup module_math
 */
#include <algorithm>
#include <array>
#include <cstddef>
#include <vector>

namespace gmx
{

/*! \libinternal
 * \brief A range within the integer hyperrectangle [begin, end).
 *
 * Enables range based looping over lattice positions, e.g.,
 *  for( auto position : MultiDimIndexRange<2>({10,10}))
 *
 * End is defined as one element past last in all dimensions -  with begin, b and end e
 * a range encompases the following lattice points (.), not including (- and |)
 *
 *                   ------------------e
 *                   . . . . . . . . . |
 *                   b . . . . . . . . |
 *
 * \tparam Rank the dimensionality of the lattice to be iterated over
 */
template <int Rank>
class MultiDimIndexRange
{
    public:
        /*! \libinternal
         * \brief Iterator visiting multidimensional lattice points within the
         *        hyperrectangle [begin, end).
         *
         * \note For performance resons, prefer direct iterator over multidimensional
         *       arrays if lattice coordinates are not required.
         *
         * The iterator visits lattice points, changing the highest, right-most index
         * first - in line with the memory access pattern of layout_right in mdspan.
         */
        class iterator
        {
            public:
                //! This iterator only implements the requirements for an input interator
                using iterator_category = std::input_iterator_tag;
                //! Value type
                using value_type = std::array<std::ptrdiff_t, Rank>;
                //! Difference type
                using difference_type = std::ptrdiff_t;
                //! Pointer type
                using pointer = std::array<std::ptrdiff_t, Rank> *;
                //! Reference type
                using reference = std::array<std::ptrdiff_t, Rank> &;
                /*! \brief Construct iterator from begin and end and set an intial position.
                 * \param[in] begin The begin of the range to be iterated over.
                 * \param[in] end   The end of the range to be iterated over
                 * \param[in] position Set a starting position of the iterator, null per default.
                 */
                constexpr iterator(const value_type &begin,
                                   const value_type &end,
                                   const value_type &position = {0}
                                   ) :
                    begin_ {begin}, end_ {end}, position_ {position}
                { }
                /*! \brief Prefix increment moves the position on the lattice with
                 *        the last rank moving the fastest.
                 */
                iterator &operator++()
                {
                    // The last rank is the fastest moving, in line with layout_right
                    for (int dim = Rank - 1; dim >= 0; --dim)
                    {
                        ++position_[dim];
                        if (position_[dim] < end_[dim])
                        {
                            return *this;
                        }
                        position_[dim] = begin_[dim];
                    }
                    // iterator hit end_ in all dimensions and is again (begin[0],...,begin[n]).
                    // Set to end instead.
                    position_ = end_;
                    return *this;
                }
                //! Equal if iterating over the same range and being at the same position
                constexpr bool operator==(const iterator &rhs) const noexcept
                {
                    return (begin_ == rhs.begin_) &&
                           (end_ == rhs.end_) &&
                           (position_ == rhs.position_);
                }
                //! not equal
                constexpr bool operator!=(const iterator &rhs) const noexcept
                {
                    return !(*this == rhs);
                }
                //! Dereferencing returns the current position on the lattice as an array.
                constexpr const value_type &operator*() const noexcept
                {
                    return position_;
                }

            private:
                //! The begin of the region to iterate over
                value_type begin_;
                //! The end of the region to be iterated over.
                value_type end_;
                //! The current position within the lattice.
                value_type position_;

        };

        /*! \brief Construct index range by setting begin and end position.
         *
         * If any elemnt of the beginPosition exceeds the endPosition, the
         * begin is set to the end, generating an empty range.
         *
         * \note Using std::less within std::equal allows element-wise checking
         *       wether beginPosition exceeds endPosition while keeping this
         *       constructor a constexpr, though substituting equal with less leads
         *       to counterintuitive syntax here.
         */
        constexpr MultiDimIndexRange(const typename iterator::value_type &beginPosition,
                                     const typename iterator::value_type &endPosition
                                     ) :
            beginIterator_
        {
            beginPosition, endPosition,
            std::equal(std::begin(beginPosition), std::end(beginPosition),
                       std::begin(endPosition), std::less<std::ptrdiff_t>()) ? beginPosition : endPosition
        },
        endIterator_ {beginPosition, endPosition, endPosition}
        {}

        /*! \brief Construct a MultiDimIndexRange from a lattice extent, setting the begin to zero.
         * \param[in] end   of the range to be iterated over
         */
        explicit constexpr MultiDimIndexRange(const typename iterator::value_type &end) :
            MultiDimIndexRange({0}, end)
        {}

        //! return iterator at begin of index range
        constexpr const iterator &begin() const noexcept { return beginIterator_; }
        //! return iterator at end of index range
        constexpr const iterator &end() const noexcept { return endIterator_; }

    private:
        //! First position on the grid to be iterated over
        const iterator beginIterator_;
        /*! \brief Past-last position to be iterated over.
         * One larger than the last position in all dimensions.
         */
        const iterator endIterator_;
};

//! Compare two index ranges for equality.
template<int Rank>
bool operator==(const MultiDimIndexRange<Rank> &firstRange,
                const MultiDimIndexRange<Rank> &otherRange)
{
    return (firstRange.begin() == otherRange.begin()) &&
           (firstRange.end()   == otherRange.end()  );
};

//! Compare two index ranges for inequality.
template<int Rank>
bool operator!=(const MultiDimIndexRange<Rank> &firstRange,
                const MultiDimIndexRange<Rank> &otherRange)
{
    return !(firstRange == otherRange);
};


/*! \libinternal
 * \brief calculate the intersection between two multidimensional index ranges.
 *
 * The element-wise maximum of input ranges determines the intersection range begin,
 * the element-wise maximum of input ranges its end.
 *
 * \tparam Rank the dimensionality of the index range
 * \returns intersection of two index ranges
 */
template <int Rank>
MultiDimIndexRange<Rank> multiDimIndexRangeIntersection(const MultiDimIndexRange<Rank> &first,
                                                        const MultiDimIndexRange<Rank> &second)
{

    typename MultiDimIndexRange<Rank>::iterator::value_type intersectionBegin;
    std::transform(std::begin(*first.begin()), std::end(*first.begin()),
                   std::begin(*second.begin()), std::begin(intersectionBegin),
                   [](auto a, auto b){return std::max(a, b); });

    //
    typename MultiDimIndexRange<Rank>::iterator::value_type intersectionEnd;
    std::transform(std::begin(*first.end()), std::end(*first.end()),
                   std::begin(*second.end()), std::begin(intersectionEnd),
                   [](auto a, auto b){return std::min(a, b); });

    return {intersectionBegin, intersectionEnd};
}
} // namespace gmx
