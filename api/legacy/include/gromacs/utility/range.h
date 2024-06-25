/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
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
 * Declares gmx::Range
 *
 * \author Berk Hess <hess@kth.se>
 * \inpublicapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_RANGE_H
#define GMX_UTILITY_RANGE_H

#include <iterator>
#include <type_traits>

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/gmxassert.h"

namespace gmx
{

/*! \brief Defines a range of integer numbers and accompanying operations
 *
 * Defines a range of integer type with begin and end.
 * Can be used in a range loop over all integers in the range as in:
 *
 * Range<int> range(2,5);
 * for (int i : range) { printf(" %d", i); }
 *
 * This will print: 2 3 4
 */
template<typename T>
class Range
{
    // TODO: Use std::is_integral_v when CUDA 11 is a requirement.
    static_assert(std::is_integral<T>::value, "Range can only be used with integral types");

    // Note: This class has as invariant: begin_ <= end_

public:
    //! An iterator that loops over a range of integers
    struct iterator
    {
        using iterator_category = std::forward_iterator_tag;
        using difference_type   = T;
        using value_type        = T;
        using pointer           = T*;
        using reference         = T&;

        //! Constructor
        iterator(T value) : value_(value) {}
        //! Value
        operator T() const { return value_; }
        //! Reference
        operator T&() { return value_; }
        //! Pointer
        T operator*() const { return value_; }
        //! Inequality comparison
        bool operator!=(const iterator other) { return value_ != other.value_; }
        //! Increment operator
        iterator& operator++()
        {
            ++value_;
            return *this;
        }
        //! Increment operator
        iterator operator++(int gmx_unused dummy)
        {
            iterator tmp(*this);
            ++value_;
            return tmp;
        }
        //! The actual value
        T value_;
    };

    //! Constructor, has to be called with \p begin <= \p end (is checked)
    Range(const T begin, const T end) : begin_(begin), end_(end)
    {
        GMX_RELEASE_ASSERT(begin_ <= end_, "A range should have begin<=end");
    }

    //! Default constructor, produces an empty range
    Range() = default;

    //! Begin iterator/value
    iterator begin() const { return begin_; }
    //! End iterator/value
    iterator end() const { return end_; }

    //! Returns the length of the range
    T size() const { return end_ - begin_; }

    //! Returns whether the range is empty
    bool empty() const { return size() == 0; }

    //! Returns whether \p value is in range
    bool isInRange(const T value) const { return (begin_ <= value && value < end_); }

private:
    //! The start of the range
    T begin_ = 0;
    //! The end of the range
    T end_ = 0;
};

} // namespace gmx

#endif
