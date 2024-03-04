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
 * Declares gmx::ListOfLists
 *
 * \author Berk Hess <hess@kth.se>
 * \inpublicapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_LISTOFLISTS_H
#define GMX_UTILITY_LISTOFLISTS_H

#include <vector>

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/exceptions.h"

namespace gmx
{

/*! \brief A list of lists, optimized for performance
 *
 * This class holds a list of \p size() lists of elements of type \p T.
 * To optimize performance, the only modification operation supporting
 * is adding a new list at the end of the list of lists.
 *
 * This implementation stores all data internally in two std::vector objects
 * and thereby avoids the overhead of managing \p size() separate objects
 * in memory.
 *
 * Internal storage consists of one std::vector<int> listRanges_ of size number
 * of lists plus one and a std::vector<T> elements_ with the elements of all
 * lists concatenated. List i is stored in entries listRanges_[i] to
 * listRanges_[i+1] in elements_.
 *
 * \tparam T value type
 */

template<typename T>
class ListOfLists
{
public:
    //! Constructs an empty list of lists
    ListOfLists() = default;

    /*! \brief Constructs a list of list from raw data in internal layout
     *
     * Does basic consistency checks and throws when one of those fail.
     *
     * \param[in] listRanges  Ranges of the lists concatenated (see above), is consumed
     * \param[in] elements    Elements for all lists concatenated, is consumed
     */
    ListOfLists(std::vector<int>&& listRanges, std::vector<T>&& elements) :
        listRanges_(std::move(listRanges)), elements_(std::move(elements))
    {
        if (listRanges_.empty() || listRanges_.at(0) != 0)
        {
            GMX_THROW(InconsistentInputError(
                    "listRanges does not have a first element with value 0"));
        }
        if (int(elements_.size()) != listRanges_.back())
        {
            GMX_THROW(InconsistentInputError(
                    "The size of elements does not match the last value in listRanges"));
        }
    }

    //! Returns the number of lists
    std::size_t size() const { return listRanges_.size() - 1; }

    /*! \brief Returns the number of lists
     *
     * \note Use ssize for any expression involving arithmetic operations
     * (including loop indices).
     */
    Index ssize() const { return Index(listRanges_.size()) - 1; }

    //! Returns whether the list holds no lists
    bool empty() const { return listRanges_.size() == 1; }

    //! Returns the sum of the number of elements over all lists
    int numElements() const { return listRanges_.back(); }

    //! Appends a new list with elements \p values, pass {} to add an empty list
    void pushBack(ArrayRef<const T> values)
    {
        elements_.insert(elements_.end(), values.begin(), values.end());
        listRanges_.push_back(int(elements_.size()));
    }

    /*! \brief Appends a new list with \p numElements elements
     *
     * \note Type T should be default constructible.
     */
    void pushBackListOfSize(int numElements)
    {
        static_assert(std::is_default_constructible_v<T>);

        elements_.resize(elements_.size() + numElements);
        listRanges_.push_back(int(elements_.size()));
    }

    //! Returns an ArrayRef to the elements of the list with the given index
    ArrayRef<const T> operator[](std::size_t listIndex) const
    {
        return ArrayRef<const T>(elements_.data() + listRanges_[listIndex],
                                 elements_.data() + listRanges_[listIndex + 1]);
    }

    //! Returns the list of elements for the list with index \p listIndex, throws an \p out_of_range exception when out of range
    ArrayRef<const T> at(std::size_t listIndex) const
    {
        return ArrayRef<const T>(elements_.data() + listRanges_.at(listIndex),
                                 elements_.data() + listRanges_.at(listIndex + 1));
    }

    /*! \brief Returns a reference to the first list
     *
     * \returns a reference to the first list
     */
    ArrayRef<T> front()
    {
        GMX_ASSERT(size() > 0, "Must contain a list if front() is called");
        auto* beginPtr = elements_.data();
        auto* endPtr   = beginPtr + listRanges_[1];
        return { beginPtr, endPtr };
    }
    /*! \brief Returns a reference to the final list
     *
     * \returns a reference to the final list
     */
    ArrayRef<T> back()
    {
        GMX_ASSERT(size() > 0, "Must contain a list if bank() is called");
        auto endIndex   = *(listRanges_.end() - 1);
        auto beginIndex = *(listRanges_.end() - 2);
        return { elements_.data() + beginIndex, elements_.data() + endIndex };
    }

    //! Clears the list
    void clear()
    {
        listRanges_.resize(1);
        elements_.clear();
    }

    //! Appends a ListOfLists at the end and increments the appended elements by \p offset
    void appendListOfLists(const ListOfLists& listOfLists, const T offset = 0)
    {
        listRanges_.insert(
                listRanges_.end(), listOfLists.listRanges_.begin() + 1, listOfLists.listRanges_.end());
        const int oldNumElements = elements_.size();
        for (std::size_t i = listRanges_.size() - listOfLists.size(); i < listRanges_.size(); i++)
        {
            listRanges_[i] += oldNumElements;
        }
        elements_.insert(elements_.end(), listOfLists.elements_.begin(), listOfLists.elements_.end());

        if (offset != 0)
        {
            for (std::size_t i = elements_.size() - listOfLists.elements_.size(); i < elements_.size(); i++)
            {
                elements_[i] += offset;
            }
        }
    }

    //! Returns concatenated ranges of the lists (see above for details)
    ArrayRef<const int> listRangesView() const { return listRanges_; }

    //! Returns the a view of the elements of all lists concatenated
    ArrayRef<const T> elementsView() const { return elements_; }

private:
    //! The ranges of the lists, list i uses range \p listRanges_[i], \p listRanges_[i+1].
    std::vector<int> listRanges_ = { 0 };
    //! The elements in all lists concatenated
    std::vector<T> elements_;
};

} // namespace gmx

#endif
