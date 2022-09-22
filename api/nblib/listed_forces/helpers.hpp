/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
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
/*! \internal \file
 * \brief
 * Helper data structures and utility functions for the nblib force calculator.
 * Intended for internal use.
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 */

#ifndef NBLIB_LISTEDFORCSES_HELPERS_HPP
#define NBLIB_LISTEDFORCSES_HELPERS_HPP

#include <unordered_map>

#include "gromacs/utility/arrayref.h"

#include "nblib/listed_forces/definitions.h"
#include "nblib/util/util.hpp"

#include "pbc.hpp"

#define NBLIB_ALWAYS_INLINE __attribute((always_inline))

namespace nblib
{

namespace detail
{
template<class T>
inline void gmxRVecZeroWorkaround([[maybe_unused]] T& value)
{
}

template<>
inline void gmxRVecZeroWorkaround<gmx::RVec>(gmx::RVec& value)
{
    for (int i = 0; i < dimSize; ++i)
    {
        value[i] = 0;
    }
}
} // namespace detail

/*! \internal \brief proxy object to access forces in an underlying buffer
 *
 * Depending on the index, either the underlying main buffer, or local
 * storage for outliers is accessed. This object does not own the main buffer.
 *
 */
template<class T>
class ForceBufferProxy
{
    using HashMap = std::unordered_map<int, T>;

public:
    ForceBufferProxy() : rangeStart_(0), rangeEnd_(0) { }

    ForceBufferProxy(int rangeStart, int rangeEnd) : rangeStart_(rangeStart), rangeEnd_(rangeEnd)
    {
    }

    void clearOutliers() { outliers.clear(); }

    inline NBLIB_ALWAYS_INLINE T& operator[](int i)
    {
        if (i >= rangeStart_ && i < rangeEnd_)
        {
            return mainForceBuffer[i];
        }
        else
        {
            if (outliers.count(i) == 0)
            {
                T zero = T();
                // if T = gmx::RVec, need to explicitly initialize it to zeros
                detail::gmxRVecZeroWorkaround(zero);
                outliers[i] = zero;
            }
            return outliers[i];
        }
    }

    typename HashMap::const_iterator begin() { return outliers.begin(); }
    typename HashMap::const_iterator end() { return outliers.end(); }

    [[nodiscard]] bool inRange(int index) const { return (index >= rangeStart_ && index < rangeEnd_); }

    void setMainBuffer(gmx::ArrayRef<T> buffer) { mainForceBuffer = buffer; }

private:
    gmx::ArrayRef<T> mainForceBuffer;
    int rangeStart_;
    int rangeEnd_;

    HashMap outliers;
};

namespace detail
{

static int computeChunkIndex(int index, int totalRange, int nSplits)
{
    if (totalRange < nSplits)
    {
        // if there's more threads than particles
        return index;
    }

    int splitLength = totalRange / nSplits;
    return std::min(index / splitLength, nSplits - 1);
}

} // namespace detail


/*! \internal \brief splits an interaction tuple into nSplits interaction tuples
 *
 * \param interactions
 * \param totalRange the number of particle sequence coordinates
 * \param nSplits number to divide the total work by
 * \return
 */
inline
std::vector<ListedInteractionData> splitListedWork(const ListedInteractionData& interactions,
                                                   int                          totalRange,
                                                   int                          nSplits)
{
    std::vector<ListedInteractionData> workDivision(nSplits);

    auto splitOneElement = [totalRange, nSplits, &workDivision](const auto& inputElement) {
        // the index of inputElement in the ListedInteractionsTuple
        constexpr int elementIndex =
                FindIndex<std::decay_t<decltype(inputElement)>, ListedInteractionData>{};

        // for now, copy all parameters to each split
        // Todo: extract only the parameters needed for this split
        for (auto& workDivisionSplit : workDivision)
        {
            std::get<elementIndex>(workDivisionSplit).parameters = inputElement.parameters;
        }

        // loop over all interactions in inputElement
        for (const auto& interactionIndex : inputElement.indices)
        {
            // each interaction has multiple coordinate indices
            // we must pick one of them to assign this interaction to one of the output index ranges
            // Todo: count indices outside the current split range in order to minimize the buffer size
            int representativeIndex =
                    *std::min_element(begin(interactionIndex), end(interactionIndex) - 1);
            int splitIndex = detail::computeChunkIndex(representativeIndex, totalRange, nSplits);

            std::get<elementIndex>(workDivision[splitIndex]).indices.push_back(interactionIndex);
        }
    };

    // split each interaction type in the input interaction tuple
    for_each_tuple(splitOneElement, interactions);

    return workDivision;
}

} // namespace nblib

#undef NBLIB_ALWAYS_INLINE

#endif // NBLIB_LISTEDFORCSES_HELPERS_HPP
