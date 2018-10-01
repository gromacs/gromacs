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
 * \brief
 * Implements gmx::UnionFinder and gmx::MappedUnionFinder.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_trajectoryanalysis
 */
#ifndef GMX_TRAJECTORYANALYSIS_UNIONFIND_H
#define GMX_TRAJECTORYANALYSIS_UNIONFIND_H

#include <algorithm>
#include <numeric>
#include <vector>

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/gmxassert.h"

namespace gmx
{

/*! \libinternal \brief
 * Union-find data structure for keeping track of disjoint sets.
 *
 * Union-find keeps track of a number of items, represented here by continuous
 * integer indices starting at zero, and supports the following operations:
 *  - Initialization puts each item into a set of its own.
 *  - Given two items, merge the sets that contain these items.
 *  - Given an item, find a representative item that is in the same set, such
 *    that queries for items in the same set yield the same value.
 * Merging and querying is supported in amortized constant time.
 *
 * Note that in order to achieve the amortized behavior, querying the structure
 * modifies the internal state, but does not alter the externally visible
 * behavior.
 *
 * \ingroup module_trajectoryanalysis
 */
class UnionFinder
{
    public:
        /*! \brief
         * Initializes `count` items, putting each in its own set.
         */
        void init(int count)
        {
            parent_.clear();
            rank_.clear();
            parent_.reserve(count);
            rank_.reserve(count);
            parent_.resize(count);
            rank_.resize(count, 0);
            std::iota(parent_.begin(), parent_.end(), 0);
        }
        /*! \brief
         * Merges sets that contain two given items.
         *
         * If the items are already in the same set, nothing happens.
         */
        void merge(int item1, int item2)
        {
            GMX_ASSERT(item1 >= 0 && item1 < count(), "Input index out of range");
            GMX_ASSERT(item2 >= 0 && item2 < count(), "Input index out of range");
            const int root1 = findRootAndCompressPath(item1);
            const int root2 = findRootAndCompressPath(item2);
            if (root1 != root2)
            {
                mergeRoots(root1, root2);
            }
        }
        /*! \brief
         * Returns a representative item from the set containing `item`.
         */
        int representativeItem(int item)
        {
            GMX_ASSERT(item >= 0 && item < count(), "Input index out of range");
            return findRootAndCompressPath(item);
        }
        /*! \brief
         * Returns the sizes of all sets (in arbitrary order).
         */
        std::vector<int> allSizes()
        {
            const int        count = parent_.size();
            std::vector<int> result(count, 0);
            for (int i = 0; i < count; ++i)
            {
                ++result[findRootAndCompressPath(i)];
            }
            result.erase(std::remove(result.begin(), result.end(), 0),
                         result.end());
            return result;
        }

    private:
        //! Number of items.
        int count() const { return parent_.size(); }
        int findRootAndCompressPath(int i)
        {
            while (parent_[i] != i)
            {
                const int prev = i;
                i              = parent_[i];
                parent_[prev]  = parent_[i];
            }
            return i;
        }
        void mergeRoots(int root1, int root2)
        {
            if (rank_[root1] > rank_[root2])
            {
                parent_[root2] = root1;
            }
            else if (rank_[root2] > rank_[root1])
            {
                parent_[root1] = root2;
            }
            else
            {
                parent_[root1] = root2;
                ++rank_[root1];
            }
        }

        /*! \brief
         * Parent item for each item in the tree representing the set.
         *
         * Root items are parents of themselves, and are the reprensentative
         * items of their sets.
         */
        std::vector<int> parent_;
        //! Worst-case height for each root (as if no compression was done).
        std::vector<int> rank_;
};

/*! \libinternal \brief
 * Extension of UnionFind that supports non-consecutive integer indices as
 * items.
 *
 * Sometimes, it is more convenient to operate on a set of integers that do not
 * start at zero and are not consecutive as UnionFind expects.  This class
 * implements a mapping on top of UnionFind such that this is possible.
 *
 * The current implementation assumes that the indices are bounded between zero
 * and some reasonably small integer, i.e., the memory usage depends on the
 * largest index number, not just on the number of items.
 *
 * \ingroup module_trajectoryanalysis
 */
class MappedUnionFinder
{
    public:
        /*! \brief
         * Initializes the finder with indices.
         *
         * The size of `indices` sets the number of input items, and each
         * unique value in `indices` maps to a single internal item.
         * If multiple indices are the same, then these items are considered
         * equivalent.
         */
        void initWithGroupIndices(ArrayRef<const int> indices)
        {
            mapping_.clear();
            int groupCount = 0;
            if (!indices.empty())
            {
                const int maxIndex = *std::max_element(indices.begin(),
                                                       indices.end());
                mapping_.resize(maxIndex + 1, -1);
                for (int item : indices)
                {
                    GMX_ASSERT(item >= 0, "Negative group numbers not supported");
                    if (mapping_[item] == -1)
                    {
                        mapping_[item] = groupCount;
                        ++groupCount;
                    }
                }
            }
            finder_.init(groupCount);
        }
        /*! \brief
         * Returns a reprensetative value for an item that is unique for each
         * set.
         *
         * `group` should be one of the values that were passed in as an index
         * to initWithGroupIndices().
         * The return value is an internal index that has no simple relation to
         * the input indices.
         */
        int representativeValue(int group)
        {
            GMX_ASSERT(group >= 0 && group < maxGroupNumber(),
                       "Input value out of range");
            GMX_ASSERT(mapping_[group] != -1,
                       "Input value not in initialization set");
            return finder_.representativeItem(mapping_[group]);
        }
        /*! \brief
         * Merges sets that contain two given items.
         *
         * If the items are already in the same set, nothing happens.
         * Each input value should be one of the values that were passed in as
         * an index to initWithGroupIndices().
         */
        void mergeGroups(int group1, int group2)
        {
            GMX_ASSERT(group1 >= 0 && group1 < maxGroupNumber(),
                       "Input value out of range");
            GMX_ASSERT(group2 >= 0 && group2 < maxGroupNumber(),
                       "Input value out of range");
            GMX_ASSERT(mapping_[group1] != -1,
                       "Input value not in initialization set");
            GMX_ASSERT(mapping_[group2] != -1,
                       "Input value not in initialization set");
            finder_.merge(mapping_[group1], mapping_[group2]);
        }
        /*! \brief
         * Returns the sizes of all sets (in arbitrary order).
         *
         * If there were multiple identical indices passed to
         * initWithGroupIndices(), these are only counted as one when
         * computing the sizes.
         */
        std::vector<int> allSizes()
        {
            return finder_.allSizes();
        }

    private:
        int maxGroupNumber() const { return mapping_.size(); }

        UnionFinder        finder_;
        //! Mapping from input indices to zero-based indices used by finder_.
        std::vector<int>   mapping_;
};

} // namespace gmx

#endif
