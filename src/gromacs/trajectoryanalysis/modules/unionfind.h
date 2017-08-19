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

class UnionFinder
{
    public:
        void init(int count)
        {
            parent_.resize(count);
            std::iota(parent_.begin(), parent_.end(), 0);
        }
        void merge(int i, int j)
        {
            const int root = findRoot(i);
            mergeToRoot(j, root);
        }
        std::vector<int> allSizes()
        {
            const int        count = parent_.size();
            std::vector<int> result(count, 0);
            for (int i = 0; i < count; ++i)
            {
                ++result[findRoot(i)];
            }
            result.erase(std::remove(result.begin(), result.end(), 0),
                         result.end());
            return result;
        }

    private:
        int count() const { return parent_.size(); }
        int findRoot(int i)
        {
            GMX_ASSERT(i >= 0 && i < count(), "Input index out of range");
            while (parent_[i] != i)
            {
                const int prev = i;
                i             = parent_[i];
                parent_[prev] = parent_[i];
            }
            return i;
        }
        void mergeToRoot(int i, int root)
        {
            GMX_ASSERT(i >= 0 && i < count(), "Input index out of range");
            GMX_ASSERT(root >= 0 && root < count(), "Internal index out of range");
            do
            {
                const int next = parent_[i];
                parent_[i] = root;
                i          = next;
            }
            while (parent_[i] != root);
        }

        std::vector<int> parent_;
};

class MappedUnionFinder
{
    public:
        void initWithGroupNumbers(ConstArrayRef<int> groupNumbers)
        {
            mapping_.clear();
            int groupCount = 0;
            if (!groupNumbers.empty())
            {
                const int maxIndex = *std::max_element(groupNumbers.begin(),
                                                       groupNumbers.end());
                mapping_.resize(maxIndex + 1, -1);
                for (int g : groupNumbers)
                {
                    GMX_ASSERT(g >= 0, "Negative group numbers not supported");
                    if (mapping_[g] == -1)
                    {
                        mapping_[g] = groupCount;
                        ++groupCount;
                    }
                }
            }
            finder_.init(groupCount);
        }
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
        std::vector<int> allSizes()
        {
            return finder_.allSizes();
        }

    private:
        int maxGroupNumber() const { return mapping_.size(); }

        UnionFinder        finder_;
        std::vector<int>   mapping_;
};

} // namespace gmx

#endif
