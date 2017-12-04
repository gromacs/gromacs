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
 * Declares gmx::ColumnMajorLattice template class.
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \inpublicapi
 * \ingroup module_griddata
 */
#ifndef GMX_MATH_COLUMNMAJORLATTICE
#define GMX_MATH_COLUMNMAJORLATTICE

#include "gromacs/utility/exceptions.h"

#include <algorithm>
#include <array>
#include <numeric>
#include <functional>

namespace gmx
{

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
template <int N> class ColumnMajorLattice
{
    public:

        //! N-dimensional integer index.
        typedef std::array<int, N> MultiIndex;

        /*! \brief
         * Constructs Lattice by setting its extend and checking extend validity.
         *
         * Indices span (0,...,0) to (extend[0]-1,...,extend[N]-1)
         * \TODO check for overflow of number of lattice points
         * \throws gmx::RangeError if any extend is not larger than zero.
         */
        ColumnMajorLattice(const MultiIndex &extend) : extend_ {extend}
        {
            const auto &smallerThanOne = [](int i) { return i < 1; };
            if (std::any_of(std::begin(extend_), std::end(extend_), smallerThanOne))
            {
                GMX_THROW(RangeError("Lattice extend must be larger zero in all dimensions."));
            }
        }
        /*! \brief
         * The lattice extend as N-dimensional index.
         */
        const MultiIndex &extend() const { return extend_; };
        /*! \brief
         * Compare lattices.
         */
        bool operator==(const ColumnMajorLattice &other)
        {
            return extend_ == other.extend_;
        }
        /*! \brief
         * Lattice are unequal if not equal.
         */
        bool operator!=(const ColumnMajorLattice &other)
        {
            return !operator==(other);
        }
        /*! \brief
         * The number of points in the lattice.
         *
         * \returns product of lattice extends
         */
        int getNumLatticePoints() const
        {
            return std::accumulate(std::begin(extend_), std::end(extend_), 1, std::multiplies<int>());
        }

        /*! \brief
         * Returns one-dimensional lattice index from N-dimensional multi index.
         *
         * i_0 + extend_0 * (i_1 + extend_1 * (i_2 + extend_2 * ( ... ))).
         * \throws gmx::RangeError if latticeIndex is not within the lattice.
         * \returns non-negative integer lattice index
         */
        int lineariseVectorIndex(const MultiIndex &latticeIndex) const
        {
            if (!inLattice(latticeIndex))
            {
                GMX_THROW(RangeError("Lattice index must be in lattice: nonnegative "
                                     "and smaller than extend along all dimensions."));
            }

            auto linearIndex            = *latticeIndex.rbegin();
            auto currentDimensionExtend = extend_.rbegin();
            for (auto currentDimensionIndex = ++latticeIndex.rbegin();
                 currentDimensionIndex != latticeIndex.rend();
                 ++currentDimensionIndex)
            {
                ++currentDimensionExtend;
                linearIndex = *currentDimensionIndex + *currentDimensionExtend * linearIndex;
            }

            return linearIndex;
        }

        /*! \brief
         * Generate the multi index from a linear lattice index.
         *
         * The inverse of lineariseVectorIndex.
         *
         * \param[in] linearIndex non-negative and smaller than number of lattice points.
         * \returns the vectorised linear index
         * \throws gmx::RangeError if linearIndex exceeds the number of lattice points.
         */
        MultiIndex vectoriseLinearIndex(int linearIndex) const
        {
            MultiIndex result;
            auto       stride = getNumLatticePoints();
            if (linearIndex < 0  || linearIndex >= stride)
            {
                GMX_THROW(RangeError("Linear index must be non-negative and smaller than number of lattice points."));
            }
            // Build a multi index from a linear index,
            // starting from the slowest moving dimension (from the back)
            // Integer division by the stride yields the vector index,
            // the remainder is used for calculating the remaining vector indices
            auto currentLatticeIndex = result.rbegin();
            for (auto currentDimensionExtend = extend_.rbegin();
                 currentDimensionExtend != extend_.rend(); ++currentDimensionExtend)
            {
                stride              /= *currentDimensionExtend;
                *currentLatticeIndex = linearIndex / stride;
                linearIndex         -= *currentLatticeIndex * stride;
                ++currentLatticeIndex;
            }
            return result;
        }

        /*! \brief
         * Check if latticeIndex is a valid index for this lattice.
         *
         * \param[in] latticeIndex Lattice index to check
         * \returns True if latticeIndex is in lattice
         */
        bool inLattice(const MultiIndex &latticeIndex) const
        {
            auto extendInCurrentDimension = extend_.begin();
            for (const auto &indexInCurrentDimension : latticeIndex)
            {
                if ((indexInCurrentDimension < 0) ||
                    (indexInCurrentDimension >= *extendInCurrentDimension))
                {
                    return false;
                }
                ++extendInCurrentDimension;
            }
            return true;
        }

    private:
        //! \brief the extend of the lattice
        MultiIndex extend_;
};

} // gmx

#endif
