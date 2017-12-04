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
/*!  \file
 * \brief
 * Defines canonical vector basis.
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \inpublicapi
 */

#ifndef GMX_MATH_CanonicalVectorBasis_H
#define GMX_MATH_CanonicalVectorBasis_H
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/real.h"
#include <array>
#include <numeric>
#include <algorithm>
#include <functional>
#include <cmath>

namespace gmx
{

/*! \brief
 * The canonical basis of an N-dimensional vector space.
 *
 * The canoncial vector basis is orthogonal and aligned with the axis.
 * It transforms vectors into and out of the
 * After an identy mapping from vector to vector, this is the next simple
 * implementation of a vector space basis.
 *
 *
 * \tparam int Number of dimensions of the vector space
 */
template <int N>
class CanonicalVectorBasis
{
    public:
        //! An N-dimensional real vector.
        typedef std::array<real, N> NdVector;

        /*! \brief
         * Construct a basis from the length of the unit vectors.
         * \param[in] basisVectorLengths length of the basis vectors.
         */
        CanonicalVectorBasis(const NdVector &basisVectorLengths)
        {
            basisVectorLengths_ = basisVectorLengths;
            const auto &isZero = [](real length){return fabs(length) <= GMX_REAL_EPS; };
            if (std::any_of(std::begin(basisVectorLengths_), std::end(basisVectorLengths_), isZero))
            {
                GMX_THROW(RangeError("Basis vector length cannot be zero."));
            }
            std::transform(std::begin(basisVectorLengths_), std::end(basisVectorLengths_), std::begin(basisVectorLengthsInverse_), [](real c){return 1/c; });
        }

        /*! \brief
         * Convert N-dimensional vector into basis representation.
         *
         * \param[in] x N-dimenisonal vector.
         * \returns N-dimensional vector in this basis representation.
         */
        NdVector transformIntoBasis(const NdVector &x)  const
        {
            NdVector result;
            std::transform(std::begin(basisVectorLengthsInverse_), std::end(basisVectorLengthsInverse_), std::begin(x), std::begin(result), std::multiplies<real>());
            return result;
        }

        /*! \brief
         * Convert N-dimensional vector from this basis representation into unit basis.
         *
         * \param[in] x N-dimenisonal vector.
         * \returns N-dimensional vector in unit representation.
         */
        NdVector transformFromBasis(const NdVector &x)  const
        {
            NdVector result;
            std::transform(std::begin(basisVectorLengths_), std::end(basisVectorLengths_), std::begin(x), std::begin(result), std::multiplies<real>());
            return result;
        }

        /*! \brief
         * The inverse of this basis.
         *
         * Swaps the effect of transformFromBasis and transformIntoBasis.
         * \returns inverted basis.
         */
        CanonicalVectorBasis<N> inverse() const
        {
            return CanonicalVectorBasis(basisVectorLengthsInverse_);
        }

        /*! \brief
         * Returns a basis with all basis vector entries multiplied by the respective scale vector entries.
         *
         * \param[in] scale Factor to multipy each basis vector with
         * \returns Scaled basis
         */
        CanonicalVectorBasis<N> scaledCopy(const NdVector &scale) const
        {
            NdVector scaledCell;
            std::transform(std::begin(basisVectorLengths_), std::end(basisVectorLengths_), std::begin(scale), std::begin(scaledCell), std::multiplies<real>());
            return CanonicalVectorBasis(scaledCell);
        }

        /*! \brief
         * The volume that is enclosed the basis vector parallelepiped.
         *
         * \returns the volume
         */
        real volume() const
        {
            return std::accumulate(std::begin(basisVectorLengths_), std::end(basisVectorLengths_), 1.0, std::multiplies<real>());
        }

        /*! \brief
         * The length of all basis vectors.
         *
         * \returns a vector with the lengths
         */
        NdVector basisVectorLengths() const
        {
            return basisVectorLengths_;
        }


        /*! \brief
         * The length of a basis vector in the specified dimension.
         *
         * \param[in] dimension the dimension of the basis vector of which the length is requested
         * \throws gmx::RangeError if dimension negative or larger than dimensionality of basis
         * \returns the length
         */
        real basisVectorLength(const int &dimension) const
        {
            if (dimension >= N || dimension < 0)
            {
                GMX_THROW(RangeError("Dimension of basis vector must be non-negative and smaller than Basis dimension."));
            }
            return basisVectorLengths_[dimension];
        }

    protected:
        //! The Lengths of the basis vector.
        NdVector basisVectorLengths_;
        //! Pre-calcualted inverse of basis vector length.
        NdVector basisVectorLengthsInverse_;
};

} // namespace gmx
#endif
