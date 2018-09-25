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
/*! \internal \file
 * \brief
 * Defines canonical vector basis.
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \ingroup module_math
 */

#ifndef GMX_MATH_CanonicalVectorBasis_H
#define GMX_MATH_CanonicalVectorBasis_H

#include <cmath>

#include <algorithm>
#include <array>
#include <functional>
#include <numeric>

#include "gromacs/utility/iserializer.h"

#include "mdarrayindexing.h"

namespace gmx
{

/*! \internal
 * \brief Multidimensional float vector.
 * Points to positions in multidimensional grids.
 * \tparam Rank The dimensionality of the vector
 */
template <size_t Rank>
class MdFloatVector
{
    public:
        //! \brief Publicly expose the rank.
        static constexpr size_t rank = Rank;
        //! \brief Reference to data elements.
        using reference              = float &;
        //! \brief Const reference to data elements.
        using const_reference        = const float &;
        //! \brief The size type for the vector.
        using size_type              = size_t;
        //! \brief Value type of the float vector.
        using value_type             = float;
        //! \brief iterator type
        using iterator = typename std::array<value_type, Rank>::iterator;
        //! \brief const iterator type
        using const_iterator = typename std::array<value_type, Rank>::const_iterator;

        /*! \brief Default constructor.
         */
        constexpr MdFloatVector() noexcept = default;
        /*! \brief Construct mdFloatVector from plain C array.
         * \param[in] mdFloatVector The mdFloatVector as plain C array.
         *
         * Allows conversion from gmx::RVec for three dimensional grids.
         */
        MdFloatVector(const value_type (&mdFloatVector)[Rank])
        {
            std::copy(std::begin(mdFloatVector), std::end(mdFloatVector), std::begin(mdFloatVector_));
        }
        /*! \brief Construct mdFloatVector from std::array.
         * \param[in] mdFloatVector The mdFloatVector as array.
         */
        MdFloatVector(const std::array<value_type, Rank> &mdFloatVector)
        {
            mdFloatVector_ = mdFloatVector;
        }
        /*! \brief Construct mdFloatVector from an offset.
         * \param[in] offset The array offset.
         *
         * Facilitates conversions from offsets to real space vectors.
         */
        MdFloatVector(const offset<Rank> &offset)
        {
            for (size_t dim = 0; dim < Rank; ++dim)
            {
                mdFloatVector_[dim] = offset[dim];
            }
        }
        /*! \brief Construct mdFloatVector from bounds.
         * \param[in] bounds The array bounds.
         *
         * Facilitates conversions from bounds to real space vectors, needed e.g.
         * in the construction of reciprocal fourier grids.
         */
        MdFloatVector(const bounds<Rank> &bounds)
        {
            for (size_t dim = 0; dim < Rank; ++dim)
            {
                mdFloatVector_[dim] = bounds[dim];
            }
        }
        /*! \brief begin iterator
         * \returns iterator to begin of multidimensional float vector.
         */
        iterator begin(){return mdFloatVector_.begin(); }
        /*! \brief begin iterator
         * \returns iterator to begin of multidimensional float vector.
         */
        const_iterator begin() const {return mdFloatVector_.begin(); }
        /*! \brief end iterator
         * \returns iterator to end of multidimensional float vector.
         */
        iterator end(){return mdFloatVector_.end(); }
        /*! \brief const end iterator
         * \returns const iterator to end of multidimensional float vector.
         */
        const_iterator end() const {return mdFloatVector_.end(); }

        /*! \brief Provide the size of the vector
         * \returns The dimension of the mdFloatVector.
         */
        constexpr size_type size() const {return Rank; }
        /*! \brief Access the constant reference to the n-th mdFloatVector.
         * \param[in] n The dimension of the mdFloatVector.
         */
        reference operator[](size_type n) { return mdFloatVector_[n]; }
        /*! \brief Access the constant reference to the n-th mdFloatVector.
         * \param[in] n The dimension of the mdFloatVector.
         */
        constexpr const_reference operator[](size_type n) const { return mdFloatVector_[n]; }
    private:
        //!\brief The mdFloatVector.
        std::array<value_type, Rank> mdFloatVector_ = {};
};

/*!\brief Check if two multidimensional vectors are equal.
 * \tparam Rank the dimension of the multidimensional vectors.
 * \param[in] lhs first multidimensional vector to compare
 * \param[in] rhs second multidimensional vector to compare
 * \returns true if two multidimensional vectors are equal.
 */
template <size_t Rank>
bool operator==(const MdFloatVector<Rank> &lhs, const MdFloatVector<Rank> &rhs) noexcept
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

/*!\brief Check if two multidimensional vectors are not equal.
 * \tparam Rank the dimension of the multidimensional vectors.
 * \param[in] lhs first multidimensional vector to compare
 * \param[in] rhs second multidimensional vector to compare
 * \returns true if two multidimensional vectors are not equal.
 */
template <size_t Rank>
constexpr bool operator!=(const MdFloatVector<Rank> &lhs, const MdFloatVector<Rank> &rhs) noexcept
{
    return !(lhs == rhs);
};

/*! \internal
 * \brief The canonical basis of an Rank-dimensional vector space.
 *
 * The canoncial vector basis is orthogonal and aligned with the axis.
 * It transforms vectors between coordinate systems.
 *
 * \tparam size_t Rank Number of dimensions of the vector space
 */
template <size_t Rank>
class CanonicalVectorBasis
{
    public:
        //! \brief Defaults to the construction of the canonical basis.
        CanonicalVectorBasis()
        {
            for (size_t dim = 0; dim < Rank; ++dim)
            {
                basisVectorLengths_[dim]        = 1;
                basisVectorLengthsInverse_[dim] = 1;
            }
        }
        /*! \brief
         * Construct a basis from the length of the unit vectors.
         * \param[in] basisVectorLengths length of the basis vectors.
         */
        CanonicalVectorBasis(const MdFloatVector<Rank> &basisVectorLengths)
        {
            basisVectorLengths_ = basisVectorLengths;
            for (size_t dim = 0; dim < Rank; ++dim)
            {
                basisVectorLengthsInverse_[dim] = 1./basisVectorLengths_[dim];
            }
        }
        /*! \brief
         * Construct a basis from the length of the unit vectors.
         * \param[in] basisVectorLengths length of the basis vectors.
         */
        CanonicalVectorBasis(const float (&basisVectorLengths)[Rank])
        {
            basisVectorLengths_ = basisVectorLengths;
            for (size_t dim = 0; dim < Rank; ++dim)
            {
                basisVectorLengthsInverse_[dim] = 1./basisVectorLengths_[dim];
            }
        }
        /*! \brief
         * Convert N-dimensional vector into basis representation.
         *
         * \param[in] x N-dimenisonal vector.
         * \returns N-dimensional vector in this basis representation.
         */
        MdFloatVector<Rank> transformIntoBasis(const MdFloatVector<Rank> &x)  const
        {
            MdFloatVector<Rank> result;
            for (size_t dim = 0; dim < Rank; ++dim)
            {
                result[dim] = basisVectorLengthsInverse_[dim] * x[dim];
            }
            return result;
        }
        /*! \brief
         * Convert N-dimensional vector from this basis representation into unit basis.
         *
         * \param[in] x N-dimenisonal vector.
         * \returns N-dimensional vector in unit representation.
         */
        MdFloatVector<Rank> transformFromBasis(const MdFloatVector<Rank> &x)  const
        {
            MdFloatVector<Rank> result;
            for (size_t dim = 0; dim < Rank; ++dim)
            {
                result[dim] = basisVectorLengths_[dim] * x[dim];
            }
            return result;
        }
        /*! \brief
         * The inverse of this basis.
         *
         * Swaps the effect of transformFromBasis and transformIntoBasis.
         * \returns inverted basis.
         */
        CanonicalVectorBasis<Rank> inverse() const
        {
            return CanonicalVectorBasis(basisVectorLengthsInverse_);
        }

        /*! \brief
         * Returns a basis with all basis vector entries multiplied by the respective scale vector entries.
         *
         * \param[in] scale Factor to multipy each basis vector with
         * \returns Scaled basis
         */
        CanonicalVectorBasis<Rank> scaledCopy(const MdFloatVector<Rank> &scale) const
        {
            MdFloatVector<Rank> scaledCell;
            for (size_t dim = 0; dim < Rank; ++dim)
            {
                scaledCell[dim] = scale[dim] * basisVectorLengths_[dim];
            }
            return CanonicalVectorBasis(scaledCell);
        }
        /*! \brief
         * The volume that is enclosed the basis vector parallelepiped.
         *
         * \returns the volume
         */
        float volume() const
        {
            float volume = 1;
            for (size_t dim = 0; dim < Rank; ++dim)
            {
                volume *= basisVectorLengths_[dim];
            }
            return volume;
        }
        /*! \brief
         * The length of all basis vectors.
         *
         * \returns a vector with the lengths
         */
        MdFloatVector<Rank> basisVectorLengths() const
        {
            return basisVectorLengths_;
        }
        /*! \brief
         * The length of a basis vector in the respective dimension.
         *
         * \param[in] dimension the dimension of the basis vector of which the length is requested
         * \returns the length
         */
        float operator[](const int &dimension) const
        {
            return basisVectorLengths_[dimension];
        }
        /*! \brief Serialize to other nodes.
         * \param[in] serializer pointer to serializing interface
         */
        void serialize(ISerializer * serializer)
        {
            for (size_t dim = 0; dim < Rank; ++dim)
            {
                serializer->doFloat(&basisVectorLengths_[dim]);
            }
            for (size_t dim = 0; dim < Rank; ++dim)
            {
                serializer->doFloat(&basisVectorLengthsInverse_[dim]);
            }
        }

    protected:
        //! The Lengths of the basis vector.
        MdFloatVector<Rank> basisVectorLengths_;
        //! Pre-calcualted inverse of basis vector length.
        MdFloatVector<Rank> basisVectorLengthsInverse_;
};

} // namespace gmx
#endif
