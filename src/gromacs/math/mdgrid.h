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
/*! \internal
 * \file
 * \brief
 * Defines a multidimensional grid.
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \ingroup module_math
 */

#ifndef GMX_MATH_GRID_H
#define GMX_MATH_GRID_H

#include <cmath>

#include <memory>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/iserializer.h"

#include "canonicalvectorbasis.h"
#include "mdarrayindexing.h"

namespace gmx
{

/*!\internal
 * \brief Defines translational coordinate system transformation.
 * Performs a translational coordinate system transformation.
 *
 * The transformation is performed such that a grid with this translation TransformationPolicy
 * has its origin at the translation vector.
 * \tparam Rank The dimension of the coordinate system
 */
template <size_t Rank>
class TranslationTransform
{
    public:
        //! The default constructor initializes the translation to zero.
        constexpr TranslationTransform() : translation_ {}
        {}
        /*! \brief Construct with translation vector.
         * \param[in] translation The translation vector to be used in the transformation.
         */
        TranslationTransform(const MdFloatVector<Rank> &translation) : translation_ {translation}
        {}
        /*! \brief Transform vector into coordinate system.
         * \param[in] x Vector to be translated from external into grid coordinates
         * \returns Translated vector.
         */
        MdFloatVector<Rank> into(const MdFloatVector<Rank> &x) const
        {
            MdFloatVector<Rank> x_translated;
            for (size_t dim = 0; dim < Rank; ++dim)
            {
                x_translated[dim] = x[dim]-translation_[dim];
            }
            return x_translated;
        }
        /*! \brief Transform vector from coordinate system.
         * \param[in] x Vector to be translated from grid into external coordinates
         * \returns the translated vector.
         */
        MdFloatVector<Rank> from(const MdFloatVector<Rank> &x) const
        {
            MdFloatVector<Rank> x_translated;
            for (size_t dim = 0; dim < Rank; ++dim)
            {
                x_translated[dim] = x[dim]+translation_[dim];
            }
            return x_translated;
        }
        /*!\brief Serialize this translation policy.
         *
         * \param[in] serializer The Serializer to read/write the properties of this class.
         */
        void serialize(ISerializer * serializer)
        {
            for (size_t dim = 0; dim < Rank; ++dim)
            {
                serializer->doFloat(&translation_[dim]);
            }
        }

    private:
        //! The external coordinates of (0,...,0)
        MdFloatVector<Rank> translation_ = {};
};

/*! \internal
 * \brief The identiy transformation between vectors.
 * \tparam Rank The dimensionality of the coordinates to be transformed.
 */
template <size_t Rank>
class IdentityTransform
{
    public:
        //! \brief Default constructor.
        constexpr IdentityTransform() = default;
        /*! \brief Implicit conversion to TranslationTransform with zero translation.
         */
        explicit operator TranslationTransform<Rank>() const {return TranslationTransform<Rank>({}); }
        /*! \brief Transform vector into coordinate system.
         * \param[in] x Vector to be translated from external into grid coordinates
         * \returns Translated vector.
         */
        constexpr MdFloatVector<Rank> into(const MdFloatVector<Rank> &x) const
        {
            return x;
        }
        /*! \brief Transform vector from coordinate system.
         * \param[in] x Vector to be translated from grid into external coordinates
         * \returns the translated vector.
         */
        constexpr MdFloatVector<Rank> from(const MdFloatVector<Rank> &x) const
        {
            return x;
        }
        /*! \brief Serialize the translation policy.
         * With no translation there is nothing to serialize
         */
        void serialize(ISerializer * /*serializer*/){}
};

/*!\internal
 * \brief A finite multi-dimensional grid describing regularly spaced points in real space.
 *
 * Addtional vector transformations are governed by transformation policies enabling
 * ,e.g., grids with translations with respect to the origin, rotated grids,
 *  and periodic grids.
 *
 * \tparam Rank The dimension of the grid
 * \tparam TransformationPolicy Coordinate transformations before placing real space vectors on grid.
 */
template <size_t Rank, typename TransformationPolicy = IdentityTransform<Rank> >
class MdGrid
{
    public:
        //! \brief Default constructor generates an empty grid.
        MdGrid() = default;
        /*! \brief Construct grid from its encompasing cell and bounds,
         * with default TranslationPolicy.
         * \param[in] extent The spacial extent of the grid
         * \param[in] bounds The number of grid points in each dimension
         */
        MdGrid(const MdFloatVector<Rank> &extent, const bounds<Rank> &bounds) :
            bounds_ (bounds), unitCell_ (unitCellFromBoundsAndExtent(bounds, extent)), transform_() {}
        /*! \brief
         * Constructs MdGrid from MdFloatVector and bounds.
         * \param[in] extent The spacial extend of the grid
         * \param[in] bounds The number of grid points in each dimension
         * \param[in] transform The coordinate transformation between grid and real space.
         */
        MdGrid(const MdFloatVector<Rank> &extent, const bounds<Rank> &bounds, const TransformationPolicy &transform ) :
            bounds_ (bounds), unitCell_ (unitCellFromBoundsAndExtent(bounds, extent)), transform_(transform) {}

        /*!
         * \brief Compares grids for equality.
         * NOTE: this does not compare if two grids describe the same set of grid points,
         * but rather if two grids have the same bounds, unitCell and transformation.
         */
        bool operator==(const MdGrid<Rank, TransformationPolicy> &other) const
        {
            return (bounds_ == other.bounds && transform_ == other.transform_ && unitCell_ == other.unitCell_);
        }
        /*! \brief
         * Vector pointing from a gridpoint to coordinate in internal grid coordinates.
         * \param[in] x untransformed input coordinate
         * \param[in] i index of bounds point
         * \returns vector from bounds point to coordinate in grid coordinates
         */
        MdFloatVector<Rank> internalVectorFromGridPointToCoordinate(const MdFloatVector<Rank> &x, const offset<Rank> &index) const
        {
            auto                realValuedIndex = unitCell_.transformIntoBasis(transform_.into(x));
            MdFloatVector<Rank> result;
            for (size_t dim = 0; dim < Rank; ++dim)
            {
                result[dim] = realValuedIndex[dim]-index[dim];
            }
            return result;
        }

        /*! \brief
         * The integer part of a point in internal Grid coordinates.
         * Designates the volume element wherein a point lies.
         * o-----o
         * | .   |  . = input vector
         * |/    |  x = returned offset
         * x-----o  o = other ridpoints
         * \param[in] x untransformed input coordinate
         *
         * \returns grid offset of coordinate
         */
        offset<Rank> coordinateToFloorOffset(const MdFloatVector<Rank> &x) const
        {
            const auto   &realValuedIndex = unitCell_.transformIntoBasis(transform_.into(x));
            offset<Rank>  flooredOffset;
            for (size_t dim = 0; dim < Rank; ++dim)
            {
                flooredOffset[dim] = std::floor(realValuedIndex[dim]);
            }
            return flooredOffset;
        }

        /*! \brief
         * The external coordinate of a bounds point.
         *
         * \param[in] i bounds point offset
         *
         * \returns the coordinate cooresponding to a bounds point index by i.
         */
        MdFloatVector<Rank> multiIndexToCoordinate(const offset<Rank> &offset) const
        {
            return transform_.from(unitCell_.transformFromBasis({offset}));
        }

        /*! \brief
         * Resets the bounds and extends or shrinks the containing grid cell correspondingly.
         * The unit cell stays the same, while grid size may change.
         * \param[in] bounds the new bounds
         */
        void setBounds(const bounds<Rank> &bounds)
        {
            bounds_ = bounds;
        }

        /*! \brief
         * Bounds accessor.
         * \returns the bounds
         */
        bounds<Rank> bounds() const
        {
            return bounds_;
        }
        /*! \brief
         * Grid cell accessor.
         * \returns the grid cell comprising all bounds
         */
        CanonicalVectorBasis<Rank> cell() const
        {
            return unitCell_.scaledCopy({bounds_});
        }
        /*! \brief
         * Unit cell accessor.
         * Unit cell times extend yields the grid cell.
         * \returns the unit cell of the grid.
         */
        CanonicalVectorBasis<Rank> unitCell() const
        {
            return unitCell_;
        }
        /*! \brief serializes grid
         * \param[in] serializer Serialzier to write/read grid.
         * TODO: silently fails if bounds::value_type to int64_t conversion fails.
         */
        void serialize(ISerializer * serializer)
        {
            for (size_t dim = 0; dim < Rank; dim++)
            {
                int64_t i;
                if (!serializer->reading())
                {
                    i = bounds_[dim];
                }
                serializer->doInt64(&i);
                if (serializer->reading())
                {
                    bounds_[dim] = i;
                }
            }

            unitCell_.serialize(serializer);
            transform_.serialize(serializer);
        }

    private:
        /*! \brief Evaluate unit cell from cell and grid extend.
         * \param[in] bounds Number of grid points in each dimension
         * \param[in] extent The full extend of the grid
         * \returns Unit cell of the grid
         */
        CanonicalVectorBasis<Rank> unitCellFromBoundsAndExtent(const gmx::bounds<Rank> &bounds, const MdFloatVector<Rank> &extent)
        {
            MdFloatVector<Rank>    unitCellLengths;
            for (size_t dim = 0; dim < Rank; ++dim)
            {
                unitCellLengths[dim] = extent[dim]/static_cast<float>(bounds[dim]);
            }
            return unitCellLengths;
        }
        //! The bounds indexing all points in the grid
        gmx::bounds<Rank>               bounds_;
        //! The pre-calculated grid unit cell.
        CanonicalVectorBasis<Rank>      unitCell_;
        //! Governs how external vectors shall be transformed into the grid.
        TransformationPolicy            transform_;
};

} // namespace gmx

#endif
