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
/*!  \file
 * \brief
 * Defines an N-dimensional grid.
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \inpublicapi
 */

#ifndef GMX_MATH_FINITEGRID_H
#define GMX_MATH_FINITEGRID_H

#include <cmath>

#include <memory>

#include "gromacs/math/vectypes.h"
#include "gromacs/math/griddata/canonicalvectorbasis.h"
#include "gromacs/math/griddata/mdarrayindexing.h"
#include "gromacs/utility/iserializer.h"

namespace gmx
{

/*!\brief serializes md array bounds
 * TODO: fails if static cast of ptrdiff_t to int fails on any bound dimension
 */
template <size_t Rank>
void serialize(ISerializer * serializer, bounds<Rank> bounds)
{
    for (size_t dim = 0; dim < Rank; dim++)
    {
        int64_t i;
        if (!serializer->reading())
        {
            i = bounds[dim];
        }
        serializer->doInt64(&i);
        if (serializer->reading())
        {
            bounds[dim] = i;
        }
    }
}

template <size_t Rank>
class GridWithTranslation
{
    public:
        constexpr GridWithTranslation() : translation_ {}
        {}
        GridWithTranslation(const MdFloatVector<Rank> &translation) : translation_ {translation}
        {}
        /*! \brief
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
        };
        /*! \brief
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
        };
        void serialize(ISerializer * serializer)
        {
            for (size_t dim = 0; dim < Rank; ++dim)
            {
                serializer->doFloat(&translation_[dim]);
            }
        }

    private:
        //! The external coordinates of lattice point (0,...,0)
        MdFloatVector<Rank> translation_ = {};
};

template <size_t Rank>
class GridWithoutTranslation
{
    public:
        constexpr GridWithoutTranslation(){};
        constexpr GridWithoutTranslation(const MdFloatVector<Rank> &translation){};
        explicit operator GridWithTranslation<Rank>() const {return GridWithTranslation<Rank>({}); };
        /*! \brief
         * \param[in] x Vector to be translated from external into grid coordinates
         * \returns Translated vector.
         */
        constexpr MdFloatVector<Rank> into(const MdFloatVector<Rank> &x) const
        {
            return x;
        };
        /*! \brief
         * \param[in] x Vector to be translated from grid into external coordinates
         * \returns the translated vector.
         */
        constexpr MdFloatVector<Rank> from(const MdFloatVector<Rank> &x) const
        {
            return x;
        };
        /*! \brief Serialize the translation policy.
         * With no translation there is nothing to serialize
         */
        void serialize(ISerializer * serializer){}
};


/*! \brief
 * A finite N-dimensional grid that is anchored at the coordinate system origin.
 *
 * Relation of the integer indices,i, to real space vector, r, are given by the unit cell.
 */
template <size_t Rank, typename TranslationPolicy = GridWithoutTranslation<Rank> >
class Grid
{
    public:
        //! short-hand notation for N-dimensional vector.
        using NdVector   = MdFloatVector<Rank>;

        /*! \brief Default constructor.
         */
        Grid() = default;
        Grid(const CanonicalVectorBasis<Rank> &cell, const bounds<Rank> &lattice) :
            cell_(cell), lattice_ (lattice), unitCell_ (calculateUnitCellFromLatticeAndCell(lattice, cell)), translator_() {}

        /*! \brief
         * Constructs Grid from CanonicalVectorBasis and bounds.
         * \param[in] cell The spacial extend of the grid
         * \param[in] lattice The lattice that describes the grid points.
         */
        Grid(const CanonicalVectorBasis<Rank> &cell, const bounds<Rank> &lattice, const TranslationPolicy &translator ) :
            cell_(cell), lattice_ (lattice), unitCell_ (calculateUnitCellFromLatticeAndCell(lattice, cell)), translator_(translator) {}
        /*! \brief
         * Copy constructor declared, because the default constructor is deleted.
         * \param[in] other Grid to be copied from
         */
        Grid(const Grid &other) : cell_ (other.cell_), lattice_ (other.lattice_), unitCell_ (other.unitCell_), translator_(other.translator_) {}

        /*! \brief
         * Copy assignment operator declared, because the default constructor is deleted.
         * \param[in]
         * \returns Copied Grid
         */
        Grid<Rank, TranslationPolicy> &operator=(Grid<Rank, TranslationPolicy> other)
        {
            std::swap(cell_, other.cell_);
            std::swap(unitCell_, other.unitCell_);
            std::swap(lattice_, other.lattice_);
            std::swap(translator_, other.translator_);
            return *this;
        }
        //! \brief compares grids for equality
        bool operator==(const Grid<Rank, TranslationPolicy> &other) const
        {
            if (other.lattice() != this->lattice_)
            {
                return false;
            }
            if (lattice_ != other.lattice())
            {
                return false;
            }
            offset<Rank> nought {};
            if (multiIndexToCoordinate(nought) != other.multiIndexToCoordinate(nought))
            {
                return false;
            }
            for (size_t dimension = 0; dimension < Rank; ++dimension)
            {
                auto basisVector = nought;
                basisVector[dimension] = 1;
                if (lattice_[dimension] > 1 && (multiIndexToCoordinate(basisVector) != other.multiIndexToCoordinate(basisVector)))
                {return false; }
            }
            return true;
        }
        /*! \brief
         * Vector pointing from a gridpoint to coordinate in internal grid coordinates.
         *
         * \param[in] x untransformed input coordinate
         * \param[in] i index of lattice point
         * \returns vector from lattice point to coordinate in grid coordinates
         */
        NdVector gridVectorFromGridPointToCoordinate(const NdVector &x, const offset<Rank> &index) const
        {
            auto     realValuedIndex = unitCell_.transformIntoBasis(translator_.into(x));
            NdVector result;
            for (size_t dim = 0; dim < Rank; ++dim)
            {
                result[dim] = realValuedIndex[dim]-index[dim];
            }
            return result;
        };

        /*! \brief
         * The integer part of a point in internal Grid coordinates.
         *
         * \param[in] x untransformed input coordinate
         *
         * \returns lattice point multi index
         */
        offset<Rank> coordinateToFloorMultiIndex(const NdVector &x) const
        {
            const auto   &realValuedIndex = unitCell_.transformIntoBasis(translator_.into(x));
            offset<Rank>  flooredMultiIndex;
            for (size_t dim = 0; dim < Rank; ++dim)
            {
                flooredMultiIndex[dim] = std::floor(realValuedIndex[dim]);
            }
            return flooredMultiIndex;
        }

        /*! \brief
         * The external coordinate of a lattice point.
         *
         * \param[in] i lattice point multi index
         *
         * \returns the coordinate cooresponding to a lattice point index by i.
         */
        NdVector multiIndexToCoordinate(const offset<Rank> &i) const
        {
            return translator_.from(unitCell_.transformFromBasis({i}));
        }

        /*! \brief
         * Resets the lattice and extends or shrinks the containing grid cell correspondingly.
         * The unit cell stays the same, while grid size may change.
         * \param[in] lattice the new lattice
         */
        void setLatticeAndRescaleCell(const bounds<Rank> &lattice)
        {
            lattice_ = lattice;
            cell_    = unitCell_.scaledCopy({lattice_});
        }


        /*! \brief
         * Lattice accessor.
         * \returns the lattice
         */
        bounds<Rank> lattice() const
        {
            return lattice_;
        }

        /*! \brief
         * grid cell accessor.
         * \returns the grid cell comprising all lattice
         */CanonicalVectorBasis<Rank> cell() const
        {
            return cell_;
        }

        /*! \brief
         * Unit cell accessor.
         * Unit cell times extend yields the grid cell.
         * \returns the unit cell of the grid.
         */
        CanonicalVectorBasis<Rank> unitCell() const
        {
            return unitCell_;
        };

        //! \brief serializes grid
        void serialize(ISerializer * serializer)
        {
            gmx::serialize<Rank>(serializer, lattice_);
            cell_.serialize(serializer);
            unitCell_.serialize(serializer);
            translator_.serialize(serializer);
        }

    private:
        /*! \brief
         * Re-evaluate unit cell from cell and grid extend
         * \param[in] lattice The lattice, the grid points are based on
         * \param[in] cell The full extend of the grid
         * \returns Unit cell of the grid
         */
        CanonicalVectorBasis<Rank> calculateUnitCellFromLatticeAndCell(const bounds<Rank> &lattice, const CanonicalVectorBasis<Rank> cell)
        {
            NdVector    cellToUnitCellScale;
            for (size_t dim = 0; dim < Rank; ++dim)
            {
                cellToUnitCellScale[dim] = 1.0/static_cast<real>(lattice[dim]);
            }
            return cell.scaledCopy(cellToUnitCellScale);
        }
        //! The cell spanning the whole grid extend
        CanonicalVectorBasis<Rank>    cell_;
        //! The lattice indexing all points in the grid
        bounds<Rank>                  lattice_;
        //! The pre-calculated grid unit cell.
        CanonicalVectorBasis<Rank>    unitCell_;
        TranslationPolicy             translator_;
};

}

#endif
