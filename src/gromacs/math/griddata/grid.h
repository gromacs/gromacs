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

#include "columnmajorlattice.h"
#include "canonicalvectorbasis.h"

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/iserializer.h"

namespace gmx
{

/*! \brief
 * Interface to a finite N-dimensional grid that relates a finite lattice to
 * N-dimensional space via a basis.
 *
 * \tparam N The dimensionality of the grid
 */
template <size_t Rank>
class IGrid
{
    public:
        //! short-hand notation for N-dimensional vector.
        using NdVector = MdFloatVector<Rank>;
        //! short-hand notation for N-dimensional index.
        using MultiIndex = offset<Rank>;
        //! short-hand notation for reference to grid interface.
        typedef std::unique_ptr < IGrid < Rank>> IGridPointer;
        /*! \brief Compares grids.
         */
        virtual bool operator==(const IGrid<Rank> &other) const = 0;

        /*! \brief
         * Vector pointing from a gridpoint to coordinate in internal grid coordinates.
         *
         * \param[in] x untransformed input coordinate
         * \param[in] i index of lattice point
         * \returns vector from lattice point to coordinate in grid coordinates
         */
        virtual NdVector gridVectorFromGridPointToCoordinate(const NdVector &x, const MultiIndex &i) const = 0;
        /*! \brief
         * The integer part of a point in internal Grid coordinates.
         *
         * \param[in] x untransformed input coordinate
         *
         * \returns lattice point multi index
         */
        virtual MultiIndex coordinateToFloorMultiIndex(const NdVector &x) const     = 0;
        /*! \brief
         * The external coordinate of a lattice point.
         *
         * \param[in] i lattice point multi index
         *
         * \returns the coordinate cooresponding to a lattice point index by i.
         */
        virtual NdVector multiIndexToCoordinate(const MultiIndex &i) const          = 0;

        /*! \brief
         * Resets the lattice and extends or shrinks the containing grid cell correspondingly.
         * The unit cell stays the same, while grid size may change.
         * \param[in] lattice the new lattice
         */
        virtual void setLatticeAndRescaleCell(const bounds<Rank> &lattice) = 0;

        /*! \brief
         * Lattice accessor.
         * \returns the lattice
         */
        virtual bounds<Rank> lattice() const      = 0;

        /*! \brief
         * grid cell accessor.
         * \returns the grid cell comprising all lattice
         */
        virtual CanonicalVectorBasis<Rank> cell() const       = 0;

        /*! \brief
         * Unit cell accessor.
         * Unit cell times extend yields the grid cell.
         * \returns the unit cell of the grid.
         */
        virtual CanonicalVectorBasis<Rank> unitCell() const   = 0;
        /*! \brief
         * Return an owning pointer to the abstract base class of a copy of a specific implementation.
         *
         * Enables copying of classes holding unique_ptr's to this interface.
         *
         * \returns IGridPointer to a copy of this.
         */
        virtual IGridPointer duplicate() const = 0;

        virtual void readwrite(ISerializer * serializer) = 0;

};


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

/*! \brief
 * A finite N-dimensional grid that is anchored at the coordinate system origin.
 *
 * Relation of the integer indices,i, to real space vector, r, are given by the unit cell.
 */
template <size_t Rank>
class Grid : public IGrid<Rank>
{
    public:
        using NdVector   = MdFloatVector<Rank>;
        using MultiIndex = offset<Rank>;
        /*! \brief Default constructor.
         */
        Grid() = default;
        /*! \brief
         * Constructs Grid from CanonicalVectorBasis and ColumnMajorLattice.
         * \param[in] cell The spacial extend of the grid
         * \param[in] lattice The lattice that describes the grid points.
         */
        Grid(const CanonicalVectorBasis<Rank> &cell, const bounds<Rank> &lattice) : IGrid<Rank>(), cell_ (cell), lattice_ (lattice), unitCell_ (calculateUnitCellFromLatticeAndCell(lattice, cell)) {}
        /*! \brief
         * Copy constructor declared, because the default constructor is deleted.
         * \param[in] other Grid to be copied from
         */
        Grid(const Grid &other) : cell_ (other.cell_), lattice_ (other.lattice_), unitCell_ (other.unitCell_) {}

        /*! \brief
         * Copy assignment operator declared, because the default constructor is deleted.
         * \param[in]
         * \returns Copied Grid
         */
        Grid<Rank> &operator=(Grid<Rank> other)
        {
            std::swap(cell_, other.cell_);
            std::swap(unitCell_, other.unitCell_);
            std::swap(lattice_, other.lattice_);
            return *this;
        }
        //! \copydoc IGrid::operator==(const IGrid<Rank> & other)
        bool operator==(const IGrid<Rank> &other) const override
        {
            if (other.lattice() != this->lattice_)
            {
                return false;
            }
            if (lattice_ != other.lattice())
            {
                return false;
            }
            MultiIndex nought;
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
        //! \copydoc IGrid::gridVectorFromGridPointToCoordinate(const NdVector &x, const MultiIndex & i)
        NdVector gridVectorFromGridPointToCoordinate(const NdVector &x, const MultiIndex &index) const override
        {
            auto     realValuedIndex = unitCell_.transformIntoBasis(x);
            NdVector result;
            for (size_t dim = 0; dim < Rank; ++dim)
            {
                result[dim] = realValuedIndex[dim]-index[dim];
            }
            return result;
        };

        //! \copydoc IGrid::coordinateToFloorMultiIndex(const NdVector &x)
        MultiIndex coordinateToFloorMultiIndex(const NdVector &x) const override
        {
            const auto &realValuedIndex = unitCell_.transformIntoBasis(x);
            MultiIndex  flooredMultiIndex;
            for (size_t dim = 0; dim < Rank; ++dim)
            {
                flooredMultiIndex[dim] = std::floor(realValuedIndex[dim]);
            }
            return flooredMultiIndex;
        }

        //! \copydoc IGrid::multiIndexToCoordinate(const MultiIndex &i)
        NdVector multiIndexToCoordinate(const MultiIndex &i) const override
        {
            return unitCell_.transformFromBasis(multiIndexToNdVector(i));
        }

        //! \copydoc IGrid::setLatticeAndRescaleCell(const ColumnMajorLattice<Rank> &lattice)
        void setLatticeAndRescaleCell(const bounds<Rank> &lattice) override
        {
            lattice_ = lattice;
            cell_    = unitCell_.scaledCopy(multiIndexToNdVector(lattice_));
        }

        //! \copydoc IGrid::lattice()
        bounds<Rank> lattice() const override
        {
            return lattice_;
        }

        //! \copydoc IGrid::cell()
        CanonicalVectorBasis<Rank> cell() const override
        {
            return cell_;
        }

        //! \copydoc IGrid::unitCell()
        CanonicalVectorBasis<Rank> unitCell() const override
        {
            return unitCell_;
        };

        //! \copydoc IGrid::duplicate()
        typename IGrid<Rank>::IGridPointer duplicate() const override
        {
            return typename IGrid<Rank>::IGridPointer(new Grid<Rank>(*this));
        }

        //! \copydoc IGrid::readwrite(ISerializer * serializer)
        void readwrite(ISerializer * serializer) override
        {
            serialize(serializer, lattice_);
            cell_.readwrite(serializer);
            unitCell_.readwrite(serializer);
        }

    private:

        /*! \brief
         * Cast MultiIndex to NdVector.
         * \param[in] i N-dimensional index.
         */
        NdVector multiIndexToNdVector(const offset<Rank> &index) const
        {
            NdVector ndVector;

            for (size_t dim = 0; dim < Rank; ++dim)
            {
                ndVector[dim] = index[dim]; // employ implicit casting
            }
            return ndVector;
        }

        /*! \brief
         * Cast MultiIndex to NdVector.
         * \param[in] i N-dimensional index.
         */
        NdVector multiIndexToNdVector(const bounds<Rank> &index) const
        {
            NdVector ndVector;

            for (size_t dim = 0; dim < Rank; ++dim)
            {
                ndVector[dim] = index[dim]; // employ implicit casting
            }
            return ndVector;
        }

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
};

/*! \brief
 * An N-dimensional grid that is shifted from the origin.
 *
 * \tparam N dimensionality of the grid
 */
template <size_t Rank>
class GridWithTranslation : public Grid<Rank>
{
    public:
        using NdVector   = MdFloatVector<Rank>;
        using MultiIndex = offset<Rank>;

        /*! \brief
         * Default construct \ref gmx::Grid with additional translation vector.
         */
        GridWithTranslation() = default;

        /*! \brief
         * Construct \ref gmx::Grid with additional translation vector.
         */
        GridWithTranslation(const CanonicalVectorBasis<Rank> &cell, const bounds<Rank> &lattice, const NdVector &translation) : Grid<Rank>(cell, lattice), translation_ (translation) { }

        /*! \brief
         * Copy constructor declared, because the default constructor is deleted.
         * \param[in] other Grid to be copied from
         */
        GridWithTranslation(const GridWithTranslation &other) : Grid<Rank>(other), translation_ (other.translation_) { }

        /*! \brief
         * Copy assignment operator declared, because the default constructor is deleted.
         * \param[in]
         * \returns Copied Grid
         */
        GridWithTranslation<Rank> &operator=(GridWithTranslation<Rank> other)
        {
            Grid<Rank>::operator=(other);
            std::swap(translation_, other.translation_);
            return *this;
        }

        //! \copydoc IGrid::gridVectorFromGridPointToCoordinate(const NdVector &x, const MultiIndex & i)
        NdVector gridVectorFromGridPointToCoordinate(const NdVector &x, const MultiIndex &i) const override
        {
            return Grid<Rank>::gridVectorFromGridPointToCoordinate(translateIntoGrid_(x), i);
        };

        //! \copydoc IGrid::coordinateToFloorMultiIndex(const NdVector &x)
        MultiIndex coordinateToFloorMultiIndex(const NdVector &x) const override
        {
            return Grid<Rank>::coordinateToFloorMultiIndex(translateIntoGrid_(x));
        }

        //! \copydoc IGrid::multiIndexToCoordinate(const MultiIndex &i)
        NdVector multiIndexToCoordinate(const MultiIndex &i) const override
        {
            return translateFromGrid_(Grid<Rank>::multiIndexToCoordinate(i));
        }

        /*! \brief
         * Set the real-space coordinate of gridpoint (0,0,0).
         * \param[in] translate Translation vector
         */
        void setTranslation(const NdVector &translate)
        {
            translation_ = translate;
        }

        //! \copydoc IGrid::duplicate()
        typename IGrid<Rank>::IGridPointer duplicate() const override
        {
            return typename IGrid<Rank>::IGridPointer(new GridWithTranslation<Rank>(*this));
        }

        //! \copydoc IGrid::readwrite(const ISerializer * serializer)
        void readwrite(ISerializer * serializer) override
        {
            Grid<Rank>::readwrite(serializer);
            for (size_t dim = 0; dim < Rank; ++dim)
            {
                serializer->doFloat(&translation_[dim]);
            }
        }

    private:
        /*! \brief
         * \param[in] x Vector to be translated from external into grid coordinates
         * \returns Translated vector.
         */
        NdVector translateIntoGrid_(const NdVector &x) const
        {
            NdVector x_translated;
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
        NdVector translateFromGrid_(const NdVector &x) const
        {
            NdVector x_translated;
            for (size_t dim = 0; dim < Rank; ++dim)
            {
                x_translated[dim] = x[dim]+translation_[dim];
            }
            return x_translated;
        };
        //! The external coordinates of lattice point (0,...,0)
        NdVector                     translation_;
};


}

#endif
