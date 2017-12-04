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
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1DIM01  USA.
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
#ifndef GMX_MATH_ROTATEDGRID_H
#define GMX_MATH_ROTATEDGRID_H

#include "grid.h"
#include "gromacs/math/quaternion.h"
namespace gmx
{


/*! \brief
 * An N-dimensional grid that is shifted from the origin.
 *
 * \tparam N dimensionality of the grid
 */

class GridWithTranslationOrientation : public Grid<DIM>
{
    public:
        typedef typename CanonicalVectorBasis<DIM>::NdVector NdVector;
        typedef typename ColumnMajorLattice<DIM>::MultiIndex MultiIndex;

        /*! \brief
         * Construct \ref gmx::Grid with additional translation and orientation.
         */
        GridWithTranslationOrientation(const CanonicalVectorBasis<DIM> &cell, const ColumnMajorLattice<DIM> &lattice, const NdVector &translation, const Quaternion orientation ) :
            Grid<DIM>{cell, lattice}, translation_ {translation}, orientation_ {orientation} { }

        /*! \brief
         * Infer translation and orientation from other grid properties
         */
        GridWithTranslationOrientation(const IGrid<DIM> &other ) : Grid<DIM>{other.cell(), other.lattice()}, orientation_ {{0., 0., 1.}, 0.}
        { translation_ = other.multiIndexToCoordinate({{0, 0, 0}}); }


        /*! \brief
         * Copy constructor declared, because the default constructor is deleted.
         * \param[in] other Grid to be copied from
         */
        GridWithTranslationOrientation(const GridWithTranslationOrientation &other) : Grid<DIM>{other}, translation_ {other.translation_}, orientation_ {other.orientation_} {}

        /*! \brief
         * Copy assignment operator declared, because the default constructor is deleted.
         * \param[in]
         * \returns Copied Grid
         */
        GridWithTranslationOrientation &operator=(GridWithTranslationOrientation other)
        {
            Grid<DIM>::operator=(other);
            std::swap(translation_, other.translation_);
            std::swap(orientation_, other.orientation_);
            return *this;
        }

        //! \copydoc IGrid::gridVectorFromGridPointToCoordinate(const NdVector &x, const MultiIndex & i)
        NdVector gridVectorFromGridPointToCoordinate(const NdVector &x, const MultiIndex &i) const override
        {
            return Grid<DIM>::gridVectorFromGridPointToCoordinate(translateIntoGrid_(x), i);
        };

        //! \copydoc IGrid::coordinateToFloorMultiIndex(const NdVector &x)
        MultiIndex coordinateToFloorMultiIndex(const NdVector &x) const override
        {
            return Grid<DIM>::coordinateToFloorMultiIndex(translateIntoGrid_(x));
        }

        //! \copydoc IGrid::multiIndexToCoordinate(const MultiIndex &i)
        NdVector multiIndexToCoordinate(const MultiIndex &i) const override
        {
            return translateFromGrid_(Grid<DIM>::multiIndexToCoordinate(i));
        }
        /*! \brief
         * Set the orientation around the grid mid-point .
         * \param[in] orientation Quaternion
         */
        void setOrientation(const Quaternion &orientation)
        {
            orientation_ = orientation;
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
        typename IGrid<DIM>::IGridPointer duplicate() const override
        {
            return typename IGrid<DIM>::IGridPointer(new GridWithTranslationOrientation(*this));
        }

    private:
        /*! \brief
         * \param[in] x Vector to be translated from external into grid coordinates
         * \returns Translated vector.
         */
        NdVector translateIntoGrid_(const NdVector &x) const
        {
            CanonicalVectorBasis<DIM>::NdVector xInternal;

            std::transform(std::begin(x), std::end(x), std::begin(translation_), std::begin(xInternal), std::minus<real>());

            auto gridCenter = this->cell().scaledCopy({{0.5, 0.5, 0.5}}).basisVectorLengths();
            std::transform(std::begin(xInternal), std::end(xInternal), std::begin(gridCenter), std::begin(xInternal), std::minus<real>());
            RVec rTranslated {xInternal.data()};
            orientation_.rotate_backwards(rTranslated);
            std::transform(std::begin(rTranslated.as_vec()), std::end(rTranslated.as_vec()), std::begin(gridCenter), std::begin(xInternal), std::plus<real>());

            return xInternal;
        };
        /*! \brief
         * \param[in] x Vector to be translated from grid into external coordinates
         * \returns the translated vector.
         */
        NdVector translateFromGrid_(const NdVector &x) const
        {
            RVec xExternal;

            // rotate around grid center: shift to center, rotate, shift back
            auto gridCenter = this->cell().scaledCopy({{0.5, 0.5, 0.5}}).basisVectorLengths();
            std::transform(std::begin(x), std::end(x), std::begin(gridCenter), std::begin(xExternal.as_vec()), std::minus<real>());
            orientation_.rotate(xExternal);
            std::transform(std::begin(xExternal.as_vec()), std::end(xExternal.as_vec()), std::begin(gridCenter), std::begin(xExternal.as_vec()), std::plus<real>());

            // add grid translation to shifted vector
            std::transform(std::begin(xExternal.as_vec()), std::end(xExternal.as_vec()), std::begin(translation_), std::begin(xExternal.as_vec()), std::plus<real>());

            return {{xExternal[0], xExternal[1], xExternal[2]}};

        };
        //! The external coordinates of lattice point (0,...,0)
        NdVector                     translation_;
        Quaternion                   orientation_;
};

}

#endif /* end of include guard: GMX_MATH_ROTATEDGRID_H */
