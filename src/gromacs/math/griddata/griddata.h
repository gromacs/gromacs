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
/*! \defgroup module_griddata Data On N-dimensional Regular Grids
 * \ingroup group_utilitymodules
 * \brief
 * Provides functionality for data manipulation on N-dimensional grids.
 *
 * To represent data on N-dimensional grids, a linear data vector is adressed
 * with a multi-dimensional index on an N-dimensional lattice (\ref gmx::ColumnMajorLattice).
 * To relate this lattice to coordinates, a vector basis is added to the lattice
 * that translates external coordinates into an internal grid representation (\ref gmx::CanonincalVectorBasis).
 * The \ref gmx::IGrid combines lattice indexing and spacial representation with a vector basis.
 * Grids have the responsibility to place N-dimensional multi-indices in space
 * and may be equipped with additional properties like translation (\ref gmx::GridWithTranslation)
 * and periodicity (not implemented yet).
 *
 * \author Christian Blau <cblau@gwdg.de>
 */
/*! \file
 * \brief
 * Declates GridData.
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \inpublicapi
 * \ingroup module_griddata
 */

#ifndef GMX_MATH_GRIDDATA_H_
#define GMX_MATH_GRIDDATA_H_

#include <vector>
#include "grid.h"
#include "gromacs/math/gmxcomplex.h"
#include "gromacs/utility/real.h"

namespace gmx
{

/*! \brief
 * Data on an N-dimensional grid.
 *
 * Ensures that it can hold as many elements as lattice points in the grid.
 * \tparam TContainer Data type container for the data on the grid
 * \tparam N Number of grid dimensions */
template <class TContainer, int N>
class GridData : public TContainer
{
    public:
        /*! \brief
         * Default constructor for convenience.
         *
         * \TODO: causes segfaults when grid functions are accessed before initialisation  */
        GridData() = default;
        /*! \brief
         * Consumes a grid to construct GridData.
         *
         * Allocates data container elements matching to the number of lattice points in the grid.
         * \param[in] grid Pointer to interface that describes the grid. */
        GridData(std::unique_ptr < IGrid < N>> grid) : grid_ {std::move(grid)}
        {
            this->resize(grid_->lattice().getNumLatticePoints());
        }
        /*! \brief
         * Copies a grid to construct GridData.
         *
         * Allocates data container elements matching to the number of lattice points in the grid.
         * \param[in] grid Pointer to interface that describes the grid. */
        GridData(const IGrid < N> &grid) : grid_ {std::move(grid.duplicate())}
        {
            this->resize(grid_->lattice().getNumLatticePoints());
        }

        /*! \brief
         * Copy constructor for GridData.
         * \param[in] other GridData to be copied from. */
        GridData(const GridData &other) : TContainer(other), grid_ {other.grid_->duplicate()}
        {
            this->resize(grid_->lattice().getNumLatticePoints());
        }

        /*! \brief
         * Copy assignment for GridData.
         * \param[in] other GridData to be copied from. */
        GridData &operator=(const GridData &other)
        {
            this->TContainer::operator=(other);
            grid_ = other.grid_->duplicate();
            this->resize(grid_->lattice().getNumLatticePoints());
            return *this;
        }
        /*! \brief
         * Access reference to the underlying grid.
         * \returns immutable handle to the grid interface */
        const IGrid<N> &getGrid() const
        {
            return *grid_;
        }
        /*! \brief
         * Set a new grid and resize data container accordingly. */
        void setGrid(std::unique_ptr < IGrid < N>> &&grid)
        {
            grid_ = std::move(grid);
            this->resize(grid_->lattice().getNumLatticePoints());
        }
        /*! \brief
         * Yield an iterator to a grid element via multi index.
         * \throws std::out_of_range if index is out of bounds for grid.
         * \returns Iterator to data element in grid at given index. */
        typename TContainer::iterator iteratorAtMultiIndex(const typename IGrid<N>::MultiIndex &index)
        {
            return std::begin(*this) + grid_->lattice().lineariseVectorIndex(index);
        }
        /*! \brief
         * Yield an iterator to a grid element via multi index.
         * \throws std::out_of_range if index is out of bounds for grid.
         * \returns Iterator to data element in grid at given index. */
        typename TContainer::const_iterator iteratorAtMultiIndex(const typename IGrid<N>::MultiIndex &index) const
        {
            return std::begin(*this) + grid_->lattice().lineariseVectorIndex(index);
        }

    private:
        //! Reference to grid.
        std::unique_ptr < IGrid < N>> grid_;
};

//! Three-dimensional real numbers on a grid.
typedef GridData<std::vector<float>, DIM> GridDataReal3D;
//! Three-dimensional complex numbers on a grid.
typedef GridData<std::vector<t_complex>, DIM> GridDataComplex3D;

}      // namespace gmx

#endif // GMX_MATH_GRIDDATA_H_
