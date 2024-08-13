/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */

/*! \internal \file
 *
 * \brief
 * Declares the GridSetData struct which holds grid data that is shared over all grids
 *
 * Also declares a struct for work data that is shared over grids.
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_nbnxm
 */

#ifndef GMX_NBNXM_GRIDSETDATA_H
#define GMX_NBNXM_GRIDSETDATA_H

#include <vector>

#include "gromacs/gpu_utils/hostallocator.h"

namespace gmx
{

/*! \internal
 * \brief Struct that holds grid data that is shared over all grids
 *
 * To enable a single coordinate and force array, a single cell range
 * is needed which covers all grids.
 */
struct GridSetData
{
    //! The cell indices for all atoms
    HostVector<int> cells;
    //! The atom indices for all atoms stored in cell order
    HostVector<int> atomIndices;
};

/*! \internal
 * \brief Working arrays for constructing a grid
 */
struct GridWork
{
    //! Number of atoms for each grid column
    std::vector<int> numAtomsPerColumn;
    //! Buffer for sorting integers
    std::vector<int> sortBuffer;
};

} // namespace gmx

#endif
