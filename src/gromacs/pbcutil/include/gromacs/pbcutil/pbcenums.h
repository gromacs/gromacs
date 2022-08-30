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
/*! \file
 * \brief
 * Defines enum classes for centering and unit cell types.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \libinternal
 * \ingroup module_pbcutil
 */
#ifndef GMX_PBCUTIL_PBCENUMS_H
#define GMX_PBCUTIL_PBCENUMS_H

namespace gmx
{

/*! \brief
 * Helper enum class to define centering types.
 */
enum class CenteringType : int
{
    Triclinic,
    Rectangular,
    Zero,
    Count
};

/*! \brief
 * Helper enum class to define Unit cell representation types.
 */
enum class UnitCellType : int
{
    Triclinic,
    Rectangular,
    Compact,
    Count
};

/*! \brief
 * Get names for the different centering types.
 *
 * \param[in] type What name needs to be provided.
 */
const char* centerTypeNames(CenteringType type);

/*! \brief
 * Get names for the different unit cell representation types.
 *
 * \param[in] type What name needs to be provided.
 */
const char* unitCellTypeNames(UnitCellType type);

} // namespace gmx

#endif
