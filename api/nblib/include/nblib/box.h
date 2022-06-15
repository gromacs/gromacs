/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
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
/*! \inpublicapi \file
 * \brief
 * Implements nblib simulation box
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */
#ifndef NBLIB_BOX_H
#define NBLIB_BOX_H

#include "nblib/basicdefinitions.h"

namespace nblib
{

/*! \brief Box specifies the size of periodic simulation systems
 *
 * \inpublicapi
 * \ingroup nblib
 *
 * Currently only cubic and rectangular boxes are supported.
 *
 */
class Box final
{
public:
    using LegacyMatrix = matrix;

    //! Construct a cubic box.
    explicit Box(real l);

    //! Construct a rectangular box.
    Box(real x, real y, real z);

    //! Return the full matrix that specifies the box. Used for gromacs setup code.
    [[nodiscard]] LegacyMatrix const& legacyMatrix() const { return legacyMatrix_; }

private:
    //! \brief check two boxes for equality
    friend bool operator==(const Box& rhs, const Box& lhs);

    //! Stores data in the GROMACS legacy data type
    LegacyMatrix legacyMatrix_;
};

} // namespace nblib
#endif // NBLIB_BOX_H
