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
 * Helper methods to place particle COM in boxes.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \inlibraryapi
 * \ingroup module_pbcutil
 */
#ifndef GMX_PBCUTIL_COM_H
#define GMX_PBCUTIL_COM_H

#include <algorithm>
#include <memory>

#include "gromacs/math/vec.h"
#include "gromacs/utility/arrayref.h"

#include "pbcenums.h"

struct gmx_mtop_t;
enum class PbcType : int;
namespace gmx
{

//! How COM shifting should be applied.
enum class COMShiftType : int
{
    Residue,
    Molecule,
    Count
};

/*! \brief
 * Shift all coordinates.
 *
 * Shift coordinates by a previously calculated value.
 *
 * Can be used to e.g. place particles in a box.
 *
 * \param[in] shift Translation that should be applied.
 * \param[in] x Coordinates to translate.
 */
void shiftAtoms(const RVec& shift, ArrayRef<RVec> x);

/*! \brief
 * Moves collection of atoms along the center of mass into a box.
 *
 * This ensures that the centre of mass (COM) of a molecule is placed
 * within a predefined coordinate space (usually a simulation box).
 *
 * \param[in]      pbcType      What kind of PBC are we handling today.
 * \param[in]      unitCellType Kind of unitcell used for the box.
 * \param[in]      centerType   How atoms should be centered.
 * \param[in]      box          The currently available box to place things into.
 * \param[in, out] x            View in coordinates to shift.
 * \param[in]      mtop         Topology with residue and molecule information.
 * \param[in]      comShiftType Whether residues or molecules are shifted.
 */
void placeCoordinatesWithCOMInBox(const PbcType&    pbcType,
                                  UnitCellType      unitCellType,
                                  CenteringType     centerType,
                                  const matrix      box,
                                  ArrayRef<RVec>    x,
                                  const gmx_mtop_t& mtop,
                                  COMShiftType      comShiftType);

} // namespace gmx

#endif
