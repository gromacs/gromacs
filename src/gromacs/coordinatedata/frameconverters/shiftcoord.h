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
/*! \file
 * \brief
 * Method to shift coordinates in coordinate frame manipulation.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \inlibraryapi
 * \ingroup module_coordinatedata
 */
#ifndef GMX_COORDINATEDATA_MODULES_SHIFTCOORD_H
#define GMX_COORDINATEDATA_MODULES_SHIFTCOORD_H

#include <algorithm>

#include "gromacs/coordinatedata/frameconverters/frameconverter.h"
#include "gromacs/math/vec.h"
#include "gromacs/selection/selection.h"

namespace gmx
{

/*!\brief
 * ShiftCoord class controls setting which coordinates are actual written.
 *
 * This method controls which coordinates are written from the original file
 * according to which flags are set in the coordinate frame, and allocates the new
 * data for coordinates, velocities and forces. The coordinates written are selected
 * through the select mechanic from user input.
 *
 * \inlibraryapi
 * \ingroup module_coordinatedata
 *
 */
class ShiftCoord : public IFrameConverter
{
    public:
        /*! \brief
         * Construct ShiftCoord object with initial selection.
         *
         * Can be used to initialize ShiftCoord from outside of trajectoryanalysis
         * framework.
         * TODO Add initializers for the remaining fields.
         */
        explicit ShiftCoord(RVec shift) : shift_(shift)
        {
        }
        /*! \brief
         *  Move constructor for ShiftCoord.
         */
        ShiftCoord(ShiftCoord &&old) noexcept : shift_(old.shift_)
        {
        }

        ~ShiftCoord() {}

        /*! \brief
         * Change coordinate frame information for output.
         *
         * Takes the previously internally stored coordinates and saves them again.
         * Applies correct number of atoms in this case.
         *
         * \param[in] input Coordinate frame to be modified later.
         */
        virtual void convertFrame(const t_trxframe &input);

    private:
        /*! \brief
         * Coordinate vector by that coordinates will be changed.
         *
         * User or program supplied vector of coordinates by that the
         * frame coordinates will be modified in a unified shift.
         */
        RVec                            shift_;
        Selection                       sel_;
};

//! Smart pointer to manage the outputselector object.
typedef std::unique_ptr<ShiftCoord>
    ShiftCoordPointer;

} // namespace gmx

#endif
