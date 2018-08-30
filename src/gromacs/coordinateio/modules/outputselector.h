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
 * Declares gmx::OutputSelector.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \inpublicapi
 * \ingroup module_coordinateio
 */
#ifndef GMX_FILEIO_OUTPUTSELECTOR_H
#define GMX_FILEIO_OUTPUTSELECTOR_H

#include <algorithm>

#include "gromacs/coordinateio/coordinateoutput.h"
#include "gromacs/selection/selectionoption.h"

namespace gmx
{

/*!\brief
 * OutputSelector class controls setting which coordinates are actual written.
 *
 * This method controls which coordinates are written from the original file
 * according to which flags are set in the coordinate frame, and allocates the new
 * data for coordinates, velocities and forces. The coordinates written are selected
 * through the select mechanic from user input.
 *
 * \inpublicapi
 * \ingroup module_coordinateio
 *
 */
class OutputSelector : public ICoordinateOutput
{
    public:
        /*! \brief
         * Construct OutputSelector object with initial selection.
         *
         * Can be used to initialize OutputSelector from outside of trajectoryanalysis
         * framework.
         */
        explicit OutputSelector(Selection *sel) : ICoordinateOutput(efCoordinateOutput),
                                                  sel_(sel)
        {
        }
        /*! \brief
         * Move constructor for OutputSelector.
         */
        OutputSelector &operator=(OutputSelector &&old) noexcept
        {
            sel_         = old.sel_;
            return *this;
        }
        /*! \brief
         *  Move assignment constructor for OutputSelector.
         */
        OutputSelector(OutputSelector &&old) noexcept : ICoordinateOutput(old.moduleFlags_),
                                                        sel_(old.sel_)
        {
        }

        ~OutputSelector() {}
        /*! \brief
         * Change coordinate frame information for output.
         *
         * Takes the previously internally stored coordinates and saves them again.
         * Applies correct number of atoms in this case.
         *
         * \param[in] input Coordinate frame to be modified later.
         */
        virtual void processFrame(int /*framenumber*/, t_trxframe *input);

    private:

        /*! \brief
         * Selection of atoms that will be written to disk.
         *
         * Internal selection of atoms chosen by the user that will be written
         * to disk during processing.
         */
        Selection                            *sel_;
};

//! Smart pointer to manage the outputselector object.
typedef std::unique_ptr<OutputSelector>
    OutputSelectorPointer;

} // namespace gmx

#endif
