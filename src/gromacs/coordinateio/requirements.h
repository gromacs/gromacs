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
/*!\file
 * \brief
 * Storage object for requirments to build outputmanager.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \ingroup module_coordinateio
 * \inlibraryapi
 */
#ifndef GMX_COORDINATEIO_OUTPUTREQUIREMENTS_H
#define GMX_COORDINATEIO_OUTPUTREQUIREMENTS_H

#include <vector>

#include "gromacs/coordinateio/enums.h"
#include "gromacs/math/vec.h"

namespace gmx
{

/*!\brief
 * Container for the user input values that will be used by the builder
 * to determine which OutputAdapters should/could/will be registered
 * to the OutputManager.
 */
struct OutputRequirements {
    //! Should velocities be written.
    ChangeSettingType           velocity           = ChangeSettingType::efUnchanged;
    //! Should forces be written.
    ChangeSettingType           force              = ChangeSettingType::efUnchanged;
    //! Should precision be changed.
    ChangeFrameInfoType         precision          = ChangeFrameInfoType::efUnchanged;
    //! Precision used in output file.
    int                         prec               = 3;
    //! If a new precision value has been set.
    bool                        setNewPrecision    = false;
    //! Should frame start time be changed.
    ChangeFrameTimeType         frameTime          = ChangeFrameTimeType::efUnchanged;
    //! Time for first frame to start.
    real                        startTimeValue     = 0;
    //! Time step to use between frames.
    real                        timeStepValue      = 0;
    //! If start time value has been assigned.
    bool                        setNewStartTime    = false;
    //! If a new time step value has been assigned.
    bool                        setNewTimeStep     = false;
    //! User supplied diagonal box vector.
    std::vector<real>           newBoxVector;
    //! If a new box vector has been set.
    bool                        setNewBox          = false;
    //! Box vector converted to matrix format.
    matrix                      newBox             = {{0}};
    //! Should frame box be changed.
    ChangeFrameInfoType         box                = ChangeFrameInfoType::efUnchanged;
    //! Should frame atom setting be changed.
    ChangeAtomsType             atoms              = ChangeAtomsType::efUnchanged;
    //! Have we done our part and checked options before going on?
    bool                        haveCheckedOptions = false;
    //! Add module for testing.
    bool                        addDummyModule = false;
    //! Module ID flag for dummy module.
    unsigned long               dummyIDFlag = 0;
    //! Module requirements flag for dummy module.
    unsigned long               dummyRequirementsFlag = 0;
    //! Check if options have been postprocessed.
    bool                        isValid = false;
};

} // namespace gmx

#endif
