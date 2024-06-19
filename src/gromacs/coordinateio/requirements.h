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
/*!\file
 * \brief
 * Storage object for requirements to build coordinate file writer.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \ingroup module_coordinateio
 * \inlibraryapi
 */
#ifndef GMX_COORDINATEIO_OUTPUTREQUIREMENTS_H
#define GMX_COORDINATEIO_OUTPUTREQUIREMENTS_H

#include <vector>

#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/utility/real.h"

#include "coordinatefileenums.h"

namespace gmx
{

/*!\brief
 * Container for the user input values that will be used by the builder
 * to determine which OutputAdapters should/could/will be registered
 * to the coordinate file writer.
 */
class OutputRequirementOptionDirector
{
public:
    /*! \brief
     * Populate requirements from option interface.
     *
     * \param[in] options Pointer to global options framework.
     */
    void initOptions(IOptionsContainer* options);

    /*! \brief
     * Provide processed requirements to create on coordinate file writing method.
     */
    struct OutputRequirements process() const;

private:
    //! Should velocities be written.
    ChangeSettingType velocity_ = ChangeSettingType::PreservedIfPresent;
    //! Should forces be written.
    ChangeSettingType force_ = ChangeSettingType::PreservedIfPresent;
    //! Precision used in output file.
    int prec_ = 3;
    //! If a new precision value has been set.
    bool setNewPrecision_ = false;
    //! Time for first frame to start.
    real startTimeValue_ = 0;
    //! Time step to use between frames.
    real timeStepValue_ = 0;
    //! If start time value has been assigned.
    bool setNewStartTime_ = false;
    //! If a new time step value has been assigned.
    bool setNewTimeStep_ = false;
    //! User supplied diagonal box vector.
    std::vector<real> newBoxVector_;
    //! If a new box vector has been set.
    bool setNewBox_ = false;
    //! Should frame atom setting be changed.
    ChangeAtomsType atoms_ = ChangeAtomsType::PreservedIfPresent;
};

/*! \brief
 * Finalized version of requirements after processing.
 */
struct OutputRequirements
{
    //! Should velocities be written.
    ChangeSettingType velocity = ChangeSettingType::PreservedIfPresent;
    //! Should forces be written.
    ChangeSettingType force = ChangeSettingType::PreservedIfPresent;
    //! Should precision be changed.
    ChangeFrameInfoType precision = ChangeFrameInfoType::PreservedIfPresent;
    //! Precision used in output file.
    int prec = 3;
    //! Should frame start time be changed.
    ChangeFrameTimeType frameTime = ChangeFrameTimeType::PreservedIfPresent;
    //! Time for first frame to start.
    real startTimeValue = 0;
    //! Time step to use between frames.
    real timeStepValue = 0;
    //! Box vector converted to matrix format.
    matrix newBox = { { 0 } };
    //! Should frame box be changed.
    ChangeFrameInfoType box = ChangeFrameInfoType::PreservedIfPresent;
    //! Should frame atom setting be changed.
    ChangeAtomsType atoms = ChangeAtomsType::PreservedIfPresent;
};

} // namespace gmx

#endif
