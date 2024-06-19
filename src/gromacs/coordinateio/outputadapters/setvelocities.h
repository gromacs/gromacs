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
 * Declares gmx:SetVelocities.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \inpublicapi
 * \ingroup module_coordinateio
 */
#ifndef GMX_COORDINATEIO_SETVELOCITIES_H
#define GMX_COORDINATEIO_SETVELOCITIES_H

#include <memory>

#include "gromacs/coordinateio/coordinatefileenums.h"
#include "gromacs/coordinateio/ioutputadapter.h"

struct t_trxframe;

namespace gmx
{

/*!\brief
 * SetVelocities class allows changing writing of velocities to file.
 *
 * This class allows the user to define if velocities should be written
 * to the output coordinate file, and checks if they are available from the
 * currently processed data.
 *
 * \inpublicapi
 * \ingroup module_coordinateio
 *
 */
class SetVelocities : public IOutputAdapter
{
public:
    /*! \brief
     * Construct SetVelocities object with choice for boolean value.
     *
     * Can be used to initialize SetVelocities from outside of trajectoryanalysis
     * with the user specified option to write coordinate velocities or not.
     */
    explicit SetVelocities(ChangeSettingType velocity) : velocity_(velocity)
    {
        if (velocity_ == ChangeSettingType::Never)
        {
            moduleRequirements_ = CoordinateFileFlags::Base;
        }
        else
        {
            moduleRequirements_ = CoordinateFileFlags::RequireVelocityOutput;
        }
    }
    /*! \brief
     *  Move constructor for SetVelocities.
     */
    SetVelocities(SetVelocities&& old) noexcept = default;

    ~SetVelocities() override {}

    /*! \brief
     * Change coordinate frame information for output.
     *
     * In this case, the correct flag for writing the velocities is applied
     * to the output frame, depending on user selection and availability
     * in the input data.
     *
     * \param[in] input Coordinate frame to be modified later.
     */
    void processFrame(int /*framenumber*/, t_trxframe* input) override;

    void checkAbilityDependencies(unsigned long abilities) const override;

private:
    //! Flag to specify if velocities should be written.
    ChangeSettingType velocity_;
    //! Local requirements determined from user input.
    CoordinateFileFlags moduleRequirements_;
};

//! Smart pointer to manage the object.
using SetVelocitiesPointer = std::unique_ptr<SetVelocities>;

} // namespace gmx

#endif
