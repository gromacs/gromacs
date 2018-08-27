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
 * Declares gmx::ICoordinateOutput interface for modifying coordinate
 * file structures before writing them to disk.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \inlibraryapi
 * \ingroup module_coordinateio
 */
#ifndef GMX_FILEIO_COORDINATEOUTPUT_H
#define GMX_FILEIO_COORDINATEOUTPUT_H

#include <algorithm>
#include <memory>

#include "gromacs/coordinateio/modulecallback.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/classhelpers.h"

namespace gmx
{

/*!\brief
 * CoordinateOutput class for handling trajectory file flag setting and processing.
 *
 * This interface provides the base point upon which modules that modify trajectory frame
 * datastructures should be build. The interface itself does not provide any direct means
 * to modify the data, but only gives the virtual method to perform work on a
 * t_trxframe object.
 *
 * All modules derived from this interface will be also based on the IModuleCallback
 * interface used to communicate if a given module is compatible with an output method
 * implemented elsewhere.
 *
 * \inlibraryapi
 * \ingroup module_coordinateio
 *
 */
class ICoordinateOutput : public IModuleCallback
{
    public:
        /*! \brief
         * Default constructor for ICoordinateOutput interface.
         *
         * When building new modules based on the interface, the module
         * needs to define its own requirements for file writing using the
         * \p moduleFlag parameter passed during construction.
         *
         * \param[in] moduleFlag Bitmask defining the requirements for
         *                       data using the given module.
         */
        ICoordinateOutput(unsigned long moduleFlag) : IModuleCallback(moduleFlag)
        {
        }
        virtual ~ICoordinateOutput()
        {
        }
        //! Move constructor for old object.
        explicit ICoordinateOutput(ICoordinateOutput &&old) : IModuleCallback(old.moduleFlags_)
        {
        }

        /*! \brief
         * Change settings in t_trxframe according to user input.
         *
         * \param[in] framenumber Frame number as reported from the
         *                        trajectoryanalysis framework or set by user.
         * \param[in] input       Pointer to trajectory analysis frame that will
         *                        be worked on.
         */
        virtual void processFrame(const int framenumber, t_trxframe *input) = 0;

        GMX_DISALLOW_COPY_AND_ASSIGN(ICoordinateOutput);

};

//! Smart pointer to manage the output manager object.
typedef std::shared_ptr<ICoordinateOutput>
    CoordinateOutputPointer;

//! Enum class for setting basic flags in a t_trxframe
enum class ChangeSettingType
{
    efUnchanged,
    efUserYes,
    efUserNo,
    efUserPossible
};
//! Mapping for enums from ChangeSettingType.
const char *const cChangeSettingTypeEnum[] = {
    "unchanged", "yes", "no", "possible"
};

//! Enum class for setting fields new or not.
enum class ChangeFrameUnchangedYesType
{
    efUnchanged,
    efUserYes
};
//! Mapping for enums from ChangeFrameUnchangedYesType.
const char *const cChangeFrameUnchangedYesTypeEnum[] = {
    "unchanged", "yes"
};

//! Enum class for setting frame time from user input.
enum class ChangeFrameTimeType
{
    efUnchanged,
    efStartTime,
    efTimeStep
};

//! Mapping for values from changing frame time.
const char *const cChangeFrameTimeTypeEnum[] = {
    "unchanged", "startime", "timestep"
};

} // namespace gmx

#endif
