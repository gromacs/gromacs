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
 * Enum class defining the different requirements that outputadapters
 * have for the output file type. OutputManager objects can only be built
 * with OutputAdapters whose requirements can be implemented with the available input.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \libinternal
 * \ingroup module_coordinateio
 */
#ifndef GMX_COORDINATEIO_ENUMS_H
#define GMX_COORDINATEIO_ENUMS_H

namespace gmx
{

/*!\brief
 * The enums here define the flags specifying the requirements
 * of different outputadapter modules.
 *
 * When building the object for deciding on the output to a new coordinate file,
 * the outputmanager needs to be able to validate that the dependencies of
 * attached IOutputAdapters are fulfilled. Classes and interfaces that use
 * the enum can check their dependencies against the information encoded in the
 * flags and can then perform an appropriate reaction if there is a mismatch.
 *
 * \todo This should be converted to on enum class, but then it won't be able to
 * perform the check of the flags as it is done now.
 *
 * \libinternal
 * \ingroup module_coordinateio
 *
 */
enum
{
    /*! \brief
     * Base setting that says that the module has no requirements.
     *
     * Sets the flags to default setting to make sure all output methods
     * are supported.
     */
    efBaseOutputManager = 1<<0,
    /*! \brief
     * Requires output method to support force output.
     *
     * If set, only output method supporting writing of forces to
     * the coordinate file will work, others will generate an invalid
     * input error.
     */
    efChangeForceModule     = 1<<1,
    /*! \brief
     * Requires output method to support velocity output.
     *
     * If set, only writing to files that support velocity output will succeed.
     * Other writing methods will generate an error.
     *
     */
    efChangeVelocityModule        = 1<<2,
    /*! \brief
     * Requires output to support changes to selection of coordinates.
     *
     * Default for most methods, will need to be able to write coordinates to
     * output file or generate an error.
     */
    efChangeCoordinateSelectionModule      = 1<<3,
    /*! \brief
     * Requires support for connection information in output format.
     *
     * If set, only file output that supports writing of connection information will succeed.
     * This means for now that only PDB and TNG files can be written. Other file writing
     * methods will fail.
     */
    efChangeConnectionModule    = 1<<4,
    /*! \brief
     * Requires that output format supports the writing of atom information to the file.
     *
     * If set, files will only be written if they can output the information from t_atoms
     * and otherwise create an error while writing.
     */
    efChangeAtomInformationModule = 1<<5,
    /*! \brief
     * Requires that output format supports custom output precision.
     *
     * If set, the coordinate files will only be written if they support the writing
     * of custom precision of the included data.
     */
    efChangeOutputPrecisionModule = 1<<6,
    /*! \brief
     * Requires that output format supports changes to time in coordinate frame
     */
    efChangeTimeModule = 1<<7,
    /*! \brief
     * Requires that output format supports changing box information.
     */
    efChangeBoxModule = 1<<8,
    /*! \brief
     * Dummy value for testing.
     */
    efDummyModule = 1<<15,
};

//! Enum class for setting basic flags in a t_trxframe
enum class ChangeSettingType
{
    efUnchanged,
    efUserYes,
    efUserNo
};
//! Mapping for enums from ChangeSettingType.
const char *const cChangeSettingTypeEnum[] = {
    "unchanged", "yes", "no"
};

//! Enum class for t_atoms settings
enum class ChangeAtomsType
{
    efUnchanged,
    efUserYes,
    efUserNo,
    efRequired
};

//! Mapping for enums from ChangeAtomsType.
const char *const cChangeAtomsTypeEnum[] = {
    "unchanged", "yes", "no", "required"
};

//! Enum class for setting fields new or not.
enum class ChangeFrameInfoType
{
    efUnchanged,
    efUserYes
};
//! Mapping for enums from ChangeFrameUnchangedYesType.
const char *const cChangeFrameInfoTypeEnum[] = {
    "unchanged", "yes"
};

//! Enum class for setting frame time from user input.
enum class ChangeFrameTimeType
{
    efUnchanged,
    efStartTime,
    efTimeStep,
    efBothTime
};

//! Mapping for values from changing frame time.
const char *const cChangeFrameTimeTypeEnum[] = {
    "unchanged", "starttime", "timestep", "both"
};


} // namespace gmx

#endif
