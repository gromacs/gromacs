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
 * the CoordinateFile object needs to be able to validate that the dependencies of
 * attached IOutputAdapters are fulfilled. Classes and interfaces that use
 * the enum can check their dependencies against the information encoded in the
 * flags and can then perform an appropriate reaction if there is a mismatch.
 * \todo Use std::bitset<16> for the entries.
 *
 * \libinternal
 * \ingroup module_coordinateio
 *
 */
enum class CoordinateFileFlags : unsigned long
{
    /*! \brief
     * Base setting that says that the module has no requirements.
     *
     * Sets the flags to default setting to make sure all output methods
     * are supported.
     */
    Base = 1 << 0,
    /*! \brief
     * Requires output method to support force output.
     *
     * If set, only output methods supporting writing of forces will work,
     * others will generate an invalid input error.
     */
    RequireForceOutput = 1 << 1,
    /*! \brief
     * Requires output method to support velocity output.
     *
     * If set, only writing to files that support velocity output will succeed.
     * Other writing methods will generate an error.
     */
    RequireVelocityOutput = 1 << 2,
    /*! \brief
     * Requires support for connection information in output format.
     *
     * If set, only file output that supports writing of connection information will succeed.
     * This means for now that only PDB and TNG files can be written. Other file writing
     * methods will fail.
     */
    RequireAtomConnections = 1 << 3,
    /*! \brief
     * Requires that output format supports the writing of atom information to the file.
     *
     * If set, files will only be written if they can output the information from t_atoms
     * and otherwise create an error while writing.
     */
    RequireAtomInformation = 1 << 4,
    /*! \brief
     * Requires that output format supports writing user-specified output precision.
     *
     * If set, output will only be written if the format supports the writing
     * of custom precision of the included data.
     */
    RequireChangedOutputPrecision = 1 << 5,
    /*! \brief
     * Requires that output format supports writing time to the file.
     */
    RequireNewFrameStartTime = 1 << 6,
    /*! \brief
     * Requires that output format supports writing time to the file.
     */
    RequireNewFrameTimeStep = 1 << 7,
    /*! \brief
     * Requires that output format supports writing box information.
     */
    RequireNewBox = 1 << 8,
    /*! \brief
     * Requires output to support changes to selection of coordinates.
     *
     * Default for most methods, will need to be able to write coordinates to
     * output file or generate an error.
     */
    RequireCoordinateSelection = 1 << 9,
    //! Needed for enumeration array.
    Count
};

//! Conversion of flag to its corresponding unsigned long value.
inline unsigned long convertFlag(CoordinateFileFlags flag)
{
    return static_cast<unsigned long>(flag);
}

//! Enum class for setting basic flags in a t_trxframe
enum class ChangeSettingType : int
{
    PreservedIfPresent,
    Always,
    Never,
    Count
};

//! Enum class for t_atoms settings
enum class ChangeAtomsType
{
    PreservedIfPresent,
    AlwaysFromStructure,
    Never,
    Always,
    Count
};

//! Enum class for setting fields new or not.
enum class ChangeFrameInfoType
{
    PreservedIfPresent,
    Always,
    Count
};

//! Enum class for setting frame time from user input.
enum class ChangeFrameTimeType
{
    PreservedIfPresent,
    StartTime,
    TimeStep,
    Both,
    Count
};

} // namespace gmx

#endif
