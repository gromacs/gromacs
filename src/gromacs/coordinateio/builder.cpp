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
 * \internal
 * \brief
 * Implements gmx::OutputManagerBuilder.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \ingroup module_coordinateio
 */

#include "gmxpre.h"

#include "builder.h"

#include <algorithm>

#include "gromacs/compat/make_unique.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

/*!\brief
 *  Get the internal file type from the \p filename.
 *
 *  \param[in] filename Filename of output file.
 *  \throws InvalidInputError When unable to work on an emoty file name.
 *  \returns integer value of file type.
 */
static int getFileType(const std::string &filename)
{
    int filetype = efNR;
    if (!filename.empty())
    {
        filetype = fn2ftp(filename.c_str());
    }
    else
    {
        GMX_THROW(InvalidInputError("Can not open file with an empty name"));
    }
    return filetype;
}

/*!\brief
 * Get the flag representing the requirements for a given file output.
 *
 * Also checks if the supplied topology is sufficient through the pointer
 * to \p mtop.
 *
 * \param[in] filetype Internal file type used to check requirements.
 * \param[in] mtop Pointer to topology to cross check requirements.
 * \throws InconsistentInputError When file writing requirements and topology don't match.
 * \returns Requirements represent by the bitmask in the return type.
 */
static unsigned long getFlagAndValidateSetting(int filetype, const gmx_mtop_t *mtop)
{
    unsigned long flag = 0;
    flag |= (efBaseOutputManager |
             efChangeBoxModule |
             efChangeCoordinateSelectionModule |
             efChangeTimeModule);
    switch (filetype)
    {
        case (efTNG):
            if (mtop == nullptr)
            {
                GMX_THROW(InconsistentInputError("Need topology to write TNG file"));
            }
            flag |= (efChangeForceModule |
                     efChangeVelocityModule |
                     efChangeConnectionModule |
                     efChangeAtomInformationModule |
                     efChangeOutputPrecisionModule);
            break;
        case (efPDB):
            if (mtop == nullptr)
            {
                GMX_THROW(InconsistentInputError("Need topology to write PDB file"));
            }
            flag |= (efChangeConnectionModule |
                     efChangeAtomInformationModule);
            break;
        case (efGRO):
            if (mtop == nullptr)
            {
                GMX_THROW(InconsistentInputError("Need topology to write PDB file"));
            }
            flag |= (efChangeAtomInformationModule |
                     efChangeVelocityModule);
            break;
        case (efTRR):
            flag |= (efChangeForceModule |
                     efChangeVelocityModule);
            break;
        case (efXTC):
            flag |= (efChangeOutputPrecisionModule);
            break;
        case (efG96):
            break;
        default:
            GMX_THROW(InvalidInputError("Invalid file type"));
    }
    return flag;
}

/*! \brief
 * Checks that new modules are consistent with requested output.
 *
 * Tests the \p newModuleId and \p newModuleRequirements against the previously
 * established abilities of the future OutputManager. If the \p abilities don't match
 * the \p newModuleRequirements, or if \p newModuleId is already found in
 * \p registeredModules, an exception is thrown and no OutputManager is created.
 *
 * \param[in] newModuleId Flag representing the new OutputAdapter to be registered.
 * \param[in] newModuleRequirements Flag indicating the necessary abilities of the outputmanager.
 * \param[in] registeredModules Flag indicating which modules are already present.
 * \param[in] abilities Flag representing what output will be possible.
 * \returns New flag representing registered modules.
 * \throws InternalError When double registering or violating registration order.
 * \throws InconsistentInputError When requirements for output are not met.
 */
static unsigned long canAddNewModule(unsigned long newModuleId,
                                     unsigned long newModuleRequirements,
                                     unsigned long registeredModules,
                                     unsigned long abilities)
{
    if ((registeredModules & newModuleId) == newModuleId)
    {
        GMX_THROW(InternalError("Trying to add module that has already been added"));
    }
    else if ( (newModuleId == efChangeAtomInformationModule) &&
              ((registeredModules & efChangeCoordinateSelectionModule)
               == efChangeCoordinateSelectionModule))
    {
        std::string message = formatString("Trying to add module that changes atom "
                                           "information after changing output selection, this is not supported.");
        GMX_THROW(InternalError(message));
    }
    if ((abilities & newModuleRequirements) == 0u)
    {
        std::string message = formatString("Adding module will violate requirements for "
                                           "output to trajectory file.");
        GMX_THROW(InconsistentInputError(message));
    }
    registeredModules |= newModuleId;
    return registeredModules;
}


OutputManagerPointer
createOutputManager(const gmx_mtop_t        *mtop,
                    const Selection         &sel,
                    const std::string       &filename,
                    OutputAdapters           adapters)
{
    int           filetype          = getFileType(filename);
    unsigned long abilities         = getFlagAndValidateSetting(filetype, mtop);
    unsigned long registeredModules = efBaseOutputManager;

    // Loop over adapter container to check if all requirements are met.
    for (const auto &adapter : adapters)
    {
        registeredModules = canAddNewModule(adapter.module_->getModuleIDFlag(),
                                            adapter.module_->getModuleRequirementFlag(),
                                            registeredModules, abilities);
    }

    OutputManagerPointer outputManager =
        compat::make_unique<OutputManager::OutputManagerBuildHelper>(filename,
                                                                     filetype,
                                                                     sel,
                                                                     mtop,
                                                                     std::move(adapters));
    return outputManager;
}

} // namespace gmx
