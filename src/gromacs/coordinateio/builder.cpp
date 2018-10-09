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
#include "gromacs/coordinateio/outputadaptercontainer.h"
#include "gromacs/coordinateio/outputadapters.h"
#include "gromacs/coordinateio/outputadapters/dummymodule.h"
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
 * \param[in] top      Pointer to topology to cross check requirements.
 * \throws InconsistentInputError When file writing requirements and topology don't match.
 * \returns Requirements represent by the bitmask in the return type.
 */
static unsigned long getSupportedOutputAdapters(int filetype, const gmx_mtop_t *top)
{
    unsigned long supportedOutputAdapters = 0;
    supportedOutputAdapters |= efBaseOutputManager;
    switch (filetype)
    {
        case (efTNG):
            if (top == nullptr)
            {
                GMX_THROW(InconsistentInputError("Need topology to write TNG file"));
            }
            supportedOutputAdapters |= (efChangeForceModule |
                                        efChangeVelocityModule |
                                        efChangeConnectionModule |
                                        efChangeAtomInformationModule |
                                        efChangeOutputPrecisionModule);
            break;
        case (efPDB):
            supportedOutputAdapters |= (efChangeConnectionModule |
                                        efChangeAtomInformationModule);
            break;
        case (efGRO):
            supportedOutputAdapters |= (efChangeAtomInformationModule |
                                        efChangeVelocityModule);
            break;
        case (efTRR):
            supportedOutputAdapters |= (efChangeForceModule |
                                        efChangeVelocityModule);
            break;
        case (efXTC):
            supportedOutputAdapters |= (efChangeOutputPrecisionModule);
            break;
        case (efG96):
            break;
        default:
            GMX_THROW(InvalidInputError("Invalid file type"));
    }
    return supportedOutputAdapters;
}

static OutputAdapterContainer
addOutputAdapters(const OutputRequirements  &requirements,
                  AtomsDataPtr                /* atoms */,
                  const Selection           &sel,
                  unsigned long              abilities)
{
    OutputAdapterContainer output(abilities);

    if (requirements.velocity != ChangeSettingType::efUnchanged)
    {
        // add adapter here
    }
    if (requirements.force != ChangeSettingType::efUnchanged)
    {
        // add adapter here
    }
    if (requirements.precision != ChangeFrameInfoType::efUnchanged)
    {
        // add adapter here
    }
    if (requirements.atoms != ChangeAtomsType::efUnchanged)
    {
        // add adapter here
    }
    if (requirements.frameTime != ChangeFrameTimeType::efUnchanged)
    {
        // add adapter here
    }
    if (requirements.box != ChangeFrameInfoType::efUnchanged)
    {
        // add adapter here
    }
    if (requirements.addDummyModule)
    {
        output.addAdapter(compat::make_unique<DummyOutputModule>(
                                  requirements.dummyRequirementsFlag,
                                  requirements.dummyIDFlag));
    }
    if (sel.isValid())
    {
        // add adapter here
    }
    return output;
}

OutputManagerPointer
createOutputManager(const gmx_mtop_t   *top,
                    const Selection    &sel,
                    const std::string  &filename,
                    AtomsDataPtr        atoms,
                    OutputRequirements  requirements)
{
    int           filetype          = getFileType(filename);
    unsigned long abilities         = getSupportedOutputAdapters(filetype, top);

    if (!requirements.isValid)
    {
        GMX_THROW(InternalError("User options have not been processed"));
    }
    // first, check if we have a special output format that needs atoms
    if ((filetype == efPDB) || (filetype == efGRO))
    {
        if (requirements.atoms == ChangeAtomsType::efUserNo)
        {
            GMX_THROW(InconsistentInputError("Can not write to PDB or GRO when"
                                             "explicitly turning atom information off"));
        }
        if (requirements.atoms != ChangeAtomsType::efUserYes)
        {
            requirements.atoms = ChangeAtomsType::efRequired;
        }
    }
    OutputAdapterContainer outputAdapters = addOutputAdapters(requirements, std::move(atoms),
                                                              sel, abilities);

    OutputManagerPointer   outputManager =
        compat::make_unique<OutputManager::OutputManagerBuildHelper>(filename,
                                                                     filetype,
                                                                     sel,
                                                                     top,
                                                                     std::move(outputAdapters));
    return outputManager;
}

} // namespace gmx
