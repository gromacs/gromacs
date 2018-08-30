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

namespace gmx
{

/*! \internal
 *  \brief
 *  Get the internal file type from the \p filename.
 *
 *  \param[in] filename Filename of output file.
 *  \throws InvalidInputError When unable to work on an emoty file name.
 *  \returns integer value of file type.
 */
static int getFileType(std::string filename)
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

/*! \internal
 *  \brief
 * Get the flag representing the requirements for a given file output.
 *
 * Also checks of the supplied topology is sufficient through the pointer
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
    flag |= ICoordinateOutput::efPositionOutput;
    switch (filetype)
    {
        case (efTNG):
            if (mtop == nullptr)
            {
                GMX_THROW(InconsistentInputError("Need topology to write TNG file"));
            }
            flag |= (ICoordinateOutput::efForceOutput |
                     ICoordinateOutput::efVelocityOutput |
                     ICoordinateOutput::efAtomOutput |
                     ICoordinateOutput::efCustomPrecision);
            break;
        case (efPDB):
            if (mtop == nullptr)
            {
                GMX_THROW(InconsistentInputError("Need topology to write PDB file"));
            }
            flag |= (ICoordinateOutput::efConnectionOutput |
                     ICoordinateOutput::efAtomOutput);
            break;
        case (efGRO):
            if (mtop == nullptr)
            {
                GMX_THROW(InconsistentInputError("Need topology to write PDB file"));
            }
            flag |= (ICoordinateOutput::efAtomOutput |
                     ICoordinateOutput::efVelocityOutput);
            break;
        case (efTRR):
            flag |= (ICoordinateOutput::efForceOutput |
                     ICoordinateOutput::efVelocityOutput);
            break;
        case (efXTC):
            flag |= (ICoordinateOutput::efCustomPrecision);
            break;
        case (efG96):
            break;
        default:
            GMX_THROW(InvalidInputError("Invalid file type"));
    }
    return flag;
}

OutputManagerPointer
createOutputManager(const gmx_mtop_t        *mtop,
                    const Selection         &sel,
                    const std::string       &filename,
                    CoordinateOutputAdapters adapters)
{
    int                  filetype      = getFileType(filename);
    unsigned long        flag          = getFlagAndValidateSetting(filetype, mtop);
    OutputManagerPointer outputManager =
        compat::make_unique<OutputManager::OutputManagerBuildHelper>(filename,
                                                                     filetype,
                                                                     flag,
                                                                     sel,
                                                                     mtop,
                                                                     std::move(adapters));
    return outputManager;
}

} // namespace gmx
