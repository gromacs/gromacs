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
 * Implements gmx::OutputAdapterContainer.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \ingroup module_coordinateio
 */

#include "gmxpre.h"

#include "outputadaptercontainer.h"

#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

void
OutputAdapterContainer::addAdapter(OutputAdapterPointer adapter)
{
    unsigned long newModuleId           = adapter->getModuleIDFlag();
    unsigned long newModuleRequirements = adapter->getModuleRequirementFlag();
    if ((registered_ & newModuleId) == newModuleId)
    {
        GMX_THROW(InternalError("Trying to add adapter that has already been added"));
    }
    else if ((registered_ & efChangeCoordinateSelectionModule)
             == efChangeCoordinateSelectionModule)
    {
        GMX_THROW(InternalError(formatString("Trying to add adapter after already setting "
                                             "the adapter for reducing output coordinates. "
                                             "This is not supported.")));
    }
    if ((abilities_ & newModuleRequirements) == 0u)
    {
        GMX_THROW(InconsistentInputError(formatString("Adding module will violate requirements for "
                                                      "output to trajectory file.")));
    }

    registered_   |= newModuleId;
    requirements_ |= newModuleRequirements;
    outputAdapters_.push_back(std::move(adapter));
}

} // namespace gmx
