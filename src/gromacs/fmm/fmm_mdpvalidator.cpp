/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2025- The GROMACS Authors
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

/*! \brief Implements validation of FMM MDP settings via WarningHandler.
 *
 * \author Muhammad Umair Sadiq <mumairsadiq1@gmail.com>
 */

#include "gmxpre.h"

#include "fmm_mdpvalidator.h"

#include "gromacs/fileio/warninp.h"
#include "gromacs/mdrunutility/mdmodulesnotifiers.h"

#include "fmmoptions.h"


gmx::FmmMdpValidator::FmmMdpValidator(const IFmmOptions& fmmOptions) : fmmOptions_(fmmOptions) {}

void gmx::FmmMdpValidator::setWarningHandler(WarningHandler* wi)
{
    wi_ = wi;
}

void gmx::FmmMdpValidator::validateFmmMdpSettings(const MdModulesCoulombTypeInfo& coulombType)
{
    if (wi_ == nullptr)
    {
        GMX_THROW(
                InconsistentInputError("FmmMdpValidator requires a valid WarningHandler to be set "
                                       "before validating MDP parameters."));
    }

    if (coulombType.coulombInteractionType != CoulombInteractionType::Fmm)
    {
        wi_->addError(formatString("FMM requires Coulomb interaction type to be FMM, but got %s",
                                   enumValueToString(coulombType.coulombInteractionType)));
    }

    fmmOptions_.validateMdpOptions(wi_);
}
