/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2020, by the GROMACS development team, led by
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
 * Implements nblib particle-types interactions
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 */

#include <sstream>

#include "gmxpre.h"

#include "interactions.h"

#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

namespace nblib
{

void ParticleTypesInteractons::add(ParticleType particleType, C6 c6, C12 c12)
{
    if (particleInteractions_.count(particleType.name()) == 0)
    {
        particleInteractions_[particleType.name()] = std::make_tuple(c6, c12);
    } else {
        std::string message = gmx::formatString("Attempting to add nonbonded interaction %s twice", particleType.name().c_str());
        GMX_THROW(gmx::InvalidInputError(message));
    }
}

void ParticleTypesInteractons::add(ParticleType particleType1, ParticleType particleType2, C6 c6, C12 c12)
{
    std::string interactionKey = particleType1.name() + "-" + particleType2.name();
    std::string possibleInteractionKey = particleType2.name() + "-" + particleType1.name();
    if (particleInteractions_.count(interactionKey) == 0)
    {
        std::string message = gmx::formatString("Attempting to add nonbonded interaction between %s %s twice", particleType1.name().c_str(), particleType2.name().c_str());
        GMX_RELEASE_ASSERT(particleInteractions_.count(possibleInteractionKey) == 0, message.c_str());

        particleInteractions_[interactionKey] = std::make_tuple(c6, c12);
        particleInteractions_[possibleInteractionKey] = std::make_tuple(c6, c12);
    }
}


}
