/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019,2021, by the GROMACS development team, led by
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

/*! \libinternal \file
 * \brief Defines atom and atom interaction locality enums
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_mdtypes
 */

#ifndef GMX_MDTYPES_LOCALITY_H
#define GMX_MDTYPES_LOCALITY_H

#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/exceptions.h"

namespace gmx
{

/*! \brief Atom locality indicator: local, non-local, all.
 *
 * Used for calls to:
 * gridding, force calculation, x/f buffer operations
 */
enum class AtomLocality : int
{
    Local    = 0, //!< Local atoms
    NonLocal = 1, //!< Non-local atoms
    All      = 2, //!< Both local and non-local atoms
    Count    = 3  //!< The number of atom locality types
};

/*! \brief Get the human-friendly name for atom localities.
 *
 * \param[in] enumValue The enum value to get the name for.
 */
[[maybe_unused]] static const char* enumValueToString(AtomLocality enumValue)
{
    static constexpr gmx::EnumerationArray<AtomLocality, const char*> atomLocalityNames = {
        "Local", "Non-local", "All"
    };
    return atomLocalityNames[enumValue];
}
/*! \brief Interaction locality indicator: local, non-local, all.
 *
 * Used for calls to:
 * pair-search, force calculation, x/f buffer operations
 */
enum class InteractionLocality : int
{
    Local    = 0, //!< Interactions between local atoms only
    NonLocal = 1, //!< Interactions between non-local and (non-)local atoms
    Count    = 2  //!< The number of interaction locality types
};

/*! \brief Get the human-friendly name for interaction localities.
 *
 * \param[in] enumValue The enum value to get the name for.
 */
[[maybe_unused]] static const char* enumValueToString(InteractionLocality enumValue)
{
    static constexpr gmx::EnumerationArray<InteractionLocality, const char*> interactionLocalityNames = {
        "Local", "Non-local"
    };
    return interactionLocalityNames[enumValue];
}

/*! \brief Convert atom locality to interaction locality.
 *
 *  In the current implementation the this is straightforward conversion:
 *  local to local, non-local to non-local.
 *
 *  \param[in] atomLocality Atom locality specifier
 *  \returns                Interaction locality corresponding to the atom locality passed.
 */
static inline InteractionLocality atomToInteractionLocality(const AtomLocality atomLocality)
{

    /* determine interaction locality from atom locality */
    if (atomLocality == AtomLocality::Local)
    {
        return InteractionLocality::Local;
    }
    else if (atomLocality == AtomLocality::NonLocal)
    {
        return InteractionLocality::NonLocal;
    }
    else
    {
        GMX_THROW(gmx::InconsistentInputError(
                "Only Local and NonLocal atom locities can be converted to interaction locality."));
    }
}

} // namespace gmx

#endif // GMX_MDTYPES_LOCALITY_H
