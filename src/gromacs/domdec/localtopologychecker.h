/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2021, by the GROMACS development team, led by
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
 *
 * \brief This file declares functionality for checking whether
 * local topologies describe all bonded interactions.
 *
 * \inlibraryapi
 * \ingroup module_domdec
 */

#ifndef GMX_DOMDEC_LOCALTOPOLOGYCHECKER_H
#define GMX_DOMDEC_LOCALTOPOLOGYCHECKER_H

#include <optional>

namespace gmx
{

struct LocalTopologyChecker
{
    /*! \brief Data to help check local topology construction
     *
     * Partitioning could incorrectly miss a bonded interaction.
     * However, checking for that requires a global communication
     * stage, which does not otherwise happen during partitioning. So,
     * for performance, we do that alongside the first global energy
     * reduction after a new DD is made. These variables handle
     * whether the check happens, its input for this domain, output
     * across all domains, and the expected value it should match. */
    /*! \{ */
    /*! \brief Number of bonded interactions found in the local
     * topology for this domain. */
    int numBondedInteractionsToReduce = 0;
    /*! \brief Whether to check at the next global communication
     * stage the total number of bonded interactions found.
     *
     * Cleared after that number is found. */
    bool shouldCheckNumberOfBondedInteractions = false;
    /*! \brief The total number of bonded interactions found in
     * the local topology across all domains.
     *
     * Only has a value after reduction across all ranks, which is
     * removed when it is again time to check after a new
     * partition. */
    std::optional<int> numBondedInteractionsOverAllDomains;
    /*! \} */
};

} // namespace gmx

#endif
