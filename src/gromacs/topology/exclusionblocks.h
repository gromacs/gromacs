/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2019, by the GROMACS development team, led by
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

#ifndef GMX_TOPOLOGY_EXCLUSIONBLOCKS_H
#define GMX_TOPOLOGY_EXCLUSIONBLOCKS_H

#include <vector>

#include "gromacs/utility/arrayref.h"

struct t_blocka;

namespace gmx
{

/*! \brief
 * Describes a single interaction exclusion.
 */
struct ExclusionBlock
{
    //! Atom numbers for exclusion.
    std::vector<int> atomNumber;
    //! Number of atoms in the exclusion.
    int              nra() const { return atomNumber.size(); }
};

/*! \brief Merge the contents of \c b2 into \c excl.
 *
 * Requires that \c b2 and \c excl describe the same number of
 * particles, if \c b2 describes a non-zero number.
 *
 * \returns New ExclusionBlock object.
 */
std::vector<ExclusionBlock> mergeExclusions(t_blocka *excl, gmx::ArrayRef<const ExclusionBlock> b2);

/*! \brief
 * Convert the exclusions expressed in \c b into ExclusionBlock form and
 * append then to \c b2.
 */
std::vector<ExclusionBlock> blockaToExclusionBlocks(const t_blocka *b, gmx::ArrayRef<const ExclusionBlock> b2);

//! Convert the exclusions expressed in \c b into t_blocka form
void exclusionBlocksToBlocka(gmx::ArrayRef<const ExclusionBlock> b2, t_blocka *b);

} // namespace gmx

#endif
