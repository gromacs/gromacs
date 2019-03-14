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

#ifndef GMX_TOPOLOGY_EXCLUSIONBLOCKS_H
#define GMX_TOPOLOGY_EXCLUSIONBLOCKS_H

struct t_blocka;

namespace gmx
{

// TODO This should be replaced by a vector of vector<int>
struct ExclusionBlocks
{
    int       nr;   /* The number of entries in the list */
    int       nra2; /* The total number of entries in a */
    int      *nra;  /* The number of entries in each a array (dim nr) */
    int     **a;    /* The atom numbers (dim nr) the length of each element */
    /* i is nra[i] */
};

//! Initialize the content of the pre-allocated blocks.
void initExclusionBlocks(ExclusionBlocks *b2, int natom);

//! Deallocate the content of the blocks.
void doneExclusionBlocks(ExclusionBlocks *b2);

/*! \brief Merge the contents of \c b2 into \c excl.
 *
 * Requires that \c b2 and \c excl describe the same number of
 * particles, if \c b2 describes a non-zero number.
 */
void mergeExclusions(t_blocka *excl, ExclusionBlocks *b2);

//! Convert the exclusions expressed in \c b into ExclusionBlock form
void blockaToExclusionBlocks(t_blocka *b, ExclusionBlocks *b2);

//! Convert the exclusions expressed in \c b into t_blocka form
void exclusionBlocksToBlocka(ExclusionBlocks *b2, t_blocka *b);

} // namespace gmx

#endif
