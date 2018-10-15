/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2018, by the GROMACS development team, led by
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

struct warninp;
struct t_blocka;

namespace gmx
{

// TODO This should be replaced by a vector of vector<int>
struct ExclusionBlocks
{
    int       nr;   /* The number of entries in the list            */
    int       nra2; /* The total number of entries in a			*/
    int      *nra;  /* The number of entries in each a array (dim nr)   */
    int     **a;    /* The atom numbers (dim nr) the length of each element	*/
    /* i is nra[i]						*/
};

void initExclusionBlocks(ExclusionBlocks *b2, int natom);

void doneExclusionBlocks(ExclusionBlocks *b2);

void mergeExclusions(t_blocka *excl, ExclusionBlocks *b2, warninp *wi);

void blockaToExclusionBlocks(t_blocka *b, ExclusionBlocks *b2);

void exclusionBlocksToBlocka(ExclusionBlocks *b2, t_blocka *b);

} // namespace gmx

#endif
