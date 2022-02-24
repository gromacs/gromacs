/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
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

#ifndef GMX_GMXPREPROCESS_GPP_NEXTNB_H
#define GMX_GMXPREPROCESS_GPP_NEXTNB_H

struct InteractionsOfType;

namespace gmx
{
template<typename>
class ArrayRef;
template<typename>
class ListOfLists;
} // namespace gmx

struct t_nextnb
{
    int nr;   /* nr atoms (0 <= i < nr) (atoms->nr)	        */
    int nrex; /* with nrex lists of neighbours		*/
    /* respectively containing zeroth, first	*/
    /* second etc. neigbours (0 <= nre < nrex)	*/
    int** nrexcl; /* with (0 <= nrx < nrexcl[i][nre]) neigbours    */
    /* per list stored in one 2d array of lists	*/
    int*** a; /* like this: a[i][nre][nrx]			*/
};

void init_nnb(t_nextnb* nnb, int nr, int nrex);
/* Initiate the arrays for nnb (see above) */

void done_nnb(t_nextnb* nnb);
/* Cleanup the nnb struct */

#ifdef DEBUG_NNB
#    define print_nnb(nnb, s) __print_nnb(nnb, s)
void print_nnb(t_nextnb* nnb, char* s);
/* Print the nnb struct */
#else
#    define print_nnb(nnb, s)
#endif

void gen_nnb(t_nextnb* nnb, gmx::ArrayRef<InteractionsOfType> plist);
/* Generate a t_nextnb structure from bond information.
 * With the structure you can either generate exclusions
 * or generate angles and dihedrals. The structure must be
 * initiated using init_nnb.
 */

void generate_excl(int nrexcl, int nratoms, gmx::ArrayRef<InteractionsOfType> plist, gmx::ListOfLists<int>* excls);
/* Generate an exclusion block from bonds and constraints in
 * plist.
 */

#endif
