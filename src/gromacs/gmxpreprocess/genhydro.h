/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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

#ifndef GMX_GMXPREPROCESS_GENHYDRO_H
#define GMX_GMXPREPROCESS_GENHYDRO_H

#include "gromacs/fileio/pdbio.h"
#include "gromacs/gmxpreprocess/hackblock.h"

#ifdef __cplusplus
extern "C" {
#endif

int add_h(t_atoms **pdbaptr, rvec *xptr[],
          int nah, t_hackblock ah[],
          int nterpairs,
          t_hackblock **ntdb, t_hackblock **ctdb,
          int *rN, int *rC, gmx_bool bMissing,
          int **nabptr, t_hack ***abptr,
          gmx_bool bUpdate_pdba, gmx_bool bKeep_old_pdba);
/* Generate hydrogen atoms and N and C terminal patches.
 * int nterpairs is the number of termini pairs in the molecule
 * ntdb[i] and ctdb[i] may be NULL, no replacement will be done then.
 * rN[i] is the residue number of the N-terminus of chain i,
 * rC[i] is the residue number of the C-terminus of chain i
 * if bMissing==TRUE, continue when atoms are not found
 * if nabptr && abptrb, the hack array will be returned in them to be used
 * a second time
 * if bUpdate_pdba, hydrogens are added to *pdbaptr, else it is unchanged
 * return the New total number of atoms
 */

int protonate(t_atoms **atoms, rvec **x, t_protonate *protdata);
/* Protonate molecule according to oplsaa.ff/aminoacids.hdb
 * when called the first time, new atoms are added to atoms,
 * second time only coordinates are generated
 * return the new total number of atoms
 */

void deprotonate(t_atoms *atoms, rvec *x);
/* Deprotonate any molecule: all atoms whose name begins with H will be
 * removed
 */

#ifdef __cplusplus
}
#endif

#endif
