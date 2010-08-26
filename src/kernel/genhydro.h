/*
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
 */

#ifndef _genhydro_h
#define _genhydro_h

#include "pdbio.h"
#include "hackblock.h"

extern int add_h(t_atoms **pdbaptr, rvec *xptr[], 
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
 * if bMissing==TRUE, conitue when atoms are not found
 * if nabptr && abptrb, the hack array will be returned in them to be used
 * a second time
 * if bUpdate_pdba, hydrogens are added to *pdbaptr, else it is unchanged
 * return the New total number of atoms 
 */

extern int protonate(t_atoms **atoms, rvec **x, t_protonate *protdata);
/* Protonate molecule according to ffgmx2 
 * when called the first time, new atoms are added to atoms, 
 * second time only coordinates are generated
 * return the New total number of atoms 
 */

extern void deprotonate(t_atoms *atoms,rvec *x);
/* Deprotonate any molecule: all atoms whose name begins with H will be 
 * removed 
 */

#endif

