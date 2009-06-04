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

#ifndef _sorting_h
#define _sorting_h

#include "typedefs.h"
typedef atom_id t_bond[2];

extern void sort_bonds(t_topology *top);
/*
 * Sort_bonds sorts all bonded force parameters in order of ascending
 * atom id of the maximum atom id specified in a bond per type bond.
 *
 * If, for example, for a specific bond type the following bonds are specified:
 *    bond1 between atoms 15, 16, 12, 18 and 20
 *    bond2 between atoms 14, 13, 12, 18 and 19
 *    bond3 between atoms 17, 11, 19, 21 and 15
 * then the maximum atom id for each bond would be:
 *    bond1: 20
 *    bond2: 19
 *    bond3: 21
 * so order in which the bonds will be sorted is bond2, bond1, bond3
 *
 * This routine is used to determine to which node a bonds should be
 * allocated. For the distribution of bonds it is necessary to keep all
 * the needed atoms when calculating a bonded force on one node. In
 * this way we prevent communication overhead in bonded force calculations.
 *
 */

extern void sort_xblock(t_block *block,rvec x[],int renum[]);
/*
 * Sort_xblock returns a renumber table which can be used to sort the 
 * blocks specified in block in an order dependent of the coordinates.
 */

extern void sort_bond_list(t_bond bonds[],int nr);
/*
 * Sort_bond_list sort the list of bonds, specified by bonds in order
 * of ascending atom id. The bonds are specified as pairs of atom ids.
 * Where nr specifies the number of bonds (the length of the array).
 */
                  
#endif	/* _sorting_h */
