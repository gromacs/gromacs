/*
 * $Id$
 * 
 *       This source code is part of
 * 
 *        G   R   O   M   A   C   S
 * 
 * GROningen MAchine for Chemical Simulations
 * 
 *               VERSION 2.0
 * 
 * Copyright (c) 1991-1999
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 * 
 * Also check out our WWW page:
 * http://md.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 * 
 * And Hey:
 * GRowing Old MAkes el Chrono Sweat
 */

#ifndef _sorting_h
#define _sorting_h

static char *SRCID_sorting_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) sorting.h 1.21 9/30/97"
#endif /* HAVE_IDENT */

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
