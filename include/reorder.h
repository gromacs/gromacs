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
 * Green Red Orange Magenta Azure Cyan Skyblue
 */

#ifndef _reorder_h
#define _reorder_h

static char *SRCID_reorder_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) reorder.h 1.4 11/23/92"
#endif /* HAVE_IDENT */

extern void reorder(t_topology *topin,t_topology topout[],
                    int nnodes,int load[],int tload[]);
     /*
      * All atoms used in topin are distributed over nnodes topologies in 
      * topout, where all atom id's are reset, start counting at zero. The 
      * bonded force parameters of topin are distributed using the highest 
      * atom id of a bond. The bonds is then placed on the node where 
      * this atom is a home particle. Due to the algorithm, a node 
      * receives the rest of its subsystem from the nodes with a lower 
      * number (modulo nnodes). The array load specifies the number of atom 
      * to be allocated to every node.
      */

#endif	/* _reorder_h */
