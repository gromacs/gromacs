/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.0
 * 
 * Copyright (c) 1991-2001
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
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
 * Do check out http://www.gromacs.org , or mail us at gromacs@gromacs.org .
 * 
 * And Hey:
 * GRoups of Organic Molecules in ACtion for Science
 */

#ifndef _reorder_h
#define _reorder_h

static char *SRCID_reorder_h = "$Id$";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

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
