/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 1.6
 * 
 * Copyright (c) 1991-1997
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 *
 * Also check out our WWW page:
 * http://rugmd0.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 *
 * And Hey:
 * GROtesk MACabre and Sinister
 */

#ifndef _reorder_h
#define _reorder_h

static char *SRCID_reorder_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) reorder.h 1.4 11/23/92"
#endif /* HAVE_IDENT */

extern void reorder(t_topology *topin,t_topology topout[],
                    int nprocs,int load[],int tload[]);
     /*
      * All atoms used in topin are distributed over nprocs topologies in 
      * topout, where all atom id's are reset, start counting at zero. The 
      * bonded force parameters of topin are distributed using the highest 
      * atom id of a bond. The bonds is then placed on the processor where 
      * this atom is a home particle. Due to the algorithm, a processor 
      * receives the rest of its subsystem from the processors with a lower 
      * number (modulo nprocs). The array load specifies the number of atom 
      * to be allocated to every processor.
      */

#endif	/* _reorder_h */
