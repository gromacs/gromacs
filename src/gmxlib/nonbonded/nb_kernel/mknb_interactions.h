/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.3.3
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
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
 * Groningen Machine for Chemical Simulation
 */
#ifndef _MKNB_INTERACTIONS_H_
#define _MKNB_INTERACTIONS_H_


/*! \file  mknb_interactions.h
 *  \brief Kernel generator (only for compile): Parwise interactions
 *
 *  \internal
 *
 *  \note This file is only used to generate the inner loop kernels
 *        at compile time, which in turn are included in Gromacs. 
 *        This code itself is NOT linked into the Gromacs library, so
 *        it does not need to be threadsafe.
 *
 *  This file is only used to generate the inner loop kernels
 *  at compile time, which in turn are included in Gromacs. 
 *  This code itself is NOT linked into the Gromacs library, so
 *  it does not need to be threadsafe.
 */

 

/*! \brief Kernel generator (only for compile): Evaluate an interaction 
 *
 *  \internal
 *
 *  \note   Only defined/used in the nonbonded kernel generator 
 *          program mknb. This program is run once at compile
 *          time to create the inner loops, and then discarded.
 *          This source is NOT linked into any Gromacs library.
 *
 * This is the routine that writes the actual interaction between
 * a pair of atoms in mknb_interaction.c. 
 * 
 * The other routines of the nonbonded kernel generator are merely
 * responsible for iterating over atoms, loading coordinates, storing
 * forces, etc - but this is where it happens.
 *
 *  When implementing a new type of interaction, you should change
 * this routine to add your interaction, and then (usually) modify
 * the loading of parameters in the outer and inner loop files.
 *
 * \param   rsq    Square distance (r*r) for the particle pair
 * \param   rinv   Inverse distance (1/r) for the particle pair
 *
 * \return  Number of floating-point operations used for this
 *          part of the nonbonded kernel. This is only used for
 *          the flopcount reporting in Gromacs, the simulation
 *          results do not depend on it being correct.
 */
int
mknb_calculate_interaction     (char *       rsq, 
								char *       rinv);



#endif /* _MKNB_INTERACTIONS_H_ */

