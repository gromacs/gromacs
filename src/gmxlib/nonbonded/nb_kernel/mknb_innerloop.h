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
#ifndef _MKNB_INNERLOOP_H_
#define _MKNB_INNERLOOP_H_


/*! \file  mknb_innerloop.h
 *  \brief Kernel generator (only for compile): Loop over neighbors in lists
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


/*! \brief Kernel generator (only for compile): Inner loop over neighbors.
 *
 * \internal
 *
 *  \note   Only defined/used in the nonbonded kernel generator 
 *          program mknb. This program is run once at compile
 *          time to create the inner loops, and then discarded.
 *          This source is NOT linked into any Gromacs library.
 *
 * This routine will open a loop, load indices, coordinates and
 * parameters for the inner (also known as 'j') particle.
 * The vectorial distance to the outer ('i') particle is
 * calculated, squared, and used to calculate 1/r through an
 * inverse square root (often heavily optimized).
 *
 * The routine mknb_calculate_interaction() is called to 
 * calculate the scalar force and potential between the two
 * particles, and when it returns we update the vectorial
 * force, subtract it on the 'j' particle and increment the
 * local force variable for the outer 'i' particle.
 *
 * The parameters to be loaded are determined by the current
 * values in the 'func' structure global to this file.
 */
void
mknb_innerloop(void);


#endif /* _MKNB_INNERLOOP_H_ */

