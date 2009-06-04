/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*- 
 *
 * 
 * This file is part of Gromacs        Copyright (c) 1991-2004
 * David van der Spoel, Erik Lindahl, University of Groningen.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org
 * 
 * And Hey:
 * Gnomes, ROck Monsters And Chili Sauce
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

