/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2012, The GROMACS development team,
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
 * Gromacs Runs On Most of All Computer Systems
 */

#ifndef _force_flags_h
#define _force_flags_h


#ifdef __cplusplus
extern "C" {
#endif


/* Flags to tell the force calculation routines what (not) to do */

/* The state has changed */
#define GMX_FORCE_STATECHANGED (1<<0)
/* The box might have changed */
#define GMX_FORCE_DYNAMICBOX   (1<<1)
/* Do neighbor searching */
#define GMX_FORCE_NS           (1<<2)
/* Update long-range neighborlists */
#define GMX_FORCE_LRNS         (1<<3)
/* Calculate bonded energies/forces */
#define GMX_FORCE_BONDED       (1<<4)
/* Store long-range forces in a separate array */
#define GMX_FORCE_SEPLRF       (1<<5)
/* Calculate non-bonded energies/forces */
#define GMX_FORCE_NONBONDED    (1<<6)
/* Calculate forces (not only energies) */
#define GMX_FORCE_FORCES       (1<<7)
/* Calculate the virial */
#define GMX_FORCE_VIRIAL       (1<<8)
/* Calculate energies */
#define GMX_FORCE_ENERGY       (1<<9)
/* Calculate dHdl */
#define GMX_FORCE_DHDL         (1<<10)
/* Calculate long-range energies/forces */
#define GMX_FORCE_DO_LR        (1<<11)

/* Normally one want all energy terms and forces */
#define GMX_FORCE_ALLFORCES    (GMX_FORCE_BONDED | GMX_FORCE_NONBONDED | GMX_FORCE_FORCES)


#ifdef __cplusplus
}
#endif

#endif	/* _force_flags_h */
