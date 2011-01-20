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

#ifndef _gen_vsite_h
#define _gen_vsite_h

#include "typedefs.h"
#include "grompp.h"
#include "gpp_atomtype.h"
#include "hackblock.h"

/* stuff for pdb2gmx */

extern void do_vsites(int nrtp, t_restp rtp[], gpp_atomtype_t atype, 
		      t_atoms *at, t_symtab *symtab, rvec *x[], 
		      t_params plist[], int *dummy_type[], int *cgnr[], 
		      real mHmult, gmx_bool bVSiteAromatics,
		      const char *ffdir);

extern void do_h_mass(t_params *psb, int vsite_type[], t_atoms *at, real mHmult,
		      gmx_bool bDeuterate);

#endif	/* _gen_vsite_h */
