/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.1
 * Copyright (c) 1991-2001, University of Groningen, The Netherlands
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
 * GROningen Mixture of Alchemy and Childrens' Stories
 */

#ifndef _gen_dum_h
#define _gen_dum_h

static char *SRCID_gen_dum_h = "$Id$";
#include "typedefs.h"
#include "grompp.h"
#include "hackblock.h"

/* stuff for pdb2gmx */

extern void do_dummies(int nrtp, t_restp rtp[], t_atomtype *atype, 
		       t_atoms *at, t_symtab *symtab, rvec *x[], 
		       t_params plist[], int *dummy_type[], int *cgnr[], 
		       real mHmult, bool bDummyAromatics);

extern void do_h_mass(t_params *psb, bool is_dum[], t_atoms *at, real mHmult,
		      bool bDeuterate);

#endif	/* _gen_dum_h */
