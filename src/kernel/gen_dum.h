/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0
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
 * GRoups of Organic Molecules in ACtion for Science
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

extern void do_h_mass(t_params *psb, bool is_dum[], t_atoms *at, real mHmult);

#endif	/* _gen_dum_h */
