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
#include "pdb2gmx.h"

/* stuff for pdb2gmx */

extern void do_dummies(int nrtp, t_restp rtp[], 
		       t_atomtype *atype, real mHmult, 
		       t_atoms *at, t_symtab *symtab, rvec *x[], 
		       t_params plist[], t_params *newbonds,
		       int *dummy_type[], int *cgnr[]);

extern void do_h_mass(t_params *psb, bool is_dum[], t_atoms *at, real mHmult);

extern void clean_dum_bonds(t_params *ps, int dummy_type[]);

extern void clean_dum_angles(t_params *ps, int natom, 
			     t_params *plist, int dummy_type[]);

extern void clean_dum_dihs(t_params *ps, int natom, char dihname[], 
			   t_params *plist, int dummy_type[]);

extern void do_dum_excl(t_block *excl, bool is_dum[]);

#endif	/* _gen_dum_h */
