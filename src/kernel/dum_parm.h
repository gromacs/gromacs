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
#ifndef _dum_parm_h
#define _dum_parm_h

static char *SRCID_dum_parm_h = "$Id$";

#include "typedefs.h"
#include "grompp.h"

extern void set_dummies(bool bVerbose, t_atoms *atoms,  t_atomtype atype,
			t_params plist[]);
/* set parameters for dummy atoms */

extern void set_dummies_ptype(bool bVerbose, t_idef *idef, t_atoms *atoms);
/* set ptype to Dummy for dummy atoms */

extern void clean_dum_bad(t_params *ps, int natoms);
/* remove all bonds, angles and (im)proper diherals that have become 
   obsolete due to dummy constructions */

#endif	/* _dum_parm_h */
