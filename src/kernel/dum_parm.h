/*
 * $Id$
 * 
 *       This source code is part of
 * 
 *        G   R   O   M   A   C   S
 * 
 * GROningen MAchine for Chemical Simulations
 * 
 *               VERSION 2.0
 * 
 * Copyright (c) 1991-1999
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 * 
 * Also check out our WWW page:
 * http://md.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 * 
 * And Hey:
 * GRowing Old MAkes el Chrono Sweat
 */

#ifndef _dum_parm_h
#define _dum_parm_h

static char *SRCID_dum_parm_h = "$Id$";

#include "typedefs.h"
#include "grompp.h"

extern int set_dummies(bool bVerbose, t_atoms *atoms,  t_atomtype atype,
		       t_params plist[]);
/* set parameters for dummy atoms, return number of dummies */

extern void set_dummies_ptype(bool bVerbose, t_idef *idef, t_atoms *atoms);
/* set ptype to Dummy for dummy atoms */

extern void clean_dum_bondeds(t_params *ps, int natoms, bool bRmDumBds);
/* remove all bonded interaction (bonds, angles and diherals) that
   have become obsolete due to dummy atom constructions */

#endif	/* _dum_parm_h */
