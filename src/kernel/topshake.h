/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 1.6
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
 * S  C  A  M  O  R  G
 */

#ifndef _topshake_h
#define _topshake_h

static char *SRCID_topshake_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) topshake.h 1.19 9/30/97"
#endif /* HAVE_IDENT */

#include "topio.h"

void make_shake (t_params plist[],t_atoms *atoms,t_atomtype *at,int nshake);

#endif	/* _topshake_h */
