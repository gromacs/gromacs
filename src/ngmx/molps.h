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
 * Great Red Oystrich Makes All Chemists Sane
 */

#ifndef _molps_h
#define _molps_h

static char *SRCID_molps_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) molps.h 1.10 9/30/97"
#endif /* HAVE_IDENT */
#include "sysstuff.h"
#include "manager.h"
	
extern void ps_draw_mol(FILE *ps,t_manager *man);
/* Draw molecules to a postscript file */

#endif	/* _molps_h */
