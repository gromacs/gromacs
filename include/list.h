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
 * Great Red Oystrich Makes All Chemists Sane
 */

#ifndef _list_h
#define _list_h

static char *SRCID_list_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) list.h 1.15 2/2/97"
#endif /* HAVE_IDENT */

#include <stdlib.h>
#include "typedefs.h"

segmptr insert_segm();

nodeptr insert_node();

void make_list ();

int count_segms ();

void read_list (FILE *statusfile, char *fn, bool bStat);

#endif	/* _list_h */
