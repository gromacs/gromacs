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
 * Good ROcking Metal Altar for Chronical Sinners
 */

#ifndef _viewit_h
#define _viewit_h

static char *SRCID_viewit_h = "$Id$";

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_IDENT
#ident "$Id$"
#endif /* HAVE_IDENT */

#ifdef CPLUSPLUS
extern "C" {
#endif

#include "typedefs.h"

extern void do_view(char *fn, char *opts);
/* forks off appropriate command to view file.
 * currently eps, xpm, xvg and pdb are supported 
 * defaults are provided, can be overriden with environment vars 
 */

extern void view_all(int nf, t_filenm fnm[]);
/* calls do_view for all viewable output files in fnm[] */
 
#ifdef CPLUSPLUS
}
#endif

#endif	/* _maths_h */
