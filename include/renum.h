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
 * Great Red Oystrich Makes All Chemists Sane
 */

#ifndef _renum_h
#define _renum_h

static char *SRCID_renum_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) renum.h 1.6 11/23/92"
#endif /* HAVE_IDENT */

#include "typedefs.h"

extern void renum_params(t_topology *top,int renum[]);
     /*
      * The atom id's in the parameters for the bonded forces will be 
      * renumbered according to the order specified in renum. 
      *
      * renum[i]=j specifies that atom i will be at postion j after 
      * renum_params.
      */

extern void renumber_top(t_topology *top,rvec *x,rvec *v,rvec *f,int renum[]);
     /*
      * All atoms in the topology top will be renumbered according to the order
      * specified in renum. The vector array's x,v & f will also be renumbered.
      *
      * renum[i]=j specifies that atom i will be at postion j after 
      * renumber_top.
      */

#endif	/* _renum_h */
