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

#ifndef _gen_ad_h
#define _gen_ad_h

static char *SRCID_gen_ad_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) gen_ad.h 1.17 9/30/97"
#endif /* HAVE_IDENT */
#include "typedefs.h"
#include "toputil.h"
#include "topio.h"
#include "topexcl.h"
#include "resall.h"

extern void gen_pad(t_nextnb *nnb,t_atoms *atoms,bool bH14,t_params plist[],
		    int nrtp,t_restp rtp[],
		    int nsang,t_angres angs[],
		    int nid,t_idihres idihs[],
		    bool bAlldih);

#endif	/* _gen_ad_h */
