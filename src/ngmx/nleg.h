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

#ifndef _nleg_h
#define _nleg_h

static char *SRCID_nleg_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) nleg.h 1.19 9/30/97"
#endif /* HAVE_IDENT */
#include "x11.h"
#include "xutil.h"
#include "writeps.h"

typedef struct {
  t_windata wd;
} t_legendwin;

extern unsigned long Type2Color(char *type);
/* Return the color for a given atomtype */

extern t_rgb *Type2RGB(char *type);
/* Return the color for a given atomtype */

extern t_legendwin *init_legw(t_x11 *x11,Window Parent,
			    int x,int y,int width,int height,
			    unsigned long fg,unsigned long bg);

extern void map_legw(t_x11 *x11,t_legendwin *lw);

extern void done_legw(t_x11 *x11,t_legendwin *lw);

#endif	/* _nleg_h */
