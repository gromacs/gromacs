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

#ifndef _logo_h
#define _logo_h

static char *SRCID_logo_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) logo.h 1.15 9/30/97"
#endif /* HAVE_IDENT */
#include <x11.h>
#include <xutil.h>

typedef struct {
  XFontStruct *bigfont;
  XFontStruct *smallfont;
  t_windata   wd;
} t_logo;

extern void show_logo(t_x11 *x11,t_logo *logo);

extern void hide_logo(t_x11 *x11,t_logo *logo);

extern t_logo *init_logo(t_x11 *x11,Window parent);

extern void done_logo(t_x11 *x11,t_logo *logo);

#endif	/* _logo_h */
