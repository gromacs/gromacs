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

#ifndef _pulldown_h
#define _pulldown_h

static char *SRCID_pulldown_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) pulldown.h 1.4 11/23/92"
#endif /* HAVE_IDENT */
#include "popup.h"

typedef struct {
  t_windata wd;
  int       nmenu;
  int       nsel;
  int       *xpos;
  t_menu    **m;
  char      **title;
} t_pulldown;

extern t_pulldown *init_pd(t_x11 *x11,Window Parent,int width,int height,
			   unsigned long fg,unsigned long bg,
			   int nmenu,int *nsub,t_mentry *ent[],char **title);
/* nmenu is the number of submenus, title are the titles of
 * the submenus, nsub are the numbers of entries in each submenu
 * ent are the entries in the pulldown menu, analogous to these in the
 * popup menu.
 * The Parent is the parent window, the width is the width of the parent
 * window. The Menu is constructed as a bar on the topside of the window
 * (as usual). It calculates it own height by the font size.
 * !!!
 * !!! Do not destroy the ent structure, or the titles, while using
 * !!! the menu.
 * !!!
 * When the menu is selected, a ClientMessage will be sent to the Parent
 * specifying the selected item in xclient.data.l[0].
 */

extern void hide_pd(t_x11 *x11,t_pulldown *pd);
/* Hides any menu that is still on the screen when it shouldn't */

extern void check_pd_item(t_pulldown *pd,int nreturn,bool bStatus);
/* Set the bChecked field in the pd item with return code
 * nreturn to bStatus. This function must always be called when
 * the bChecked flag has to changed.
 */

extern void done_pd(t_x11 *x11,t_pulldown *pd);
/* This routine destroys the menu pd, and unregisters it with x11 */

extern int pd_width(t_pulldown *pd);
/* Return the width of the window */

extern int pd_height(t_pulldown *pd);
/* Return the height of the window */

#endif	/* _pulldown_h */
