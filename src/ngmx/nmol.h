/*
 *       @(#) copyrgt.c 1.12 9/30/97
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0b
 * 
 * Copyright (c) 1990-1997,
 * BIOSON Research Institute, Dept. of Biophysical Chemistry,
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
 * GRoups of Organic Molecules in ACtion for Science
 */

#ifndef	_nmol_h
#define	_nmol_h

#ifdef HAVE_IDENT
#ident	"@(#) nmol.h 1.1 11/19/92"
#endif /* HAVE_IDENT */

#include "x11.h"
#include "xutil.h"

extern t_molwin *init_mw(t_x11 *x11,Window Parent,
			 int x,int y,int width,int height,
			 ulong fg,ulong bg);
/* Create the molecule window using the x,y etc. */

extern void map_mw(t_x11 *x11,t_molwin *mw);

extern void z_fill(t_manager *man, real *zz);
extern void create_visibility(t_manager *man);
extern int compare_obj(const void *a,const void *b);
extern int filter_vis(t_manager *man);
extern void set_sizes(t_manager *man,real sx,real sy);

extern bool toggle_hydrogen(t_x11 *x11,t_molwin *mw);
/* Toggle the state of the hydrogen drawing,
 * return the current state
 */

extern void set_bond_type(t_x11 *x11,t_molwin *mw,int bt);
/* Set the state of the atoms drawing. */

extern bool toggle_box (t_x11 *x11,t_molwin *mw);
/* Toggle the state of the box drawing,
 * return the current state
 */

extern void done_mw(t_x11 *x11,t_molwin *mw);

#endif	/* _nmol_h */
