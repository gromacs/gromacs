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

#ifndef _nload_h
#define _nload_h

static char *SRCID_nload_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) nload.h 1.19 9/30/97"
#endif /* HAVE_IDENT */
#include "x11.h"
#include "xutil.h"

typedef struct {
  t_windata wd;
  int       nprocs;
  int       *load;
} t_loadwin;

extern t_loadwin *init_lw(t_x11 *x11,Window Parent,
			    int x,int y,int width,int height,
			    unsigned long fg,unsigned long bg);

extern void map_lw(t_x11 *x11,t_loadwin *lw);

extern void set_load(t_x11 *x11,t_loadwin *lw,int nprocs,int load[]);

extern void done_lw(t_x11 *x11,t_loadwin *lw);

#endif	/* _nload_h */
