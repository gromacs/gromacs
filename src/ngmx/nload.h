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

#ifndef	_nload_h
#define	_nload_h

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
			    ulong fg,ulong bg);

extern void map_lw(t_x11 *x11,t_loadwin *lw);

extern void set_load(t_x11 *x11,t_loadwin *lw,int nprocs,int load[]);

extern void done_lw(t_x11 *x11,t_loadwin *lw);

#endif	/* _nload_h */
