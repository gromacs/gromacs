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

#ifndef	_topdirs_h
#define	_topdirs_h

#ifdef HAVE_IDENT
#ident	"@(#) topdirs.h 1.30 9/30/97"
#endif /* HAVE_IDENT */

#include "grompp.h"

typedef struct tagDirStack {
  directive d;
  struct tagDirStack *prev;
} DirStack;

extern int ifunc_index(directive d,int type);

extern char *dir2str (directive d);

extern directive str2dir (char *dstr);

extern void DS_Init (DirStack **DS);

extern void DS_Done (DirStack **DS);

extern void DS_Push (DirStack **DS, directive d);

extern int  DS_Search (DirStack *DS, directive d);

extern int  DS_Check_Order (DirStack *DS, directive d);

#endif	/* _topdirs_h */
