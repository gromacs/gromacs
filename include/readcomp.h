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
 * Green Red Orange Magenta Azure Cyan Skyblue
 */

#ifndef _readcomp_h
#define _readcomp_h

static char *SRCID_readcomp_h = "$Id$";

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_IDENT
#ident	"@(#) readcomp.h 1.12 25 Jan 1993"
#endif /* HAVE_IDENT */

#include "typedefs.h"

/*
 * readcompiled.h
 *
 */

typedef struct
{
  int natom;
  int *atomid;
} t_noded;

typedef struct
{
  int nnode;
  t_noded *nodes;
} t_nodedl;

/* structs are going to be used to determine the load of a node */

extern void read_compiled(char *compiled,t_nodedl **dpl,t_atoms *atoms,
                          t_nbs *nb,t_exclrec *excl);

#endif	/* _readcomp_h */
