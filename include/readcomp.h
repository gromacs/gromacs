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

#ifndef _readcomp_h
#define _readcomp_h

static char *SRCID_readcomp_h = "$Id$";

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
} t_procd;

typedef struct
{
  int nproc;
  t_procd *procs;
} t_procdl;

/* structs are going to be used to determine the load of a processor */

extern void read_compiled(char *compiled,t_procdl **dpl,t_atoms *atoms,
                          t_nbs *nb,t_exclrec *excl);

#endif	/* _readcomp_h */
