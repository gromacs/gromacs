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
 * GROwing Monsters And Cloning Shrimps
 */

#ifndef _typedefs_h
#define _typedefs_h

static char *SRCID_typedefs_h = "$Id$";

#ifdef CPLUSPLUS
extern "C" {
#endif

#define STRLEN 4096
#define NOTSET -12345

#include <sys/types.h>
#include "sysstuff.h"
#include "types/simple.h"
#include "types/enums.h"
#include "types/block.h"
#include "types/symtab.h"
#include "types/idef.h"
#include "types/atoms.h"
#include "types/topology.h"
#include "types/energy.h"
#include "types/inputrec.h"
#include "types/nsborder.h"
#include "types/ishift.h"
#include "types/graph.h"
#include "types/nrnb.h"
#include "types/nblist.h"
#include "types/nsgrid.h"
#include "types/commrec.h"
#include "types/forcerec.h"
#include "types/mdatom.h"
#include "types/ifunc.h"
#include "types/filenm.h"
#include "types/group.h"
#include "types/drblock.h"
#include "types/parm.h"
#include "types/matrix.h"
#include "types/edsams.h"

/* Functions to initiate and delete structures */
extern void init_block(t_block *block);
extern void init_atom (t_atoms *at);
extern void init_top (t_topology *top);
extern void init_inputrec(t_inputrec *ir);
extern void done_block(t_block *block);
extern void done_atom (t_atoms *at);
extern void done_symtab(t_symtab *symtab);
extern void done_top(t_topology *top);
extern void done_inputrec(t_inputrec *ir);
/* These functions defined in gmxlib/typedefs.c */

#ifdef CPLUSPLUS
}
#endif


#endif	/* _typedefs_h */
