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
 * Good ROcking Metal Altar for Chronical Sinners
 */

#ifndef _mshift_h
#define _mshift_h

static char *SRCID_mshift_h = "$Id$";

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_IDENT
#ident	"@(#) mshift.h 1.11 2/2/97"
#endif /* HAVE_IDENT */

#include "typedefs.h"

extern t_graph *mk_graph(t_idef *idef,int natoms,bool bShakeOnly,bool bSettle);
/* Build a graph from an idef description. The graph can be used
 * to generate mol-shift indices.
 * If bShakeOnly, only the connections in the shake list are used.
 * If bSettle && bShakeOnly the settles are used too.
 */
extern void done_graph(t_graph *g);
/* Free the memory in g */
 
extern void p_graph(FILE *log,char *title,t_graph *g);
/* Print a graph to log */

extern void mk_mshift(FILE *log,t_graph *g,matrix box,rvec x[]);
/* Calculate the mshift codes, based on the connection graph in g. */

extern void shift_atom(t_graph *g,matrix box,rvec x[],rvec x_s[],atom_id ai);
/* Shift single atom ai */

extern void shift_x(t_graph *g,matrix box,rvec x[],rvec x_s[]);
/* Add the shift vector to x, and store in x_s (may be same array as x) */

extern void shift_self(t_graph *g,matrix box,rvec x[]);
/* Id. but in place */

extern void unshift_atom(t_graph *g,matrix box,rvec x[],rvec x_s[],atom_id ai);
/* Unshift single atom ai */

extern void unshift_x(t_graph *g,matrix box,rvec x[],rvec x_s[]);
/* Subtract the shift vector from x_s, and store in x (may be same array) */

extern void unshift_self(t_graph *g,matrix box,rvec x[]);
/* Id, but in place */

#endif	/* _mshift_h */
