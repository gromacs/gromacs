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
 * GRowing Old MAkes el Chrono Sweat
 */

#ifndef _toputil_h
#define _toputil_h

static char *SRCID_toputil_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) toputil.h 1.25 19 Nov 1995"
#endif /* HAVE_IDENT */

#include "grompp.h"

/* UTILITIES */

extern int at2type(char *str, t_atomtype *at);

extern char *type2nm(int nt, t_atomtype *at);

extern void pr_alloc (int extra, t_params *pr);

extern void set_p_string(t_param *p,char *s);

/* INITIATE */

extern void init_plist(t_params plist[]);

extern void init_atomtype (t_atomtype *at);

extern void init_molinfo(t_molinfo *mol);

extern void init_top  (t_topology *top);

extern void done_top(t_topology *top);

/* FREE */
extern void done_block(t_block *block);

extern void done_top(t_topology *top);

extern void done_atom (t_atoms *at);

extern void done_mi(t_molinfo *mi);

/* PRINTING */

extern void print_block(FILE *out,char *szName,char *szIndex, 
			char *szA,t_block *block);

extern void print_atoms(FILE *out,t_atomtype *atype,t_atoms *at,int *cgnr);

extern void print_bondeds(FILE *out,int natoms,directive d,
			  int ftype,int fsubtype,t_params plist[]);

extern void print_excl(FILE *out, int natoms, t_excls excls[]);

#endif	/* _toputil_h */
