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
 * GROtesk MACabre and Sinister
 */

#ifndef	_topio_h
#define	_topio_h

#ifdef HAVE_IDENT
#ident	"@(#) topio.h 1.46 9/30/97"
#endif /* HAVE_IDENT */

#include "typedefs.h"
#include "readir.h"
#include "grompp.h"

typedef struct {
  int whichmol;
  int nrcopies;
} t_simsystem;

extern void stupid_fill(t_block *grp, int maxf);

extern void preprocess(char *infile, 
		       char *outfile,
		       char *cpp,
		       char *define,
		       char *include);

extern char **do_top(bool         bVerbose,
		     char         *topol,
		     t_gromppopts *opts,
		     t_symtab     *symtab,
		     t_params     plist[],
		     t_atomtype   *atype,
		     int          *nrmols,
		     t_molinfo    **molinfo,
		     t_inputrec   *ir,
		     int          *nsim,
		     t_simsystem  **sims);

#endif	/* _topio_h */
