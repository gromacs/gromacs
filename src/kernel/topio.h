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

#ifndef _topio_h
#define _topio_h

static char *SRCID_topio_h = "$Id$";

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

extern real check_mol(t_atoms *atoms);
/* Check mass and charge */

extern void preprocess(char *infile, 
		       char *outfile,
		       char *cpp,
		       char *define,
		       char *include);

extern char **do_top(bool         bVerbose,
		     char         *topfile,
		     char         *topppfile,
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
