/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 1.6
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
 * S  C  A  M  O  R  G
 */

#ifndef _ter_db_h
#define _ter_db_h

static char *SRCID_ter_db_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) ter_db.h 1.16 9/30/97"
#endif /* HAVE_IDENT */

#include "sysstuff.h"
#include "pdb2gmx.h"

extern int read_ter_db(char *inf,t_terblock **tbptr,t_atomtype *atype);
/* Read database for N&C terminal hacking */

extern t_terblock *choose_ter(int nb,t_terblock tb[],char *title);
/* Interactively select one.. */

extern void print_ter_db(FILE *out,int nb,t_terblock tb[],t_atomtype *atype);
/* Print the stuff */

#endif	/* _ter_db_h */
