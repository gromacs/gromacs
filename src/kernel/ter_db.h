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

#ifndef _ter_db_h
#define _ter_db_h

static char *SRCID_ter_db_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) ter_db.h 1.16 9/30/97"
#endif /* HAVE_IDENT */

#include "sysstuff.h"
#include "hackblock.h"

extern int read_ter_db(char *inf,t_hackblock **tbptr,t_atomtype *atype);
/* Read database for N&C terminal hacking */

extern t_hackblock *choose_ter(int nb,t_hackblock tb[],char *title);
/* Interactively select one.. */

extern void print_ter_db(FILE *out,int nb,t_hackblock tb[],t_atomtype *atype);
/* Print the stuff */

#endif	/* _ter_db_h */
