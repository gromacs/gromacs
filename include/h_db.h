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
 * Grunge ROck MAChoS
 */

#ifndef	_h_db_h
#define	_h_db_h

#ifdef HAVE_IDENT
#ident	"@(#) h_db.h 1.10 2/2/97"
#endif /* HAVE_IDENT */
#include "sysstuff.h"
#include "pdb2gmx.h"

/* functions for the h-database */

extern void read_ab(FILE *in,t_add_block *ab);
/* Read one add block */

extern int read_h_db(char *fn,t_addh **ah);
/* Read the database from file */

extern void print_ab(FILE *out,t_add_block *ab);
/* print one add block */

extern void print_h_db(FILE *out,int nh,t_addh ah[]);
/* Print the database to file */

extern int compaddh(const void *a,const void *b);

extern t_addh *search_h_db(int nh,t_addh ah[],char *key);
/* Search for an entry in the database */



#endif	/* _h_db_h */
