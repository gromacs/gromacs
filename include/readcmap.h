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
 * GRowing Old MAkes el Chrono Sweat
 */

#ifndef	_readcmap_h
#define	_readcmap_h

#ifdef HAVE_IDENT
#ident	"@(#) readcmap.h 1.17 8/25/97"
#endif /* HAVE_IDENT */
#include "typedefs.h"
#include "writeps.h"
	
extern int searchcmap(int n,t_mapping map[],char c);
/* Seach in the map for code 'c' and return entry number. 
 * return -1 if not found
 */

extern int readcmap(char *fn,t_mapping **map);
/* Read the mapping table from fn, return number of entries */
extern void read2cmap(char *fn,
		      t_mapping **map1,int *n1,t_mapping **map2,int *n2);
/* Read 2 mapping tables from fn */

extern void printcmap(FILE *out,int n,t_mapping map[]);
/* print mapping table to out */

extern void writecmap(char *fn,int n,t_mapping map[]);
/* print mapping table to fn */

#endif	/* _readcmap_h */
