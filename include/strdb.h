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
 * Good gRace! Old Maple Actually Chews Slate
 */

#ifndef _strdb_h
#define _strdb_h

static char *SRCID_strdb_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) strdb.h 1.9 2/2/97"
#endif /* HAVE_IDENT */

extern int fget_lines(FILE *in,char ***strings);
/* Read an array of lines from file in. strings should be
 * the address of an array of strings (to be malloced by this routine)
 * return the number of strings.
 */
extern int get_lines(char *db,char ***strings);
/* Open file db, or if non-existant file $GMXLIB/db and read strings 
 * return the number of strings.
 */

extern int search_str(int nstr,char **str,char *key);
/* Search an array of strings for key, return the index if found
 * -1 if not found.
 */

extern int get_strings(char *db,char ***strings);
/* Read an array of strings from file db or $GMXLIB/db. strings should be
 * the address of an array of strings (to be malloced by this routine)
 * return the number of strings.
 */
extern int get_file(char *db,char ***strings);
/* Read an array of strings from file db or $GMXLIB/db. strings should be
 * the address of an array of strings (to be malloced by this routine)
 * Does not need number of lines as first line in the file. 
 * return the number of strings.
 */

#endif	/* _strdb_h */
