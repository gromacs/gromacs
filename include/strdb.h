/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.0
 * 
 * Copyright (c) 1991-2001
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * Do check out http://www.gromacs.org , or mail us at gromacs@gromacs.org .
 * 
 * And Hey:
 * Giving Russians Opium May Alter Current Situation
 */

#ifndef _strdb_h
#define _strdb_h

static char *SRCID_strdb_h = "$Id$";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_IDENT
#ident	"@(#) strdb.h 1.9 2/2/97"
#endif /* HAVE_IDENT */

#include "typedefs.h"

extern bool get_a_line(FILE *fp,char line[],int n);
/* Read a line of at most n characters form *fp to line. 
 * Comment ';...' and leading spaces are removed, empty lines are skipped.
 * Return FALSE when eof. 
 */

extern bool get_header(char line[],char header[]);
/* Read a header between '[' and ']' from line to header.
 * Returns FALSE no header is found.
 */

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
