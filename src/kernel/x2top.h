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
 * Giant Rising Ordinary Mutants for A Clerical Setup
 */

#ifndef _x2top_h
#define _x2top_h

static char *SRCID_x2top_h = "$Id$";
	
#include <stdio.h>
	
typedef struct {
  char *elem,*type;
  int  nbonds;
  char **bond;
} t_nm2type;

extern t_nm2type *rd_nm2type(char *ff,int *nnm);
/* Read the name 2 type database. nnm is the number of entries 
 * ff is the force field.
 */

extern void dump_nm2type(FILE *fp,int nnm,t_nm2type nm2t[]);
/* Dump the database for debugging. Can be reread by the program */

extern char *nm2type(int nnm,t_nm2type nm2t[],char *nm,int nbonds);
/* Try to determine the atomtype (force field dependent) for an element */

#endif
