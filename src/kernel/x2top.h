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

extern char *nm2type(int nnm,t_nm2type nm2t[],char *elem,int nbonds);
/* Try to determine the atomtype (force field dependent) for an element */

#endif
