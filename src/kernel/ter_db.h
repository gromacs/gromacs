/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.1
 * Copyright (c) 1991-2001, University of Groningen, The Netherlands
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
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Gromacs Runs One Microsecond At Cannonball Speeds
 */

#ifndef _ter_db_h
#define _ter_db_h


#include "sysstuff.h"
#include "hackblock.h"
#include "grompp.h"

extern int read_ter_db(char *FF,char ter,
		       t_hackblock **tbptr,t_atomtype *atype);
/* Read database for N&C terminal hacking */

extern t_hackblock **filter_ter(int nb,t_hackblock tb[],char *resname,int *nret);
/* Return a list of pointers to blocks that match residue name */

extern t_hackblock *choose_ter(int nb,t_hackblock **tb,char *title);
/* Interactively select one.. */

extern void print_ter_db(FILE *out,int nb,t_hackblock tb[],t_atomtype *atype);
/* Print the stuff */

#endif	/* _ter_db_h */
