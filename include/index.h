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
 * Getting the Right Output Means no Artefacts in Calculating Stuff
 */
static char *SRCID_index_h = "$Id$";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

extern t_block *new_block(void);
/* allocate new block */

extern void write_index(char *outf, t_block *b,char **gnames);
/* Writes index blocks to outf (writes an indexfile) */

void add_grp(t_block *b,char ***gnames,int nra,atom_id a[],char *name);
/* Ads group a with name name to block b and namelist gnames */ 

extern void analyse(t_atoms *atoms,t_block *gb,char ***gn,
                    bool bASK,bool bVerb);
/* Makes index groups gb with names gn for atoms in atoms.
 * bASK=FALSE gives default groups.
 */

extern bool is_protein(char *resnm);
/* gives true if resnm occurs in aminoacids.dat */
