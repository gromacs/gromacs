/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

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
 * Gromacs Runs On Most of All Computer Systems
 */

#ifndef _index_h
#define _index_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <typedefs.h>

#ifdef CPLUSPLUS
extern "C" { 
#endif

extern void check_index(char *gname,int n,atom_id index[],
			char *traj,int natoms);
/* Checks if any index is smaller than zero or larger than natoms,
 * if so a fatal_error is given with the gname (if gname=NULL, "Index" is used)
 * and traj (if traj=NULL, "the trajectory" is used).
 */

t_block *init_index(char *gfile, char ***grpname);
/* Lower level routine than the next */

void rd_index(char *statfile,int ngrps,int isize[],
	      atom_id *index[],char *grpnames[]);
/* Assume the group file is generated, so the
 * format need not be user-friendly. The format is:
 * nr of groups, total nr of atoms
 * for each group: name nr of element, elements
 * The function opens a file, reads ngrps groups, puts the
 * sizes in isize, the atom_id s in index and the names of
 * the groups in grpnames.
 *
 * It is also assumed, that when ngrps groups are requested
 * memory has been allocated for ngrps index arrays, and that
 * the dimension of the isize and grpnames arrays are ngrps.
 */
 
void rd_index_nrs(char *statfile,int ngrps,int isize[],
		  atom_id *index[],char *grpnames[],int grpnr[]);
/* the same but also reads the number of the selected group*/

void get_index(t_atoms *atoms, char *fnm, int ngrps,
	       int isize[], atom_id *index[],char *grpnames[]);
/* Does the same as rd_index, but if the fnm pointer is NULL it
 * will not read from fnm, but it will make default index groups
 * for the atoms in *atoms.
 */ 

typedef struct {
  int     maxframe;
  char    **grpname;
  t_block *clust;
  atom_id *inv_clust;
} t_cluster_ndx;

extern t_cluster_ndx *cluster_index(char *ndx);
  
typedef struct {
  int n;
  char **aa;
} t_aa_names;

extern t_aa_names *get_aa_names(void);
/* Read the database in aminoacids.dat */

extern bool is_protein(t_aa_names *aan,char *resnm);
/* gives true if resnm occurs in aminoacids.dat */

extern void done_aa_names(t_aa_names **aan);
/* Free memory. Pass address of the pointer youget from get_aa_names */

extern t_block *new_block(void);
/* allocate new block */

extern void write_index(char *outf, t_block *b,char **gnames);
/* Writes index blocks to outf (writes an indexfile) */

void add_grp(t_block *b,char ***gnames,int nra,atom_id a[],const char *name);
/* Ads group a with name name to block b and namelist gnames */ 

extern void analyse(t_atoms *atoms,t_block *gb,char ***gn,
                    bool bASK,bool bVerb);
/* Makes index groups gb with names gn for atoms in atoms.
 * bASK=FALSE gives default groups.
 */

extern int find_group(char s[], int ngrps, char **grpname);


#ifdef CPLUSPLUS
}
#endif

#endif	/* _index_h */
