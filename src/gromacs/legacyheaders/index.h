/*
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

#include "typedefs.h"

#ifdef __cplusplus
extern "C" { 
#endif

void check_index(char *gname,int n,atom_id index[],
			char *traj,int natoms);
/* Checks if any index is smaller than zero or larger than natoms,
 * if so a fatal_error is given with the gname (if gname=NULL, "Index" is used)
 * and traj (if traj=NULL, "the trajectory" is used).
 */

t_blocka *init_index(const char *gfile, char ***grpname);
/* Lower level routine than the next */

void rd_index(const char *statfile,int ngrps,int isize[],
	      atom_id *index[],char *grpnames[]);
/* Assume the group file is generated, so the
 * format need not be user-friendly. The format is:
 * nr of groups, total nr of atoms
 * for each group: name nr of element, elements.
 *
 * The function opens a file, reads ngrps groups, asks the 
 * user for group numbers, and puts the resulting sizes in 
 * isize, the atom_id s in index and the names of
 * the groups in grpnames.
 *
 * It is also assumed, that when ngrps groups are requested
 * memory has been allocated for ngrps index arrays, and that
 * the dimension of the isize and grpnames arrays are ngrps.
 */
 
void rd_index_nrs(char *statfile,int ngrps,int isize[],
		  atom_id *index[],char *grpnames[],int grpnr[]);
/* the same but also reads the number of the selected group*/

void get_index(t_atoms *atoms, const char *fnm, int ngrps,
	       int isize[], atom_id *index[],char *grpnames[]);
/* Does the same as rd_index, but if the fnm pointer is NULL it
 * will not read from fnm, but it will make default index groups
 * for the atoms in *atoms.
 */ 

typedef struct {
  int      maxframe;
  char     **grpname;
  t_blocka *clust;
  atom_id  *inv_clust;
} t_cluster_ndx;

t_cluster_ndx *cluster_index(FILE *fplog,const char *ndx);
  
typedef struct {
  int n;
  char **name;
} t_names;

typedef struct gmx_residuetype *
gmx_residuetype_t;

int
gmx_residuetype_init(gmx_residuetype_t *rt);

int
gmx_residuetype_destroy(gmx_residuetype_t rt);

int
gmx_residuetype_get_type(gmx_residuetype_t rt,const char * resname, const char ** p_restype);

int
gmx_residuetype_add(gmx_residuetype_t rt,const char *newresname, const char *newrestype);

int
gmx_residuetype_get_alltypes(gmx_residuetype_t    rt,
                             const char ***       p_typenames,
                             int *                ntypes);

gmx_bool 
gmx_residuetype_is_protein(gmx_residuetype_t rt, const char *resnm);

gmx_bool 
gmx_residuetype_is_dna(gmx_residuetype_t rt, const char *resnm);

gmx_bool 
gmx_residuetype_is_rna(gmx_residuetype_t rt, const char *resnm);

int
gmx_residuetype_get_size(gmx_residuetype_t rt);

int
gmx_residuetype_get_index(gmx_residuetype_t rt, const char *resnm);

const char *
gmx_residuetype_get_name(gmx_residuetype_t rt, int index);






t_blocka *new_blocka(void);
/* allocate new block */

void write_index(const char *outf, t_blocka *b,char **gnames);
/* Writes index blocks to outf (writes an indexfile) */

void add_grp(t_blocka *b,char ***gnames,int nra,atom_id a[],const char *name);
/* Ads group a with name name to block b and namelist gnames */ 

void analyse(t_atoms *atoms,t_blocka *gb,char ***gn,
                    gmx_bool bASK,gmx_bool bVerb);
/* Makes index groups gb with names gn for atoms in atoms.
 * bASK=FALSE gives default groups.
 */

int find_group(char s[], int ngrps, char **grpname);


#ifdef __cplusplus
}
#endif

#endif	/* _index_h */
