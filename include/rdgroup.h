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
 * Green Red Orange Magenta Azure Cyan Skyblue
 */

#ifndef _rdgroup_h
#define _rdgroup_h

static char *SRCID_rdgroup_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) rdgroup.h 1.15 2/2/97"
#endif /* HAVE_IDENT */
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

#ifdef CPLUSPLUS
}
#endif

#endif	/* _rdgroup_h */
