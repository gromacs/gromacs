/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2010,2014, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
#ifndef GMX_TOPOLOGY_INDEX_H
#define GMX_TOPOLOGY_INDEX_H

#include <stdio.h>

#include "gromacs/legacyheaders/types/simple.h"

#ifdef __cplusplus
extern "C" {
#endif

struct t_atoms;
struct t_blocka;

void check_index(char *gname, int n, atom_id index[],
                 char *traj, int natoms);
/* Checks if any index is smaller than zero or larger than natoms,
 * if so a fatal_error is given with the gname (if gname=NULL, "Index" is used)
 * and traj (if traj=NULL, "the trajectory" is used).
 */

struct t_blocka *init_index(const char *gfile, char ***grpname);
/* Lower level routine than the next */

void rd_index(const char *statfile, int ngrps, int isize[],
              atom_id *index[], char *grpnames[]);
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

void rd_index_nrs(char *statfile, int ngrps, int isize[],
                  atom_id *index[], char *grpnames[], int grpnr[]);
/* the same but also reads the number of the selected group*/

void get_index(struct t_atoms *atoms, const char *fnm, int ngrps,
               int isize[], atom_id *index[], char *grpnames[]);
/* Does the same as rd_index, but if the fnm pointer is NULL it
 * will not read from fnm, but it will make default index groups
 * for the atoms in *atoms.
 */

typedef struct {
    int               maxframe;
    char            **grpname;
    struct t_blocka  *clust;
    atom_id          *inv_clust;
} t_cluster_ndx;

t_cluster_ndx *cluster_index(FILE *fplog, const char *ndx);


void write_index(const char *outf, struct t_blocka *b, char **gnames, gmx_bool bDuplicate, int natoms);
/* Writes index blocks to outf (writes an indexfile) */

void add_grp(struct t_blocka *b, char ***gnames, int nra, atom_id a[], const char *name);
/* Ads group a with name name to block b and namelist gnames */

void analyse(struct t_atoms *atoms, struct t_blocka *gb, char ***gn,
             gmx_bool bASK, gmx_bool bVerb);
/* Makes index groups gb with names gn for atoms in atoms.
 * bASK=FALSE gives default groups.
 */

int find_group(char s[], int ngrps, char **grpname);


#ifdef __cplusplus
}
#endif

#endif
