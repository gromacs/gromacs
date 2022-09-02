/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
#ifndef GMX_TOPOLOGY_INDEX_H
#define GMX_TOPOLOGY_INDEX_H

#include <cstdio>

#include <vector>

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"

struct t_atoms;
struct t_blocka;

void check_index(const char* gname, int n, int index[], const char* traj, int natoms);
/* Checks if any index is smaller than zero or larger than natoms,
 * if so a fatal_error is given with the gname (if gname=NULL, "Index" is used)
 * and traj (if traj=NULL, "the trajectory" is used).
 */

struct t_blocka* init_index(const char* gfile, char*** grpname);
/* Lower level routine than the next */

void rd_index(const char* statfile, int ngrps, int isize[], int* index[], char* grpnames[]);
/* Assume the group file is generated, so the
 * format need not be user-friendly. The format is:
 * nr of groups, total nr of atoms
 * for each group: name nr of element, elements.
 *
 * The function opens a file, reads ngrps groups, asks the
 * user for group numbers, and puts the resulting sizes in
 * isize, the int s in index and the names of
 * the groups in grpnames.
 *
 * It is also assumed, that when ngrps groups are requested
 * memory has been allocated for ngrps index arrays, and that
 * the dimension of the isize and grpnames arrays are ngrps.
 */

void get_index(const t_atoms* atoms, const char* fnm, int ngrps, int isize[], int* index[], char* grpnames[]);
/* Does the same as rd_index, but if the fnm pointer is NULL it
 * will not read from fnm, but it will make default index groups
 * for the atoms in *atoms.
 */

struct t_cluster_ndx
{
    int              maxframe = -1;
    char**           grpname  = nullptr;
    struct t_blocka* clust    = nullptr;
    std::vector<int> inv_clust;
};

t_cluster_ndx cluster_index(FILE* fplog, const char* ndx);

//! Write index blocks to file.
void write_index(const char* outf, struct t_blocka* b, char** gnames, gmx_bool bDuplicate, int natoms);

/*! \brief
 * Add a new group with \p name to \p b.
 *
 * \param[in] b Block struct to add group to.
 * \param[in] gnames Names of groups.
 * \param[in] a Group to add to Block.
 * \param[in] name Group name.
 */
void add_grp(struct t_blocka* b, char*** gnames, gmx::ArrayRef<const int> a, const std::string& name);

/*! \brief
 * Builds index group \p gb form input \p atoms and \p gn names.
 *
 * \param[in] atoms Topology atoms to generate index groups for.
 * \param[out] gb Index datastructures to populate.
 * \param[out] gn Names of index groups.
 * \param[in] bASK If user should be prompted for groups.
 *                 If false, generates default groups.
 * \param[in] bVerb If output should be verbose.
 */
void analyse(const t_atoms* atoms, struct t_blocka* gb, char*** gn, gmx_bool bASK, gmx_bool bVerb);

/*! \brief Look up a group in a list.
 *
 * \param[inout] s    The string to look up
 * \param[in] ngrps   The number of groups
 * \param[in] grpname The names of the groups
 * \return the group number or -1 if not found.
 */
int find_group(const char* s, int ngrps, char** grpname);


#endif
