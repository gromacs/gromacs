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

#include <filesystem>
#include <optional>
#include <string>
#include <vector>

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"

struct t_atoms;

//! An index group consisting or a name and list of atom indices
struct IndexGroup
{
    //! The name of the group, should not be empty
    std::string name;
    //! List of atom indices
    std::vector<int> particleIndices;
};

void check_index(const char* gname, int n, int index[], const char* traj, int natoms);
/* Checks if any index is smaller than zero or larger than natoms,
 * if so a fatal_error is given with the gname (if gname=NULL, "Index" is used)
 * and traj (if traj=NULL, "the trajectory" is used).
 */

/* Returns a list of atom index groups read from gfile */
std::vector<IndexGroup> init_index(const std::filesystem::path& gfile);

// DEPRECATED version of the above
std::vector<IndexGroup> init_index(const char* gfile);

void rd_index(const std::filesystem::path& statfile, int ngrps, int isize[], int* index[], char* grpnames[]);
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

// DEPRECATED version of the above
void rd_index(const char* statfile, int ngrps, int isize[], int* index[], char* grpnames[]);

void get_index(const t_atoms*                              atoms,
               const std::optional<std::filesystem::path>& fnm,
               int                                         ngrps,
               int                                         isize[],
               int*                                        index[],
               char*                                       grpnames[]);
/* Does the same as rd_index, but if the fnm has no value it
 * will not read from fnm, but it will make default index groups
 * for the atoms in *atoms.
 */

void get_index(const t_atoms* atoms, const char* fnm, int ngrps, int isize[], int* index[], char* grpnames[]);
/* DEPRECATED, use the above function instead.
 *
 * Does the same as rd_index, but if the fnm pointer is NULL it
 * will not read from fnm, but it will make default index groups
 * for the atoms in *atoms.
 */

struct t_cluster_ndx
{
    int                     maxframe = -1;
    std::vector<IndexGroup> clusters;
    std::vector<int>        inv_clust;
};

t_cluster_ndx cluster_index(FILE* fplog, const char* ndx);


/*! \brief Writes index groups to outf (writes an indexfile)
 *
 * \param[in] outf         Name of file to write to
 * \param[in] indexGroups  The index groups to write
 * \param[in] duplicate    Whether to write a duplicate of the groups after the normal groups, with
 * indices offset by \p numAtoms \param[in] numAtoms     The offset used with \p duplicate == true
 */
void write_index(const char* outf, gmx::ArrayRef<const IndexGroup> indexGroups, bool duplicate, int numAtoms);

//! Generates and returns standard index groups, bASK=FALSE gives default groups.
std::vector<IndexGroup> analyse(const t_atoms* atoms, gmx_bool bASK, gmx_bool bVerb);

/*! \brief Look up a group in a list of index groups.
 *
 * \param[inout] s         The string to look up
 * \param[in] indexGroups  The index groups
 */
int find_group(const char* s, gmx::ArrayRef<const IndexGroup> indexGroups);

/*! \brief Look up a group in a legacy list of index groups.
 *
 * \param[inout] s    The string to look up
 * \param[in] ngrps   The number of groups
 * \param[in] grpname The names of the groups
 * \return the group number or -1 if not found.
 */
int find_group(const char* s, int ngrps, char** grpname);

#endif
