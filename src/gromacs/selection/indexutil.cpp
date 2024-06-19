/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2009- The GROMACS Authors
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
/*! \internal \file
 * \brief
 * Implements functions in indexutil.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_selection
 */
#include "gmxpre.h"

#include "gromacs/selection/indexutil.h"

#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <iterator>
#include <numeric>
#include <string>
#include <vector>

#include "gromacs/topology/atoms.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/index.h"
#include "gromacs/topology/mtop_lookup.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textwriter.h"

namespace gmx
{

IndexGroupsAndNames::IndexGroupsAndNames(gmx::ArrayRef<const IndexGroup> indexGroups) :
    indexGroups_(indexGroups.begin(), indexGroups.end())
{
}

bool IndexGroupsAndNames::containsGroupName(const std::string& groupName) const
{
    return std::any_of(
            std::begin(indexGroups_), std::end(indexGroups_), [&groupName](const IndexGroup& indexGroup) {
                return equalCaseInsensitive(groupName, indexGroup.name);
            });
}

std::vector<Index> IndexGroupsAndNames::indices(const std::string& groupName) const
{
    if (!containsGroupName(groupName))
    {
        GMX_THROW(
                InconsistentInputError(
                        std::string("Group ") + groupName
                        + " referenced in the .mdp file was not found in the list of index "
                          "groups.\n"
                          "Group names must match either [moleculetype] names or custom index "
                          "group\n"
                          "names, in which case you must supply an index file to the '-n' option\n"
                          "of grompp."));
    }
    const auto groupNamePosition = std::find_if(
            std::begin(indexGroups_), std::end(indexGroups_), [&groupName](const IndexGroup& indexGroup) {
                return equalCaseInsensitive(groupName, indexGroup.name);
            });
    const auto         groupIndex = std::distance(std::begin(indexGroups_), groupNamePosition);
    std::vector<Index> groupIndices(indexGroups_[groupIndex].particleIndices.begin(),
                                    indexGroups_[groupIndex].particleIndices.end());
    return groupIndices;
}

} // namespace gmx

/********************************************************************
 * gmx_ana_indexgrps_t functions
 ********************************************************************/

/*! \internal \brief
 * Stores a set of index groups.
 */
struct gmx_ana_indexgrps_t
{
    //! Initializes an empty set of groups.
    explicit gmx_ana_indexgrps_t(int nr) : g(nr) { names.reserve(nr); }
    ~gmx_ana_indexgrps_t()
    {
        for (auto& indexGrp : g)
        {
            gmx_ana_index_deinit(&indexGrp);
        }
    }


    /** Array of index groups. */
    std::vector<gmx_ana_index_t> g;
    /** Group names. */
    std::vector<std::string> names;
};

/*!
 * \param[out] g     Index group structure.
 * \param[in]  top   Topology structure.
 * \param[in]  fnm   File name for the index file.
 *   Memory is automatically allocated.
 *
 * One or both of \p top or \p fnm can be NULL.
 * If \p top is NULL, an index file is required and the groups are read
 * from the file (uses Gromacs routine init_index()).
 * If \p fnm is NULL, default groups are constructed based on the
 * topology (uses Gromacs routine analyse()).
 * If both are null, the index group structure is initialized empty.
 */
void gmx_ana_indexgrps_init(gmx_ana_indexgrps_t** g, gmx_mtop_t* top, const char* fnm)
{
    std::vector<IndexGroup> indexGroups;

    if (fnm)
    {
        indexGroups = init_index(fnm);
    }
    else if (top)
    {
        // TODO: Propagate mtop further.
        t_atoms atoms = gmx_mtop_global_atoms(*top);
        indexGroups   = analyse(&atoms, FALSE, FALSE);
        done_atom(&atoms);
    }
    else
    {
        *g = new gmx_ana_indexgrps_t(0);
        return;
    }

    *g = new gmx_ana_indexgrps_t(gmx::ssize(indexGroups));
    for (int i = 0; i < gmx::ssize(indexGroups); i++)
    {
        gmx::ArrayRef<const int> indexGroup = indexGroups[i].particleIndices;

        gmx_ana_index_t* grp = &(*g)->g[i];

        grp->isize = indexGroup.ssize();
        snew(grp->index, grp->isize);
        for (int j = 0; j < grp->isize; ++j)
        {
            grp->index[j] = indexGroup[j];
        }
        grp->nalloc_index = grp->isize;
        (*g)->names.emplace_back(indexGroups[i].name);
    }
}

/*!
 * \param[in] g  Index groups structure.
 *
 * The pointer \p g is invalid after the call.
 */
void gmx_ana_indexgrps_free(gmx_ana_indexgrps_t* g)
{
    delete g;
}


/*!
 * \param[out] dest     Output structure.
 * \param[out] destName Receives the name of the group if found.
 * \param[in]  src      Input index groups.
 * \param[in]  n        Number of the group to extract.
 * \returns true if \p n is a valid group in \p src, false otherwise.
 */
bool gmx_ana_indexgrps_extract(gmx_ana_index_t* dest, std::string* destName, gmx_ana_indexgrps_t* src, int n)
{
    destName->clear();
    if (n < 0 || n >= gmx::Index(src->g.size()))
    {
        dest->isize = 0;
        return false;
    }
    *destName = src->names[n];
    gmx_ana_index_copy(dest, &src->g[n], true);
    return true;
}

/*!
 * \param[out] dest     Output structure.
 * \param[out] destName Receives the name of the group if found.
 * \param[in]  src      Input index groups.
 * \param[in]  name     Name (or part of the name) of the group to extract.
 * \returns true if \p name is a valid group in \p src, false otherwise.
 *
 * Uses the Gromacs routine find_group() to find the actual group;
 * the comparison is case-insensitive.
 */
bool gmx_ana_indexgrps_find(gmx_ana_index_t* dest, std::string* destName, gmx_ana_indexgrps_t* src, const char* name)
{
    const char** names;

    destName->clear();
    snew(names, src->g.size());
    for (size_t i = 0; i < src->g.size(); ++i)
    {
        names[i] = src->names[i].c_str();
    }
    int n = find_group(const_cast<char*>(name), src->g.size(), const_cast<char**>(names));
    sfree(names);
    if (n < 0)
    {
        dest->isize = 0;
        return false;
    }

    return gmx_ana_indexgrps_extract(dest, destName, src, n);
}

/*!
 * \param[in]  writer Writer to use for output.
 * \param[in]  g      Index groups to print.
 * \param[in]  maxn   Maximum number of indices to print
 *      (-1 = print all, 0 = print only names).
 */
void gmx_ana_indexgrps_print(gmx::TextWriter* writer, gmx_ana_indexgrps_t* g, int maxn)
{
    for (gmx::Index i = 0; i < gmx::ssize(g->g); ++i)
    {
        writer->writeString(gmx::formatString(" Group %2zd \"%s\" ", i, g->names[i].c_str()));
        gmx_ana_index_dump(writer, &g->g[i], maxn);
    }
}

/********************************************************************
 * gmx_ana_index_t functions
 ********************************************************************/

/*!
 * \param[in,out] g      Index group structure.
 * \param[in]     isize  Maximum number of atoms to reserve space for.
 */
void gmx_ana_index_reserve(gmx_ana_index_t* g, int isize)
{
    if (g->nalloc_index < isize)
    {
        srenew(g->index, isize);
        g->nalloc_index = isize;
    }
}

/*!
 * \param[in,out] g      Index group structure.
 *
 * Resizes the memory allocated for holding the indices such that the
 * current contents fit.
 */
void gmx_ana_index_squeeze(gmx_ana_index_t* g)
{
    srenew(g->index, g->isize);
    g->nalloc_index = g->isize;
}

/*!
 * \param[out] g      Output structure.
 *
 * Any contents of \p g are discarded without freeing.
 */
void gmx_ana_index_clear(gmx_ana_index_t* g)
{
    g->isize        = 0;
    g->index        = nullptr;
    g->nalloc_index = 0;
}

/*!
 * \param[out] g      Output structure.
 * \param[in]  isize  Number of atoms in the new group.
 * \param[in]  index  Array of \p isize atoms (can be NULL if \p isize is 0).
 * \param[in]  nalloc Number of elements allocated for \p index
 *   (if 0, \p index is not freed in gmx_ana_index_deinit())
 *
 * No copy if \p index is made.
 */
void gmx_ana_index_set(gmx_ana_index_t* g, int isize, int* index, int nalloc)
{
    g->isize        = isize;
    g->index        = index;
    g->nalloc_index = nalloc;
}

/*!
 * \param[out] g      Output structure.
 * \param[in]  natoms Number of atoms.
 */
void gmx_ana_index_init_simple(gmx_ana_index_t* g, int natoms)
{
    int i;

    g->isize = natoms;
    snew(g->index, natoms);
    for (i = 0; i < natoms; ++i)
    {
        g->index[i] = i;
    }
    g->nalloc_index = natoms;
}

/*!
 * \param[in] g  Index group structure.
 *
 * The pointer \p g is not freed.
 */
void gmx_ana_index_deinit(gmx_ana_index_t* g)
{
    if (g->nalloc_index > 0)
    {
        sfree(g->index);
    }
    gmx_ana_index_clear(g);
}

/*!
 * \param[out] dest   Destination index group.
 * \param[in]  src    Source index group.
 * \param[in]  bAlloc If true, memory is allocated at \p dest; otherwise,
 *   it is assumed that enough memory has been allocated for index.
 */
void gmx_ana_index_copy(gmx_ana_index_t* dest, gmx_ana_index_t* src, bool bAlloc)
{
    dest->isize = src->isize;
    if (bAlloc)
    {
        snew(dest->index, dest->isize);
        dest->nalloc_index = dest->isize;
    }
    if (dest->isize > 0)
    {
        std::memcpy(dest->index, src->index, dest->isize * sizeof(*dest->index));
    }
}

/*!
 * \param[in]  writer Writer to use for output.
 * \param[in]  g      Index group to print.
 * \param[in]  maxn   Maximum number of indices to print (-1 = print all).
 */
void gmx_ana_index_dump(gmx::TextWriter* writer, gmx_ana_index_t* g, int maxn)
{
    writer->writeString(gmx::formatString("(%d atoms)", g->isize));
    if (maxn != 0)
    {
        writer->writeString(":");
        int n = g->isize;
        if (maxn >= 0 && n > maxn)
        {
            n = maxn;
        }
        for (int j = 0; j < n; ++j)
        {
            writer->writeString(gmx::formatString(" %d", g->index[j] + 1));
        }
        if (n < g->isize)
        {
            writer->writeString(" ...");
        }
    }
    writer->ensureLineBreak();
}

int gmx_ana_index_get_max_index(gmx_ana_index_t* g)
{
    if (g->isize == 0)
    {
        return 0;
    }
    else
    {
        return *std::max_element(g->index, g->index + g->isize);
    }
}

/*!
 * \param[in]  g      Index group to check.
 * \returns    true if the index group is sorted and has no duplicates,
 *   false otherwise.
 */
bool gmx_ana_index_check_sorted(gmx_ana_index_t* g)
{
    int i;

    for (i = 0; i < g->isize - 1; ++i)
    {
        if (g->index[i + 1] <= g->index[i])
        {
            return false;
        }
    }
    return true;
}

bool gmx_ana_index_check_range(gmx_ana_index_t* g, int natoms)
{
    for (int i = 0; i < g->isize; ++i)
    {
        if (g->index[i] < 0 || g->index[i] >= natoms)
        {
            return false;
        }
    }
    return true;
}

/********************************************************************
 * Set operations
 ********************************************************************/

/*!
 * \param[in,out] g  Index group to be sorted.
 */
void gmx_ana_index_sort(gmx_ana_index_t* g)
{
    std::sort(g->index, g->index + g->isize);
}

void gmx_ana_index_remove_duplicates(gmx_ana_index_t* g)
{
    int j = 0;
    for (int i = 0; i < g->isize; ++i)
    {
        if (i == 0 || g->index[i - 1] != g->index[i])
        {
            g->index[j] = g->index[i];
            ++j;
        }
    }
    g->isize = j;
}

/*!
 * \param[in]  a      Index group to check.
 * \param[in]  b      Index group to check.
 * \returns    true if \p a and \p b are equal, false otherwise.
 */
bool gmx_ana_index_equals(gmx_ana_index_t* a, gmx_ana_index_t* b)
{
    int i;

    if (a->isize != b->isize)
    {
        return false;
    }
    for (i = 0; i < a->isize; ++i)
    {
        if (a->index[i] != b->index[i])
        {
            return false;
        }
    }
    return true;
}

/*!
 * \param[in]  a      Index group to check against.
 * \param[in]  b      Index group to check.
 * \returns    true if \p b is contained in \p a,
 *   false otherwise.
 *
 * If the elements are not in the same order in both groups, the function
 * fails. However, the groups do not need to be sorted.
 */
bool gmx_ana_index_contains(gmx_ana_index_t* a, gmx_ana_index_t* b)
{
    int i, j;

    for (i = j = 0; j < b->isize; ++i, ++j)
    {
        while (i < a->isize && a->index[i] != b->index[j])
        {
            ++i;
        }
        if (i == a->isize)
        {
            return false;
        }
    }
    return true;
}

/*!
 * \param[out] dest Output index group (the intersection of \p a and \p b).
 * \param[in]  a    First index group.
 * \param[in]  b    Second index group.
 *
 * \p dest can be the same as \p a or \p b.
 */
void gmx_ana_index_intersection(gmx_ana_index_t* dest, gmx_ana_index_t* a, gmx_ana_index_t* b)
{
    int i, j, k;

    for (i = j = k = 0; i < a->isize && j < b->isize; ++i)
    {
        while (j < b->isize && b->index[j] < a->index[i])
        {
            ++j;
        }
        if (j < b->isize && b->index[j] == a->index[i])
        {
            dest->index[k++] = b->index[j++];
        }
    }
    dest->isize = k;
}

/*!
 * \param[out] dest Output index group (the difference \p a - \p b).
 * \param[in]  a    First index group.
 * \param[in]  b    Second index group.
 *
 * \p dest can equal \p a, but not \p b.
 */
void gmx_ana_index_difference(gmx_ana_index_t* dest, gmx_ana_index_t* a, gmx_ana_index_t* b)
{
    int i, j, k;

    for (i = j = k = 0; i < a->isize; ++i)
    {
        while (j < b->isize && b->index[j] < a->index[i])
        {
            ++j;
        }
        if (j == b->isize || b->index[j] != a->index[i])
        {
            dest->index[k++] = a->index[i];
        }
    }
    dest->isize = k;
}

/*!
 * \param[in]  a    First index group.
 * \param[in]  b    Second index group.
 * \returns    Size of the difference \p a - \p b.
 */
int gmx_ana_index_difference_size(gmx_ana_index_t* a, gmx_ana_index_t* b)
{
    int i, j, k;

    for (i = j = k = 0; i < a->isize; ++i)
    {
        while (j < b->isize && b->index[j] < a->index[i])
        {
            ++j;
        }
        if (j == b->isize || b->index[j] != a->index[i])
        {
            ++k;
        }
    }
    return k;
}

/*!
 * \param[out] dest1 Output group 1 (will equal \p g).
 * \param[out] dest2 Output group 2 (will equal \p src - \p g).
 * \param[in]  src   Group to be partitioned.
 * \param[in]  g     One partition.
 *
 * \pre \p g is a subset of \p src and both sets are sorted
 * \pre \p dest1 has allocated storage to store \p src
 * \post \p dest1 == \p g
 * \post \p dest2 == \p src - \p g
 *
 * No storage should be allocated for \p dest2; after the call,
 * \p dest2->index points to the memory allocated for \p dest1
 * (to a part that is not used by \p dest1).
 *
 * The calculation can be performed in-place by setting \p dest1 equal to
 * \p src.
 */
void gmx_ana_index_partition(gmx_ana_index_t* dest1, gmx_ana_index_t* dest2, gmx_ana_index_t* src, gmx_ana_index_t* g)
{
    int i, j, k;

    dest2->index = dest1->index + g->isize;
    dest2->isize = src->isize - g->isize;
    for (i = g->isize - 1, j = src->isize - 1, k = dest2->isize - 1; i >= 0; --i, --j)
    {
        while (j >= 0 && src->index[j] != g->index[i])
        {
            dest2->index[k--] = src->index[j--];
        }
    }
    while (j >= 0)
    {
        dest2->index[k--] = src->index[j--];
    }
    gmx_ana_index_copy(dest1, g, false);
}

/*!
 * \param[out] dest Output index group (the union of \p a and \p b).
 * \param[in]  a    First index group.
 * \param[in]  b    Second index group.
 *
 * \p a and \p b can have common items.
 * \p dest can equal \p a or \p b.
 *
 * \see gmx_ana_index_merge()
 */
void gmx_ana_index_union(gmx_ana_index_t* dest, gmx_ana_index_t* a, gmx_ana_index_t* b)
{
    int dsize;
    int i, j, k;

    dsize       = gmx_ana_index_difference_size(b, a);
    i           = a->isize - 1;
    j           = b->isize - 1;
    dest->isize = a->isize + dsize;
    for (k = dest->isize - 1; k >= 0; k--)
    {
        if (i < 0 || (j >= 0 && a->index[i] < b->index[j]))
        {
            dest->index[k] = b->index[j--];
        }
        else
        {
            if (j >= 0 && a->index[i] == b->index[j])
            {
                --j;
            }
            dest->index[k] = a->index[i--];
        }
    }
}

void gmx_ana_index_union_unsorted(gmx_ana_index_t* dest, gmx_ana_index_t* a, gmx_ana_index_t* b)
{
    if (gmx_ana_index_check_sorted(b))
    {
        gmx_ana_index_union(dest, a, b);
    }
    else
    {
        gmx_ana_index_t tmp;
        gmx_ana_index_copy(&tmp, b, true);
        gmx_ana_index_sort(&tmp);
        gmx_ana_index_remove_duplicates(&tmp);
        gmx_ana_index_union(dest, a, &tmp);
        gmx_ana_index_deinit(&tmp);
    }
}

/*!
 * \param[out] dest Output index group (the union of \p a and \p b).
 * \param[in]  a    First index group.
 * \param[in]  b    Second index group.
 *
 * \p a and \p b should not have common items.
 * \p dest can equal \p a or \p b.
 *
 * \see gmx_ana_index_union()
 */
void gmx_ana_index_merge(gmx_ana_index_t* dest, gmx_ana_index_t* a, gmx_ana_index_t* b)
{
    int i, j, k;

    i           = a->isize - 1;
    j           = b->isize - 1;
    dest->isize = a->isize + b->isize;
    for (k = dest->isize - 1; k >= 0; k--)
    {
        if (i < 0 || (j >= 0 && a->index[i] < b->index[j]))
        {
            dest->index[k] = b->index[j--];
        }
        else
        {
            dest->index[k] = a->index[i--];
        }
    }
}

/********************************************************************
 * gmx_ana_indexmap_t and related things
 ********************************************************************/

/*! \brief
 * Helper for splitting a sequence of atom indices into groups.
 *
 * \param[in]     atomIndex  Index of the next atom in the sequence.
 * \param[in]     top        Topology structure.
 * \param[in]     type       Type of group to split into.
 * \param[in,out] id         Variable to receive the next group id.
 * \returns  `true` if \p atomIndex starts a new group in the sequence, i.e.,
 *     if \p *id was changed.
 *
 * \p *id should be initialized to `-1` before first call of this function, and
 * then each atom index in the sequence passed to the function in turn.
 *
 * \ingroup module_selection
 */
static bool next_group_index(int atomIndex, const gmx_mtop_t* top, e_index_t type, int* id)
{
    int prev = *id;
    switch (type)
    {
        case INDEX_ATOM: *id = atomIndex; break;
        case INDEX_RES:
        {
            int resind, molb = 0;
            mtopGetAtomAndResidueName(*top, atomIndex, &molb, nullptr, nullptr, nullptr, &resind);
            *id = resind;
            break;
        }
        case INDEX_MOL:
        {
            int molb = 0;
            *id      = mtopGetMoleculeIndex(*top, atomIndex, &molb);
            break;
        }
        case INDEX_UNKNOWN:
        case INDEX_ALL: *id = 0; break;
    }
    return prev != *id;
}

/*!
 * \param[in,out] t    Output block.
 * \param[in]  top  Topology structure
 *   (only used if \p type is \ref INDEX_RES or \ref INDEX_MOL, can be NULL
 *   otherwise).
 * \param[in]  g    Index group
 *   (can be NULL if \p type is \ref INDEX_UNKNOWN).
 * \param[in]  type Type of partitioning to make.
 * \param[in]  bComplete
 *   If true, the index group is expanded to include any residue/molecule
 *   (depending on \p type) that is partially contained in the group.
 *   If \p type is not INDEX_RES or INDEX_MOL, this has no effect.
 *
 * \p m should have been initialized somehow (calloc() is enough).
 * \p g should be sorted.
 */
void gmx_ana_index_make_block(t_blocka* t, const gmx_mtop_t* top, gmx_ana_index_t* g, e_index_t type, bool bComplete)
{
    if (type == INDEX_UNKNOWN)
    {
        sfree(t->a);
        srenew(t->index, 2);
        t->nr           = 1;
        t->nalloc_index = 2;
        t->index[0]     = 0;
        t->index[1]     = 0;
        t->nra          = 0;
        t->a            = nullptr;
        t->nalloc_a     = 0;
        return;
    }

    // TODO: Check callers and either check these there as well, or turn these
    // into exceptions.
    GMX_RELEASE_ASSERT(top != nullptr || (type != INDEX_RES && type != INDEX_MOL),
                       "Topology must be provided for residue or molecule blocks");
    GMX_RELEASE_ASSERT(type != INDEX_MOL || top->haveMoleculeIndices,
                       "Molecule information must be present for molecule blocks");

    /* bComplete only does something for INDEX_RES or INDEX_MOL, so turn it
     * off otherwise. */
    if (type != INDEX_RES && type != INDEX_MOL)
    {
        bComplete = false;
    }
    /* Allocate memory for the atom array and fill it unless we are using
     * completion. */
    if (bComplete)
    {
        t->nra = 0;
        /* We may allocate some extra memory here because we don't know in
         * advance how much will be needed. */
        if (t->nalloc_a < top->natoms)
        {
            srenew(t->a, top->natoms);
            t->nalloc_a = top->natoms;
        }
    }
    else
    {
        t->nra = g->isize;
        if (t->nalloc_a < g->isize)
        {
            srenew(t->a, g->isize);
            t->nalloc_a = g->isize;
        }
        if (t->nra > 0)
        {
            std::memcpy(t->a, g->index, g->isize * sizeof(*(t->a)));
        }
    }

    /* Allocate memory for the block index. We don't know in advance
     * how much will be needed, so we allocate some extra and free it in the
     * end. */
    if (t->nalloc_index < g->isize + 1)
    {
        srenew(t->index, g->isize + 1);
        t->nalloc_index = g->isize + 1;
    }
    /* Clear counters */
    t->nr    = 0;
    int id   = -1;
    int molb = 0;
    for (int i = 0; i < g->isize; ++i)
    {
        const int ai = g->index[i];
        /* Find the ID number of the atom/residue/molecule corresponding to
         * the atom. */
        if (next_group_index(ai, top, type, &id))
        {
            /* If this is the first atom in a new block, initialize the block. */
            if (bComplete)
            {
                /* For completion, we first set the start of the block. */
                t->index[t->nr++] = t->nra;
                /* And then we find all the atoms that should be included. */
                switch (type)
                {
                    case INDEX_RES:
                    {
                        int molnr, atnr_mol;
                        mtopGetMolblockIndex(*top, ai, &molb, &molnr, &atnr_mol);
                        const t_atoms& mol_atoms    = top->moltype[top->molblock[molb].type].atoms;
                        int            last_atom    = atnr_mol + 1;
                        const int      currentResid = mol_atoms.atom[atnr_mol].resind;
                        while (last_atom < mol_atoms.nr && mol_atoms.atom[last_atom].resind == currentResid)
                        {
                            ++last_atom;
                        }
                        int first_atom = atnr_mol - 1;
                        while (first_atom >= 0 && mol_atoms.atom[first_atom].resind == currentResid)
                        {
                            --first_atom;
                        }
                        const MoleculeBlockIndices& molBlock = top->moleculeBlockIndices[molb];
                        int                         first_mol_atom = molBlock.globalAtomStart;
                        first_mol_atom += molnr * molBlock.numAtomsPerMolecule;
                        first_atom = first_mol_atom + first_atom + 1;
                        last_atom  = first_mol_atom + last_atom - 1;
                        for (int j = first_atom; j <= last_atom; ++j)
                        {
                            t->a[t->nra++] = j;
                        }
                        break;
                    }
                    case INDEX_MOL:
                    {
                        int molnr, atnr_mol;
                        mtopGetMolblockIndex(*top, ai, &molb, &molnr, &atnr_mol);
                        const MoleculeBlockIndices& blockIndices = top->moleculeBlockIndices[molb];
                        const int                   atomStart    = blockIndices.globalAtomStart
                                              + (id - blockIndices.moleculeIndexStart)
                                                        * blockIndices.numAtomsPerMolecule;
                        for (int j = 0; j < blockIndices.numAtomsPerMolecule; ++j)
                        {
                            t->a[t->nra++] = atomStart + j;
                        }
                        break;
                    }
                    default: /* Should not be reached */
                        GMX_RELEASE_ASSERT(false, "Unreachable code was reached");
                        break;
                }
            }
            else
            {
                /* If not using completion, simply store the start of the block. */
                t->index[t->nr++] = i;
                if (type == INDEX_ATOM)
                {
                    id = -1;
                }
            }
        }
    }
    /* Set the end of the last block */
    t->index[t->nr] = t->nra;
    /* Free any unnecessary memory */
    srenew(t->index, t->nr + 1);
    t->nalloc_index = t->nr + 1;
    if (bComplete)
    {
        srenew(t->a, t->nra);
        t->nalloc_a = t->nra;
    }
}

/*!
 * \param[in] g   Index group to check.
 * \param[in] b   Block data to check against.
 * \returns   true if \p g consists of one or more complete blocks from \p b,
 *   false otherwise.
 *
 * The atoms in \p g are assumed to be sorted.
 */
bool gmx_ana_index_has_full_blocks(const gmx_ana_index_t* g, const gmx::RangePartitioning* b)
{
    int i, j, bi;

    i = bi = 0;
    /* Each round in the loop matches one block */
    while (i < g->isize)
    {
        /* Find the block that begins with the first unmatched atom */
        while (bi < b->numBlocks() && *b->block(bi).begin() != g->index[i])
        {
            ++bi;
        }
        /* If not found, or if too large, return */
        if (bi == b->numBlocks() || i + b->block(bi).size() > g->isize)
        {
            return false;
        }
        /* Check that the block matches the index */
        for (j = *b->block(bi).begin(); j < *b->block(bi).end(); ++j, ++i)
        {
            if (g->index[i] != j)
            {
                return false;
            }
        }
        /* Move the search to the next block */
        ++bi;
    }
    return true;
}

/*!
 * \param[in] g   Index group to check.
 * \param[in] b   Block data to check against.
 * \returns   true if \p g consists of one or more complete blocks from \p b,
 *   false otherwise.
 *
 * The atoms in \p g and \p b->a are assumed to be in the same order.
 */
bool gmx_ana_index_has_full_ablocks(gmx_ana_index_t* g, t_blocka* b)
{
    int i, j, bi;

    i = bi = 0;
    /* Each round in the loop matches one block */
    while (i < g->isize)
    {
        /* Find the block that begins with the first unmatched atom */
        while (bi < b->nr && b->a[b->index[bi]] != g->index[i])
        {
            ++bi;
        }
        /* If not found, or if too large, return */
        if (bi == b->nr || i + b->index[bi + 1] - b->index[bi] > g->isize)
        {
            return false;
        }
        /* Check that the block matches the index */
        for (j = b->index[bi]; j < b->index[bi + 1]; ++j, ++i)
        {
            if (b->a[j] != g->index[i])
            {
                return false;
            }
        }
        /* Move the search to the next block */
        ++bi;
    }
    return true;
}

/*!
 * \brief Returns if an atom is at a residue boundary.
 *
 * \param[in]     top   Topology data.
 * \param[in]     a     Atom index to check, should be -1 <= \p a < top->natoms.
 * \param[in,out] molb  The molecule block of atom a
 * \returns       true if atoms \p a and \p a + 1 are in different residues, false otherwise.
 */
static bool is_at_residue_boundary(const gmx_mtop_t& top, int a, int* molb)
{
    if (a == -1 || a + 1 == top.natoms)
    {
        return true;
    }
    int resindA;
    mtopGetAtomAndResidueName(top, a, molb, nullptr, nullptr, nullptr, &resindA);
    int resindAPlusOne;
    mtopGetAtomAndResidueName(top, a + 1, molb, nullptr, nullptr, nullptr, &resindAPlusOne);
    return resindAPlusOne != resindA;
}

/*!
 * \param[in] g     Index group to check.
 * \param[in] type  Block data to check against.
 * \param[in] top   Topology data.
 * \returns   true if \p g consists of one or more complete elements of type
 *   \p type, false otherwise.
 *
 * \p g is assumed to be sorted, otherwise may return false negatives.
 *
 * If \p type is \ref INDEX_ATOM, the return value is always true.
 * If \p type is \ref INDEX_UNKNOWN or \ref INDEX_ALL, the return value is
 * always false.
 */
bool gmx_ana_index_has_complete_elems(gmx_ana_index_t* g, e_index_t type, const gmx_mtop_t* top)
{
    if (g->isize == 0)
    {
        return true;
    }

    // TODO: Consider whether unsorted groups need to be supported better.
    switch (type)
    {
        case INDEX_UNKNOWN:
        case INDEX_ALL: return false;

        case INDEX_ATOM: return true;

        case INDEX_RES:
        {
            int molb  = 0;
            int aPrev = -1;
            for (int i = 0; i < g->isize; ++i)
            {
                const int a = g->index[i];
                // Check if a is consecutive or on a residue boundary
                if (a != aPrev + 1)
                {
                    if (!is_at_residue_boundary(*top, aPrev, &molb))
                    {
                        return false;
                    }
                    if (!is_at_residue_boundary(*top, a - 1, &molb))
                    {
                        return false;
                    }
                }
                aPrev = a;
            }
            GMX_ASSERT(g->isize > 0, "We return above when isize=0");
            const int a = g->index[g->isize - 1];
            if (!is_at_residue_boundary(*top, a, &molb))
            {
                return false;
            }
            break;
        }

        case INDEX_MOL:
        {
            auto molecules = gmx_mtop_molecules(*top);
            return gmx_ana_index_has_full_blocks(g, &molecules);
        }
    }
    return true;
}

/*!
 * \param[out] m      Output structure.
 *
 * Any contents of \p m are discarded without freeing.
 */
void gmx_ana_indexmap_clear(gmx_ana_indexmap_t* m)
{
    m->type              = INDEX_UNKNOWN;
    m->refid             = nullptr;
    m->mapid             = nullptr;
    m->mapb.nr           = 0;
    m->mapb.index        = nullptr;
    m->mapb.nalloc_index = 0;
    m->mapb.nra          = 0;
    m->mapb.a            = nullptr;
    m->mapb.nalloc_a     = 0;
    m->orgid             = nullptr;
    m->b.nr              = 0;
    m->b.index           = nullptr;
    m->b.nra             = 0;
    m->b.a               = nullptr;
    m->b.nalloc_index    = 0;
    m->b.nalloc_a        = 0;
    m->bStatic           = true;
}

/*!
 * \param[in,out] m      Mapping structure.
 * \param[in]     nr     Maximum number of blocks to reserve space for.
 * \param[in]     isize  Maximum number of atoms to reserve space for.
 */
void gmx_ana_indexmap_reserve(gmx_ana_indexmap_t* m, int nr, int isize)
{
    if (m->mapb.nalloc_index < nr + 1)
    {
        srenew(m->refid, nr);
        srenew(m->mapid, nr);
        srenew(m->orgid, nr);
        srenew(m->mapb.index, nr + 1);
        srenew(m->b.index, nr + 1);
        m->mapb.nalloc_index = nr + 1;
        m->b.nalloc_index    = nr + 1;
    }
    if (m->b.nalloc_a < isize)
    {
        srenew(m->b.a, isize);
        m->b.nalloc_a = isize;
    }
}

/*!
 * \param[in,out] m    Mapping structure to initialize.
 * \param[in]     g    Index group to map
 *   (can be NULL if \p type is \ref INDEX_UNKNOWN).
 * \param[in]     top  Topology structure
 *   (can be NULL if \p type is not \ref INDEX_RES or \ref INDEX_MOL).
 * \param[in]     type Type of mapping to construct.
 *
 * Initializes a new index group mapping.
 * The index group provided to gmx_ana_indexmap_update() should always be a
 * subset of the \p g given here.
 *
 * \p m should have been initialized somehow (calloc() is enough).
 */
void gmx_ana_indexmap_init(gmx_ana_indexmap_t* m, gmx_ana_index_t* g, const gmx_mtop_t* top, e_index_t type)
{
    m->type = type;
    gmx_ana_index_make_block(&m->b, top, g, type, false);
    gmx_ana_indexmap_reserve(m, m->b.nr, m->b.nra);
    int id = -1;
    for (int i = 0; i < m->b.nr; ++i)
    {
        const int ii = (type == INDEX_UNKNOWN ? 0 : m->b.a[m->b.index[i]]);
        next_group_index(ii, top, type, &id);
        m->refid[i] = i;
        m->mapid[i] = id;
        m->orgid[i] = id;
    }
    m->mapb.nr  = m->b.nr;
    m->mapb.nra = m->b.nra;
    m->mapb.a   = m->b.a;
    std::memcpy(m->mapb.index, m->b.index, (m->b.nr + 1) * sizeof(*(m->mapb.index)));
    m->bStatic = true;
}

int gmx_ana_indexmap_init_orgid_group(gmx_ana_indexmap_t* m, const gmx_mtop_t* top, e_index_t type)
{
    GMX_RELEASE_ASSERT(m->bStatic,
                       "Changing original IDs is not supported after starting "
                       "to use the mapping");
    GMX_RELEASE_ASSERT(top != nullptr || (type != INDEX_RES && type != INDEX_MOL),
                       "Topology must be provided for residue or molecule blocks");
    // Check that all atoms in each block belong to the same group.
    // This is a separate loop for better error handling (no state is modified
    // if there is an error.
    if (type == INDEX_RES || type == INDEX_MOL)
    {
        int id = -1;
        for (int i = 0; i < m->b.nr; ++i)
        {
            const int ii = m->b.a[m->b.index[i]];
            if (next_group_index(ii, top, type, &id))
            {
                for (int j = m->b.index[i] + 1; j < m->b.index[i + 1]; ++j)
                {
                    if (next_group_index(m->b.a[j], top, type, &id))
                    {
                        std::string message("Grouping into residues/molecules is ambiguous");
                        GMX_THROW(gmx::InconsistentInputError(message));
                    }
                }
            }
        }
    }
    // Do a second loop, where things are actually set.
    int id    = -1;
    int group = -1;
    for (int i = 0; i < m->b.nr; ++i)
    {
        const int ii = (type == INDEX_UNKNOWN ? 0 : m->b.a[m->b.index[i]]);
        if (next_group_index(ii, top, type, &id))
        {
            ++group;
        }
        m->mapid[i] = group;
        m->orgid[i] = group;
    }
    // Count also the last group.
    ++group;
    return group;
}

/*!
 * \param[in,out] m    Mapping structure to initialize.
 * \param[in]     b    Block information to use for data.
 *
 * Frees some memory that is not necessary for static index group mappings.
 * Internal pointers are set to point to data in \p b; it is the responsibility
 * of the caller to ensure that the block information matches the contents of
 * the mapping.
 * After this function has been called, the index group provided to
 * gmx_ana_indexmap_update() should always be the same as \p g given here.
 *
 * This function breaks modularity of the index group mapping interface in an
 * ugly way, but allows reducing memory usage of static selections by a
 * significant amount.
 */
void gmx_ana_indexmap_set_static(gmx_ana_indexmap_t* m, t_blocka* b)
{
    sfree(m->mapid);
    sfree(m->mapb.index);
    sfree(m->b.index);
    sfree(m->b.a);
    m->mapb.nalloc_index = 0;
    m->mapb.nalloc_a     = 0;
    m->b.nalloc_index    = 0;
    m->b.nalloc_a        = 0;
    m->mapid             = m->orgid;
    m->mapb.index        = b->index;
    m->mapb.a            = b->a;
    m->b.index           = b->index;
    m->b.a               = b->a;
}

/*!
 * \param[in,out] dest Destination data structure.
 * \param[in]     src  Source mapping.
 * \param[in]     bFirst If true, memory is allocated for \p dest and a full
 *   copy is made; otherwise, only variable parts are copied, and no memory
 *   is allocated.
 *
 * \p dest should have been initialized somehow (calloc() is enough).
 */
void gmx_ana_indexmap_copy(gmx_ana_indexmap_t* dest, gmx_ana_indexmap_t* src, bool bFirst)
{
    if (bFirst)
    {
        gmx_ana_indexmap_reserve(dest, src->b.nr, src->b.nra);
        dest->type  = src->type;
        dest->b.nr  = src->b.nr;
        dest->b.nra = src->b.nra;
        std::memcpy(dest->orgid, src->orgid, dest->b.nr * sizeof(*dest->orgid));
        std::memcpy(dest->b.index, src->b.index, (dest->b.nr + 1) * sizeof(*dest->b.index));
        if (dest->b.nra > 0)
        {
            std::memcpy(dest->b.a, src->b.a, dest->b.nra * sizeof(*dest->b.a));
        }
    }
    dest->mapb.nr  = src->mapb.nr;
    dest->mapb.nra = src->mapb.nra;
    if (src->mapb.nalloc_a > 0)
    {
        if (bFirst)
        {
            snew(dest->mapb.a, src->mapb.nalloc_a);
            dest->mapb.nalloc_a = src->mapb.nalloc_a;
        }
        std::memcpy(dest->mapb.a, src->mapb.a, dest->mapb.nra * sizeof(*dest->mapb.a));
    }
    else
    {
        dest->mapb.a = src->mapb.a;
    }
    std::memcpy(dest->refid, src->refid, dest->mapb.nr * sizeof(*dest->refid));
    std::memcpy(dest->mapid, src->mapid, dest->mapb.nr * sizeof(*dest->mapid));
    std::memcpy(dest->mapb.index, src->mapb.index, (dest->mapb.nr + 1) * sizeof(*dest->mapb.index));
    dest->bStatic = src->bStatic;
}

/*! \brief
 * Helper function to set the source atoms in an index map.
 *
 * \param[in,out] m     Mapping structure.
 * \param[in]     isize Number of atoms in the \p index array.
 * \param[in]     index List of atoms.
 */
static void set_atoms(gmx_ana_indexmap_t* m, int isize, int* index)
{
    m->mapb.nra = isize;
    if (m->mapb.nalloc_a == 0)
    {
        m->mapb.a = index;
    }
    else
    {
        for (int i = 0; i < isize; ++i)
        {
            m->mapb.a[i] = index[i];
        }
    }
}

/*!
 * \param[in,out] m         Mapping structure.
 * \param[in]     g         Current index group.
 * \param[in]     bMaskOnly true if the unused blocks should be masked with
 *   -1 instead of removing them.
 *
 * Updates the index group mapping with the new index group \p g.
 *
 * \see gmx_ana_indexmap_t
 */
void gmx_ana_indexmap_update(gmx_ana_indexmap_t* m, gmx_ana_index_t* g, bool bMaskOnly)
{
    int i, j, bi, bj;

    /* Process the simple cases first */
    if (m->type == INDEX_UNKNOWN && m->b.nra == 0)
    {
        return;
    }
    if (m->type == INDEX_ALL)
    {
        set_atoms(m, g->isize, g->index);
        if (m->b.nr > 0)
        {
            m->mapb.index[1] = g->isize;
        }
        return;
    }
    /* Reset the reference IDs and mapping if necessary */
    const bool bToFull  = (g->isize == m->b.nra);
    const bool bWasFull = (m->mapb.nra == m->b.nra);
    if (bToFull || bMaskOnly)
    {
        if (!m->bStatic)
        {
            for (bj = 0; bj < m->b.nr; ++bj)
            {
                m->refid[bj] = bj;
            }
        }
        if (!bWasFull)
        {
            for (bj = 0; bj < m->b.nr; ++bj)
            {
                m->mapid[bj] = m->orgid[bj];
            }
            for (bj = 0; bj <= m->b.nr; ++bj)
            {
                m->mapb.index[bj] = m->b.index[bj];
            }
        }
        set_atoms(m, m->b.nra, m->b.a);
        m->mapb.nr = m->b.nr;
    }
    /* Exit immediately if the group is static */
    if (bToFull)
    {
        m->bStatic = true;
        return;
    }

    if (bMaskOnly)
    {
        for (i = j = bj = 0; i < g->isize; ++i, ++j)
        {
            /* Find the next atom in the block */
            while (m->b.a[j] != g->index[i])
            {
                ++j;
            }
            /* Mark blocks that did not contain any atoms */
            while (bj < m->b.nr && m->b.index[bj + 1] <= j)
            {
                m->refid[bj++] = -1;
            }
            /* Advance the block index if we have reached the next block */
            if (m->b.index[bj] <= j)
            {
                ++bj;
            }
        }
        /* Mark the last blocks as not accessible */
        while (bj < m->b.nr)
        {
            m->refid[bj++] = -1;
        }
    }
    else
    {
        set_atoms(m, g->isize, g->index);
        for (i = j = bi = 0, bj = -1; i < g->isize; ++i)
        {
            /* Find the next atom in the block */
            while (m->b.a[j] != g->index[i])
            {
                ++j;
            }
            /* If we have reached a new block, add it */
            if (m->b.index[bj + 1] <= j)
            {
                /* Skip any blocks in between */
                while (bj < m->b.nr && m->b.index[bj + 1] <= j)
                {
                    ++bj;
                }
                m->refid[bi]      = bj;
                m->mapid[bi]      = m->orgid[bj];
                m->mapb.index[bi] = i;
                bi++;
            }
        }
        /* Update the number of blocks */
        m->mapb.index[bi] = g->isize;
        m->mapb.nr        = bi;
    }
    m->bStatic = false;
}

/*!
 * \param[in,out] m         Mapping structure to free.
 *
 * All the memory allocated for the mapping structure is freed, and
 * the pointers set to NULL.
 * The pointer \p m is not freed.
 */
void gmx_ana_indexmap_deinit(gmx_ana_indexmap_t* m)
{
    sfree(m->refid);
    if (m->mapid != m->orgid)
    {
        sfree(m->mapid);
    }
    if (m->mapb.nalloc_index > 0)
    {
        sfree(m->mapb.index);
    }
    if (m->mapb.nalloc_a > 0)
    {
        sfree(m->mapb.a);
    }
    sfree(m->orgid);
    if (m->b.nalloc_index > 0)
    {
        sfree(m->b.index);
    }
    if (m->b.nalloc_a > 0)
    {
        sfree(m->b.a);
    }
    gmx_ana_indexmap_clear(m);
}
