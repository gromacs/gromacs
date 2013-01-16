/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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
/*! \internal \file
 * \brief Implementation of functions in indexutil.h.
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <index.h>
#include <smalloc.h>
#include <string2.h>
#include <typedefs.h>
#include <gmx_fatal.h>

#include <indexutil.h>

/********************************************************************
 * gmx_ana_indexgrps_t functions
 ********************************************************************/

/*! \internal \brief
 * Stores a set of index groups.
 */
struct gmx_ana_indexgrps_t
{
    /** Number of index groups. */
    int                 nr;
    /** Array of index groups. */
    gmx_ana_index_t    *g;
};

/*!
 * \param[out] g     Index group structure.
 * \param[in]  ngrps Number of groups for which memory is allocated.
 */
void
gmx_ana_indexgrps_alloc(gmx_ana_indexgrps_t **g, int ngrps)
{
    snew(*g, 1);
    (*g)->nr = ngrps;
    snew((*g)->g,    ngrps);
}

/*!
 * \param[out] g     Index group structure.
 * \param[in]  ngrps Number of index groups.
 * \param[in]  isize Array of index group sizes.
 * \param[in]  index Array of pointers to indices of each group.
 * \param[in]  name  Array of names of the groups.
 * \param[in]  bFree If TRUE, the \p isize, \p index and \p name arrays
 *   are freed after they have been copied.
 */
void
gmx_ana_indexgrps_set(gmx_ana_indexgrps_t **g, int ngrps, int *isize,
                      atom_id **index, char **name, gmx_bool bFree)
{
    int  i;

    gmx_ana_indexgrps_alloc(g, ngrps);
    for (i = 0; i < ngrps; ++i)
    {
        gmx_ana_index_set(&(*g)->g[i], isize[i], index[i], name[i], isize[i]);
    }
    if (bFree)
    {
        sfree(isize);
        sfree(index);
        sfree(name);
    }
}

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
void
gmx_ana_indexgrps_init(gmx_ana_indexgrps_t **g, t_topology *top,
                       const char *fnm)
{
    t_blocka *block = NULL;
    char    **names = NULL;
    int       i, j;

    if (fnm)
    {
        block = init_index(fnm, &names);
    }
    else if (top)
    {
        block = new_blocka();
        analyse(&top->atoms, block, &names, FALSE, FALSE);
    }
    else
    {
        snew(*g, 1);
        (*g)->nr    = 0;
        (*g)->g     = NULL;
        return;
    }

    gmx_ana_indexgrps_alloc(g, block->nr);
    for (i = 0; i < block->nr; ++i)
    {
        gmx_ana_index_t *grp = &(*g)->g[i];

        grp->isize = block->index[i+1] - block->index[i];
        snew(grp->index, grp->isize);
        for (j = 0; j < grp->isize; ++j)
        {
            grp->index[j] = block->a[block->index[i]+j];
        }
        grp->name         = names[i];
        grp->nalloc_index = grp->isize;
    }

    done_blocka(block);
    sfree(block);
    sfree(names);
}

/*!
 * \param[out] g     Index group structure.
 * \param[in]  top   Topology structure.
 * \param[in]  fnm   File name for the index file.
 * \param[in]  ngrps Number of required groups.
 *   Memory is automatically allocated.
 *
 * One of \p top or \p fnm can be NULL, but not both.
 * If \p top is NULL, an index file is required and the groups are read
 * from the file (uses Gromacs routine rd_index()).
 * If \p fnm is NULL, default groups are constructed based on the
 * topology (uses Gromacs routine get_index()).
 */
void
gmx_ana_indexgrps_get(gmx_ana_indexgrps_t **g, t_topology *top,
                      const char *fnm, int ngrps)
{
    int      *isize;
    atom_id **index;
    char    **name;

    snew(isize, ngrps);
    snew(index, ngrps);
    snew(name,  ngrps);
    if (!top)
    {
        rd_index(fnm, ngrps, isize, index, name);
    }
    else
    {
        get_index(&(top->atoms), fnm, ngrps, isize, index, name);
    }
    gmx_ana_indexgrps_set(g, ngrps, isize, index, name, TRUE);
}

/*!
 * \param[out] g     Index group structure.
 * \param[in]  fnm   File name for the index file.
 * \param[in]  ngrps Number of required groups.
 *   Memory is automatically allocated.
 *
 * This is a convenience function for calling the Gromacs routine
 * rd_index().
 */
void
gmx_ana_indexgrps_rd(gmx_ana_indexgrps_t **g, const char *fnm, int ngrps)
{
    gmx_ana_indexgrps_get(g, NULL, fnm, ngrps);
}

/*!
 * \param[in] g  Index groups structure.
 *
 * The pointer \p g is invalid after the call.
 */
void
gmx_ana_indexgrps_free(gmx_ana_indexgrps_t *g)
{
    int  i;

    if (g->nr == 0)
    {
        sfree(g);
        return;
    }
    for (i = 0; i < g->nr; ++i)
    {
        gmx_ana_index_deinit(&g->g[i]);
    }
    sfree(g->g);
    g->nr    = 0;
    g->g     = NULL;
    sfree(g);
}

/*!
 * \param[out] dest Destination index groups.
 * \param[in]  src  Source index groups.
 *
 * A deep copy is made for all fields, including the group names.
 */
void
gmx_ana_indexgrps_clone(gmx_ana_indexgrps_t **dest, gmx_ana_indexgrps_t *src)
{
    int g;

    gmx_ana_indexgrps_alloc(dest, src->nr);
    for (g = 0; g < src->nr; ++g)
    {
        gmx_ana_index_copy(&(*dest)->g[g], &src->g[g], TRUE);
    }
}

/*!
 * \param[out] g     Index group structure.
 * \returns    TRUE if \p g is empty, i.e., has 0 index groups.
 */
gmx_bool
gmx_ana_indexgrps_is_empty(gmx_ana_indexgrps_t *g)
{
    return g->nr == 0;
}

/*!
 * \param[in]  g     Index groups structure.
 * \param[in]  n     Index group number to get.
 * \returns    Pointer to the \p n'th index group in \p g.
 *
 * The returned pointer should not be freed.
 */
gmx_ana_index_t *
gmx_ana_indexgrps_get_grp(gmx_ana_indexgrps_t *g, int n)
{
    if (n < 0 || n >= g->nr)
    {
        return NULL;
    }
    return &g->g[n];
}

/*!
 * \param[out] dest Output structure.
 * \param[in]  src  Input index groups.
 * \param[in]  n    Number of the group to extract.
 * \returns TRUE if \p n is a valid group in \p src, FALSE otherwise.
 */
gmx_bool
gmx_ana_indexgrps_extract(gmx_ana_index_t *dest, gmx_ana_indexgrps_t *src, int n)
{
    if (n < 0 || n >= src->nr)
    {
        dest->isize = 0;
        return FALSE;
    }

    gmx_ana_index_copy(dest, &src->g[n], TRUE);
    return TRUE;
}

/*!
 * \param[out] dest Output structure.
 * \param[in]  src  Input index groups.
 * \param[in]  name Name (or part of the name) of the group to extract.
 * \returns TRUE if \p name is a valid group in \p src, FALSE otherwise.
 *
 * Uses the Gromacs routine find_group() to find the actual group;
 * the comparison is case-insensitive.
 */
gmx_bool
gmx_ana_indexgrps_find(gmx_ana_index_t *dest, gmx_ana_indexgrps_t *src, char *name)
{
    int    i;
    char **names;

    snew(names, src->nr);
    for (i = 0; i < src->nr; ++i)
    {
        names[i] = src->g[i].name;
    }
    i = find_group(name, src->nr, names);
    sfree(names);
    if (i == NOTSET)
    {
        dest->isize = 0;
        return FALSE;
    }

    return gmx_ana_indexgrps_extract(dest, src, i);
}

/*!
 * \param[in]  g      Index groups to print.
 * \param[in]  maxn   Maximum number of indices to print
 *      (-1 = print all, 0 = print only names).
 */
void
gmx_ana_indexgrps_print(gmx_ana_indexgrps_t *g, int maxn)
{
    int  i;

    for (i = 0; i < g->nr; ++i)
    {
        fprintf(stderr, " %2d: ", i);
        gmx_ana_index_dump(&g->g[i], i, maxn);
    }
}

/********************************************************************
 * gmx_ana_index_t functions
 ********************************************************************/

/*!
 * \param[in,out] g      Index group structure.
 * \param[in]     isize  Maximum number of atoms to reserve space for.
 */
void
gmx_ana_index_reserve(gmx_ana_index_t *g, int isize)
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
void
gmx_ana_index_squeeze(gmx_ana_index_t *g)
{
    srenew(g->index, g->isize);
    g->nalloc_index = g->isize;
}

/*!
 * \param[out] g      Output structure.
 *
 * Any contents of \p g are discarded without freeing.
 */
void
gmx_ana_index_clear(gmx_ana_index_t *g)
{
    g->isize        = 0;
    g->index        = NULL;
    g->name         = NULL;
    g->nalloc_index = 0;
}

/*!
 * \param[out] g      Output structure.
 * \param[in]  isize  Number of atoms in the new group.
 * \param[in]  index  Array of \p isize atoms (can be NULL if \p isize is 0).
 * \param[in]  name   Name for the new group (can be NULL).
 * \param[in]  nalloc Number of elements allocated for \p index
 *   (if 0, \p index is not freed in gmx_ana_index_deinit())
 *
 * No copy if \p index is made.
 */
void
gmx_ana_index_set(gmx_ana_index_t *g, int isize, atom_id *index, char *name,
                  int nalloc)
{
    g->isize        = isize;
    g->index        = index;
    g->name         = name;
    g->nalloc_index = nalloc;
}

/*!
 * \param[out] g      Output structure.
 * \param[in]  natoms Number of atoms.
 * \param[in]  name   Name for the new group (can be NULL).
 */
void
gmx_ana_index_init_simple(gmx_ana_index_t *g, int natoms, char *name)
{
    int  i;

    g->isize = natoms;
    snew(g->index, natoms);
    for (i = 0; i < natoms; ++i)
    {
        g->index[i] = i;
    }
    g->name         = name;
    g->nalloc_index = natoms;
}

/*!
 * \param[in] g  Index group structure.
 *
 * The pointer \p g is not freed.
 */
void
gmx_ana_index_deinit(gmx_ana_index_t *g)
{
    if (g->nalloc_index > 0)
    {
        sfree(g->index);
    }
    sfree(g->name);
    gmx_ana_index_clear(g);
}

/*!
 * \param[out] dest   Destination index group.
 * \param[in]  src    Source index group.
 * \param[in]  bAlloc If TRUE, memory is allocated at \p dest; otherwise,
 *   it is assumed that enough memory has been allocated for index.
 *
 * A deep copy of the name is only made if \p bAlloc is TRUE.
 */
void
gmx_ana_index_copy(gmx_ana_index_t *dest, gmx_ana_index_t *src, gmx_bool bAlloc)
{
    dest->isize = src->isize;
    if (dest->isize > 0)
    {
        if (bAlloc)
        {
            snew(dest->index, dest->isize);
            dest->nalloc_index = dest->isize;
        }
        memcpy(dest->index, src->index, dest->isize*sizeof(*dest->index));
    }
    if (bAlloc && src->name)
    {
        dest->name = strdup(src->name);
    }
    else if (bAlloc || src->name)
    {
        dest->name = src->name;
    }
}

/*!
 * \param[in]  g      Index group to print.
 * \param[in]  i      Group number to use if the name is NULL.
 * \param[in]  maxn   Maximum number of indices to print (-1 = print all).
 */
void
gmx_ana_index_dump(gmx_ana_index_t *g, int i, int maxn)
{
    int  j, n;

    if (g->name)
    {
        fprintf(stderr, "\"%s\"", g->name);
    }
    else
    {
        fprintf(stderr, "Group %d", i+1);
    }
    fprintf(stderr, " (%d atoms)", g->isize);
    if (maxn != 0)
    {
        fprintf(stderr, ":");
        n = g->isize;
        if (maxn >= 0 && n > maxn)
        {
            n = maxn;
        }
        for (j = 0; j < n; ++j)
        {
            fprintf(stderr, " %d", g->index[j]+1);
        }
        if (n < g->isize)
        {
            fprintf(stderr, " ...");
        }
    }
    fprintf(stderr, "\n");
}

/*!
 * \param[in]  g      Input index group.
 * \param[in]  natoms Number of atoms to check against.
 *
 * If any atom index in the index group is less than zero or >= \p natoms,
 * gmx_fatal() is called.
 */
void
gmx_ana_index_check(gmx_ana_index_t *g, int natoms)
{
    int  j;

    for (j = 0; j < g->isize; ++j)
    {
        if (g->index[j] >= natoms)
        {
            gmx_fatal(FARGS, "Atom index (%d) in index group %s (%d atoms) "
                      "larger than number of atoms in trajectory (%d atoms)",
                      g->index[j], g->name, g->isize, natoms);
        }
        else if (g->index[j] < 0)
        {
            gmx_fatal(FARGS, "Atom index (%d) in index group %s (%d atoms) "
                      "is less than zero",
                      g->index[j], g->name, g->isize);
        }
    }
}

/*!
 * \param[in]  g      Index group to check.
 * \returns    TRUE if the index group is sorted and has no duplicates,
 *   FALSE otherwise.
 */
gmx_bool
gmx_ana_index_check_sorted(gmx_ana_index_t *g)
{
    int  i;

    for (i = 0; i < g->isize-1; ++i)
    {
        if (g->index[i+1] <= g->index[i])
        {
            return FALSE;
        }
    }
    return TRUE;
}

/********************************************************************
 * Set operations
 ********************************************************************/

/** Helper function for gmx_ana_index_sort(). */
static int
cmp_atomid(const void *a, const void *b)
{
    if (*(atom_id *)a < *(atom_id *)b)
    {
        return -1;
    }
    if (*(atom_id *)a > *(atom_id *)b)
    {
        return 1;
    }
    return 0;
}

/*!
 * \param[in,out] g  Index group to be sorted.
 */
void
gmx_ana_index_sort(gmx_ana_index_t *g)
{
    qsort(g->index, g->isize, sizeof(*g->index), cmp_atomid);
}

/*!
 * \param[in]  a      Index group to check.
 * \param[in]  b      Index group to check.
 * \returns    TRUE if \p a and \p b are equal, FALSE otherwise.
 */
gmx_bool
gmx_ana_index_equals(gmx_ana_index_t *a, gmx_ana_index_t *b)
{
    int  i;

    if (a->isize != b->isize)
    {
        return FALSE;
    }
    for (i = 0; i < a->isize; ++i)
    {
        if (a->index[i] != b->index[i])
        {
            return FALSE;
        }
    }
    return TRUE;
}

/*!
 * \param[in]  a      Index group to check against.
 * \param[in]  b      Index group to check.
 * \returns    TRUE if \p b is contained in \p a,
 *   FALSE otherwise.
 *
 * If the elements are not in the same order in both groups, the function
 * fails. However, the groups do not need to be sorted.
 */
gmx_bool
gmx_ana_index_contains(gmx_ana_index_t *a, gmx_ana_index_t *b)
{
    int  i, j;

    for (i = j = 0; j < b->isize; ++i, ++j)
    {
        while (i < a->isize && a->index[i] != b->index[j])
        {
            ++i;
        }
        if (i == a->isize)
        {
            return FALSE;
        }
    }
    return TRUE;
}

/*!
 * \param[out] dest Output index group (the intersection of \p a and \p b).
 * \param[in]  a    First index group.
 * \param[in]  b    Second index group.
 *
 * \p dest can be the same as \p a or \p b.
 */
void
gmx_ana_index_intersection(gmx_ana_index_t *dest,
                           gmx_ana_index_t *a, gmx_ana_index_t *b)
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
void
gmx_ana_index_difference(gmx_ana_index_t *dest,
                         gmx_ana_index_t *a, gmx_ana_index_t *b)
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
int
gmx_ana_index_difference_size(gmx_ana_index_t *a, gmx_ana_index_t *b)
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
void
gmx_ana_index_partition(gmx_ana_index_t *dest1, gmx_ana_index_t *dest2,
                        gmx_ana_index_t *src, gmx_ana_index_t *g)

{
    int i, j, k;

    dest2->index = dest1->index + g->isize;
    dest2->isize = src->isize - g->isize;
    for (i = g->isize-1, j = src->isize-1, k = dest2->isize-1; i >= 0; --i, --j)
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
    gmx_ana_index_copy(dest1, g, FALSE);
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
void
gmx_ana_index_union(gmx_ana_index_t *dest,
                    gmx_ana_index_t *a, gmx_ana_index_t *b)
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
void
gmx_ana_index_merge(gmx_ana_index_t *dest,
                    gmx_ana_index_t *a, gmx_ana_index_t *b)
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

/*!
 * \param[in,out] t    Output block.
 * \param[in]  top  Topology structure
 *   (only used if \p type is \ref INDEX_RES or \ref INDEX_MOL, can be NULL
 *   otherwise).
 * \param[in]  g    Index group
 *   (can be NULL if \p type is \ref INDEX_UNKNOWN).
 * \param[in]  type Type of partitioning to make.
 * \param[in]  bComplete
 *   If TRUE, the index group is expanded to include any residue/molecule
 *   (depending on \p type) that is partially contained in the group.
 *   If \p type is not INDEX_RES or INDEX_MOL, this has no effect.
 *
 * \p m should have been initialized somehow (calloc() is enough) unless
 * \p type is INDEX_UNKNOWN.
 * \p g should be sorted.
 */
void
gmx_ana_index_make_block(t_blocka *t, t_topology *top, gmx_ana_index_t *g,
                         e_index_t type, gmx_bool bComplete)
{
    int      i, j, ai;
    int      id, cur;

    if (type == INDEX_UNKNOWN)
    {
        t->nr           = 1;
        snew(t->index, 2);
        t->nalloc_index = 2;
        t->index[0]     = 0;
        t->index[1]     = 0;
        t->nra          = 0;
        t->a            = NULL;
        t->nalloc_a     = 0;
        return;
    }

    /* bComplete only does something for INDEX_RES or INDEX_MOL, so turn it
     * off otherwise. */
    if (type != INDEX_RES && type != INDEX_MOL)
    {
        bComplete = FALSE;
    }
    /* Allocate memory for the atom array and fill it unless we are using
     * completion. */
    if (bComplete)
    {
        t->nra = 0;
        /* We may allocate some extra memory here because we don't know in
         * advance how much will be needed. */
        if (t->nalloc_a < top->atoms.nr)
        {
            srenew(t->a, top->atoms.nr);
            t->nalloc_a = top->atoms.nr;
        }
    }
    else
    {
        t->nra      = g->isize;
        if (t->nalloc_a < g->isize)
        {
            srenew(t->a, g->isize);
            t->nalloc_a = g->isize;
        }
        memcpy(t->a, g->index, g->isize*sizeof(*(t->a)));
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
    t->nr = 0;
    j     = 0; /* j is used by residue completion for the first atom not stored */
    id    = cur = -1;
    for (i = 0; i < g->isize; ++i)
    {
        ai = g->index[i];
        /* Find the ID number of the atom/residue/molecule corresponding to
         * atom ai. */
        switch (type)
        {
            case INDEX_ATOM:
                id = ai;
                break;
            case INDEX_RES:
                id = top->atoms.atom[ai].resind;
                break;
            case INDEX_MOL:
                while (ai >= top->mols.index[id+1])
                {
                    id++;
                }
                break;
            case INDEX_UNKNOWN: /* Should not occur */
            case INDEX_ALL:
                id = 0;
                break;
        }
        /* If this is the first atom in a new block, initialize the block. */
        if (id != cur)
        {
            if (bComplete)
            {
                /* For completion, we first set the start of the block. */
                t->index[t->nr++] = t->nra;
                /* And then we find all the atoms that should be included. */
                switch (type)
                {
                    case INDEX_RES:
                        while (top->atoms.atom[j].resind != id)
                        {
                            ++j;
                        }
                        while (j < top->atoms.nr && top->atoms.atom[j].resind == id)
                        {
                            t->a[t->nra++] = j;
                            ++j;
                        }
                        break;

                    case INDEX_MOL:
                        for (j = top->mols.index[id]; j < top->mols.index[id+1]; ++j)
                        {
                            t->a[t->nra++] = j;
                        }
                        break;

                    default: /* Should not be reached */
                        gmx_bug("internal error");
                        break;
                }
            }
            else
            {
                /* If not using completion, simply store the start of the block. */
                t->index[t->nr++] = i;
            }
            cur = id;
        }
    }
    /* Set the end of the last block */
    t->index[t->nr] = t->nra;
    /* Free any unnecessary memory */
    srenew(t->index, t->nr+1);
    t->nalloc_index = t->nr+1;
    if (bComplete)
    {
        srenew(t->a, t->nra);
        t->nalloc_a = t->nra;
    }
}

/*!
 * \param[in] g   Index group to check.
 * \param[in] b   Block data to check against.
 * \returns   TRUE if \p g consists of one or more complete blocks from \p b,
 *   FALSE otherwise.
 *
 * The atoms in \p g are assumed to be sorted.
 */
gmx_bool
gmx_ana_index_has_full_blocks(gmx_ana_index_t *g, t_block *b)
{
    int  i, j, bi;

    i = bi = 0;
    /* Each round in the loop matches one block */
    while (i < g->isize)
    {
        /* Find the block that begins with the first unmatched atom */
        while (bi < b->nr && b->index[bi] != g->index[i])
        {
            ++bi;
        }
        /* If not found, or if too large, return */
        if (bi == b->nr || i + b->index[bi+1] -  b->index[bi] > g->isize)
        {
            return FALSE;
        }
        /* Check that the block matches the index */
        for (j = b->index[bi]; j < b->index[bi+1]; ++j, ++i)
        {
            if (g->index[i] != j)
            {
                return FALSE;
            }
        }
        /* Move the search to the next block */
        ++bi;
    }
    return TRUE;
}

/*!
 * \param[in] g   Index group to check.
 * \param[in] b   Block data to check against.
 * \returns   TRUE if \p g consists of one or more complete blocks from \p b,
 *   FALSE otherwise.
 *
 * The atoms in \p g and \p b->a are assumed to be in the same order.
 */
gmx_bool
gmx_ana_index_has_full_ablocks(gmx_ana_index_t *g, t_blocka *b)
{
    int  i, j, bi;

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
        if (bi == b->nr || i + b->index[bi+1] -  b->index[bi] > g->isize)
        {
            return FALSE;
        }
        /* Check that the block matches the index */
        for (j = b->index[bi]; j < b->index[bi+1]; ++j, ++i)
        {
            if (b->a[j] != g->index[i])
            {
                return FALSE;
            }
        }
        /* Move the search to the next block */
        ++bi;
    }
    return TRUE;
}

/*!
 * \param[in] g     Index group to check.
 * \param[in] type  Block data to check against.
 * \param[in] top   Topology data.
 * \returns   TRUE if \p g consists of one or more complete elements of type
 *   \p type, FALSE otherwise.
 *
 * If \p type is \ref INDEX_ATOM, the return value is always TRUE.
 * If \p type is \ref INDEX_UNKNOWN or \ref INDEX_ALL, the return value is
 * always FALSE.
 */
gmx_bool
gmx_ana_index_has_complete_elems(gmx_ana_index_t *g, e_index_t type,
                                 t_topology *top)
{
    switch (type)
    {
        case INDEX_UNKNOWN:
        case INDEX_ALL:
            return FALSE;

        case INDEX_ATOM:
            return TRUE;

        case INDEX_RES:
        {
            int      i, ai;
            int      id, prev;

            prev = -1;
            for (i = 0; i < g->isize; ++i)
            {
                ai = g->index[i];
                id = top->atoms.atom[ai].resind;
                if (id != prev)
                {
                    if (ai > 0 && top->atoms.atom[ai-1].resind == id)
                    {
                        return FALSE;
                    }
                    if (i > 0 && g->index[i-1] < top->atoms.nr - 1
                        && top->atoms.atom[g->index[i-1]+1].resind == prev)
                    {
                        return FALSE;
                    }
                }
                prev = id;
            }
            if (g->index[i-1] < top->atoms.nr - 1
                && top->atoms.atom[g->index[i-1]+1].resind == prev)
            {
                return FALSE;
            }
            break;
        }

        case INDEX_MOL:
            return gmx_ana_index_has_full_blocks(g, &top->mols);
    }
    return TRUE;
}

/*!
 * \param[out] m      Output structure.
 *
 * Any contents of \p m are discarded without freeing.
 */
void
gmx_ana_indexmap_clear(gmx_ana_indexmap_t *m)
{
    m->type              = INDEX_UNKNOWN;
    m->nr                = 0;
    m->refid             = NULL;
    m->mapid             = NULL;
    m->mapb.nr           = 0;
    m->mapb.index        = NULL;
    m->mapb.nalloc_index = 0;
    m->orgid             = NULL;
    m->b.nr              = 0;
    m->b.index           = NULL;
    m->b.nra             = 0;
    m->b.a               = NULL;
    m->b.nalloc_index    = 0;
    m->b.nalloc_a        = 0;
    m->bStatic           = TRUE;
    m->bMapStatic        = TRUE;
}

/*!
 * \param[in,out] m      Mapping structure.
 * \param[in]     nr     Maximum number of blocks to reserve space for.
 * \param[in]     isize  Maximum number of atoms to reserve space for.
 */
void
gmx_ana_indexmap_reserve(gmx_ana_indexmap_t *m, int nr, int isize)
{
    if (m->mapb.nalloc_index < nr + 1)
    {
        srenew(m->refid,      nr);
        srenew(m->mapid,      nr);
        srenew(m->orgid,      nr);
        srenew(m->mapb.index, nr + 1);
        srenew(m->b.index,    nr + 1);
        m->mapb.nalloc_index = nr + 1;
        m->b.nalloc_index    = nr + 1;
    }
    if (m->b.nalloc_a < isize)
    {
        srenew(m->b.a,        isize);
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
void
gmx_ana_indexmap_init(gmx_ana_indexmap_t *m, gmx_ana_index_t *g,
                      t_topology *top, e_index_t type)
{
    int      i, ii, mi;

    m->type   = type;
    gmx_ana_index_make_block(&m->b, top, g, type, FALSE);
    gmx_ana_indexmap_reserve(m, m->b.nr, m->b.nra);
    m->nr = m->b.nr;
    for (i = mi = 0; i < m->nr; ++i)
    {
        ii = (type == INDEX_UNKNOWN ? 0 : m->b.a[m->b.index[i]]);
        switch (type)
        {
            case INDEX_ATOM:
                m->orgid[i] = ii;
                break;
            case INDEX_RES:
                m->orgid[i] = top->atoms.atom[ii].resind;
                break;
            case INDEX_MOL:
                while (top->mols.index[mi+1] <= ii)
                {
                    ++mi;
                }
                m->orgid[i] = mi;
                break;
            case INDEX_ALL:
            case INDEX_UNKNOWN:
                m->orgid[i] = 0;
                break;
        }
    }
    for (i = 0; i < m->nr; ++i)
    {
        m->refid[i] = i;
        m->mapid[i] = m->orgid[i];
    }
    m->mapb.nr = m->nr;
    memcpy(m->mapb.index, m->b.index, (m->nr+1)*sizeof(*(m->mapb.index)));
    m->bStatic    = TRUE;
    m->bMapStatic = TRUE;
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
void
gmx_ana_indexmap_set_static(gmx_ana_indexmap_t *m, t_blocka *b)
{
    sfree(m->mapid);
    m->mapid = m->orgid;
    sfree(m->b.index);
    m->b.nalloc_index = 0;
    m->b.index        = b->index;
    sfree(m->mapb.index);
    m->mapb.nalloc_index = 0;
    m->mapb.index        = m->b.index;
    sfree(m->b.a);
    m->b.nalloc_a = 0;
    m->b.a        = b->a;
}

/*!
 * \param[in,out] dest Destination data structure.
 * \param[in]     src  Source mapping.
 * \param[in]     bFirst If TRUE, memory is allocated for \p dest and a full
 *   copy is made; otherwise, only variable parts are copied, and no memory
 *   is allocated.
 *
 * \p dest should have been initialized somehow (calloc() is enough).
 */
void
gmx_ana_indexmap_copy(gmx_ana_indexmap_t *dest, gmx_ana_indexmap_t *src, gmx_bool bFirst)
{
    if (bFirst)
    {
        gmx_ana_indexmap_reserve(dest, src->b.nr, src->b.nra);
        dest->type       = src->type;
        dest->b.nr       = src->b.nr;
        dest->b.nra      = src->b.nra;
        memcpy(dest->orgid,      src->orgid,      dest->b.nr*sizeof(*dest->orgid));
        memcpy(dest->b.index,    src->b.index,   (dest->b.nr+1)*sizeof(*dest->b.index));
        memcpy(dest->b.a,        src->b.a,        dest->b.nra*sizeof(*dest->b.a));
    }
    dest->nr         = src->nr;
    dest->mapb.nr    = src->mapb.nr;
    memcpy(dest->refid,      src->refid,      dest->nr*sizeof(*dest->refid));
    memcpy(dest->mapid,      src->mapid,      dest->nr*sizeof(*dest->mapid));
    memcpy(dest->mapb.index, src->mapb.index, (dest->mapb.nr+1)*sizeof(*dest->mapb.index));
    dest->bStatic    = src->bStatic;
    dest->bMapStatic = src->bMapStatic;
}

/*!
 * \param[in,out] m         Mapping structure.
 * \param[in]     g         Current index group.
 * \param[in]     bMaskOnly TRUE if the unused blocks should be masked with
 *   -1 instead of removing them.
 *
 * Updates the index group mapping with the new index group \p g.
 *
 * \see gmx_ana_indexmap_t
 */
void
gmx_ana_indexmap_update(gmx_ana_indexmap_t *m, gmx_ana_index_t *g,
                        gmx_bool bMaskOnly)
{
    int      i, j, bi, bj;
    gmx_bool bStatic;

    /* Process the simple cases first */
    if (m->type == INDEX_UNKNOWN && m->b.nra == 0)
    {
        return;
    }
    if (m->type == INDEX_ALL)
    {
        if (m->b.nr > 0)
        {
            m->mapb.index[1] = g->isize;
        }
        return;
    }
    /* Reset the reference IDs and mapping if necessary */
    bStatic = (g->isize == m->b.nra && m->nr == m->b.nr);
    if (bStatic || bMaskOnly)
    {
        if (!m->bStatic)
        {
            for (bj = 0; bj < m->b.nr; ++bj)
            {
                m->refid[bj] = bj;
            }
        }
        if (!m->bMapStatic)
        {
            for (bj = 0; bj < m->b.nr; ++bj)
            {
                m->mapid[bj] = m->orgid[bj];
            }
            for (bj = 0; bj <= m->b.nr; ++bj)
            {
                m->mapb.index[bj] = m->b.index[bj];
            }
            m->bMapStatic = TRUE;
        }
    }
    /* Exit immediately if the group is static */
    if (bStatic)
    {
        m->bStatic = TRUE;
        return;
    }

    if (bMaskOnly)
    {
        m->nr = m->b.nr;
        for (i = j = bj = 0; i < g->isize; ++i, ++j)
        {
            /* Find the next atom in the block */
            while (m->b.a[j] != g->index[i])
            {
                ++j;
            }
            /* Mark blocks that did not contain any atoms */
            while (bj < m->b.nr && m->b.index[bj+1] <= j)
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
        for (i = j = bi = 0, bj = -1; i < g->isize; ++i)
        {
            /* Find the next atom in the block */
            while (m->b.a[j] != g->index[i])
            {
                ++j;
            }
            /* If we have reached a new block, add it */
            if (m->b.index[bj+1] <= j)
            {
                /* Skip any blocks in between */
                while (bj < m->b.nr && m->b.index[bj+1] <= j)
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
        m->nr             = bi;
        m->bMapStatic     = FALSE;
    }
    m->mapb.nr = m->nr;
    m->bStatic = FALSE;
}

/*!
 * \param[in,out] m         Mapping structure to free.
 *
 * All the memory allocated for the mapping structure is freed, and
 * the pointers set to NULL.
 * The pointer \p m is not freed.
 */
void
gmx_ana_indexmap_deinit(gmx_ana_indexmap_t *m)
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
