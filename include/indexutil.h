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
/*! \file
 * \brief API for handling index files and index groups.
 *
 * The API contains functions and data structures for handling index
 * files more conveniently than as several separate variables.
 * In addition to basic functions for initializing the data structures and
 * making copies, functions are provided for performing (most) set operations
 * on sorted index groups.
 * There is also a function for partitioning a index group based on
 * topology information such as residues or molecules.
 * Finally, there is a set of functions for constructing mappings between
 * an index group and its subgroups such.
 * These can be used with dynamic index group in calculations if one
 * needs to have a unique ID for each possible atom/residue/molecule in the
 * selection, e.g., for analysis of dynamics or for look-up tables.
 *
 * Mostly, these functions are used internally by the library and the
 * selection engine.
 * However, some of the checking functions can be useful in user code to
 * check the validity of input groups.
 * Also, the mapping functions are useful when dealing with dynamic index
 * groups.
 */
#ifndef INDEXUTIL_H
#define INDEXUTIL_H
#include "visibility.h"
#include "typedefs.h"

#ifdef __cplusplus
extern "C"
{
#endif

/** Stores a set of index groups. */
typedef struct gmx_ana_indexgrps_t gmx_ana_indexgrps_t;

/*! \brief
 * Specifies the type of index partition or index mapping in several contexts.
 *
 * \see gmx_ana_index_make_block(), gmx_ana_indexmap_init()
 */
typedef enum
{
    INDEX_UNKNOWN, /**< Unknown index type.*/
    INDEX_ATOM,    /**< Each atom in a separate block.*/
    INDEX_RES,     /**< Each residue in a separate block.*/
    INDEX_MOL,     /**< Each molecule in a separate block.*/
    INDEX_ALL      /**< All atoms in a single block.*/
} e_index_t;

/*! \brief
 * Stores a single index group.
 */
typedef struct gmx_ana_index_t
{
    /** Number of atoms. */
    int                 isize;
    /** List of atoms. */
    atom_id            *index;
    /** Group name. */
    char               *name;
    /** Number of items allocated for \p index. */
    int                 nalloc_index;
} gmx_ana_index_t;

/*! \brief
 * Data structure for calculating index group mappings.
 */
typedef struct gmx_ana_indexmap_t
{
    /** Type of the mapping. */
    e_index_t           type;
    /*! \brief
     * Current number of mapped values.
     *
     * This is the current number of values in the \p refid and \p mapid
     * arrays.
     * If \p bMaskOnly is provided to gmx_ana_indexmap_update(), this
     * is always equal to \p b.nr, i.e., the number of blocks in the
     * original index group.
     */
    int                 nr;
    /*! \brief
     * Current reference IDs.
     *
     * This array provides a mapping from the current index group (last given
     * to gmx_ana_indexmap_update()) to the blocks in \p b, i.e., the
     * original index group used in gmx_ana_indexmap_init().
     * The mapping is zero-based.
     * If \p bMaskOnly is provided to gmx_ana_indexmap_update(), the indices
     * for blocks not present in the current group are set to -1, otherwise
     * they are removed completely and the \p nr field updated.
     */
    int                *refid;
    /*! \brief
     * Current mapped IDs.
     *
     * This array provides an arbitrary mapping from the current index group
     * to the original index group. Instead of a zero-based mapping, the
     * values from the \p orgid array are used. That is,
     * \c mapid[i]=orgid[refid[i]].
     * If \p bMaskOnly is provided to gmx_ana_indexmap_update(), this array
     * equals \p orgid.
     */
    int                *mapid;
    /*! \brief
     * Mapped block structure.
     *
     * A block structure that corresponds to the current index group.
     */
    t_block             mapb;

    /*! \brief
     * Arbitrary ID numbers for the blocks.
     *
     * This array has \p b.nr elements, each defining an ID number for a
     * block in \p b.
     * These are initialized in gmx_ana_indexmap_init() based on the type:
     *  - \ref INDEX_ATOM : the atom indices
     *  - \ref INDEX_RES :  the residue numbers
     *  - \ref INDEX_MOL :  the molecule numbers
     *
     * All the above numbers are zero-based.
     * After gmx_ana_indexmap_init(), the user is free to change these values
     * if the above are not appropriate.
     * The mapped values can be read through \p mapid.
     */
    int                *orgid;

    /*! \brief
     * Block data that defines the mapping (internal use only).
     *
     * The data is initialized by gmx_ana_indexmap_init() and is not changed
     * after that.
     * Hence, it cannot be directly applied to the index group passed to
     * gmx_ana_indexmap_update() unless \p bMaskOnly was specified or the
     * index group is identical to the one provided to gmx_ana_indexmap_init().
     */
    t_blocka            b;
    /*! \brief
     * TRUE if the current reference IDs are for the whole group (internal use only).
     *
     * This is used internally to optimize the evaluation such that
     * gmx_ana_indexmap_update() does not take any time if the group is
     * actually static.
     */
    gmx_bool                bStatic;
    /*! \brief
     * TRUE if the current mapping is for the whole group (internal use only).
     *
     * This is used internally to optimize the evaluation such that
     * gmx_ana_indexmap_update() does not take any time if the group is
     * actually static.
     */
    gmx_bool                bMapStatic;
} gmx_ana_indexmap_t;


/*! \name Functions for handling gmx_ana_indexgrps_t
 */
/*@{*/
/** Allocate memory for index groups. */
void
gmx_ana_indexgrps_alloc(gmx_ana_indexgrps_t **g, int ngrps);
/** Initializes index groups from arrays. */
void
gmx_ana_indexgrps_set(gmx_ana_indexgrps_t **g, int ngrps, int *isize,
                      atom_id **index, char **name, gmx_bool bFree);
/** Reads index groups from a file or constructs them from topology. */
void
gmx_ana_indexgrps_init(gmx_ana_indexgrps_t **g, t_topology *top,
                       const char *fnm);
/** Ask user to select index groups, possibly constructing groups from
 *  topology. */
void
gmx_ana_indexgrps_get(gmx_ana_indexgrps_t **g, t_topology *top,
                      const char *fnm, int ngrps);
/** Ask user to select index groups from those specified in a file. */
void
gmx_ana_indexgrps_rd(gmx_ana_indexgrps_t **g, const char *fnm, int ngrps);
/** Frees memory allocated for index groups. */
void
gmx_ana_indexgrps_free(gmx_ana_indexgrps_t *g);
/** Create a deep copy of \c gmx_ana_indexgrps_t. */
void
gmx_ana_indexgrps_clone(gmx_ana_indexgrps_t **dest, gmx_ana_indexgrps_t *src);
/** Returns TRUE if the index group structure is emtpy. */
gmx_bool
gmx_ana_indexgrps_is_empty(gmx_ana_indexgrps_t *g);

/** Returns a pointer to an index group. */
gmx_ana_index_t *
gmx_ana_indexgrps_get_grp(gmx_ana_indexgrps_t *g, int n);
/** Extracts a single index group. */
gmx_bool
gmx_ana_indexgrps_extract(gmx_ana_index_t *dest, gmx_ana_indexgrps_t *src, int n);
/** Finds and extracts a single index group by name. */
gmx_bool
gmx_ana_indexgrps_find(gmx_ana_index_t *dest, gmx_ana_indexgrps_t *src, char *name);

/** Writes out a list of index groups. */
void
gmx_ana_indexgrps_print(gmx_ana_indexgrps_t *g, int maxn);
/*@}*/

/*! \name Functions for handling gmx_ana_index_t
 */
/*@{*/
/** Reserves memory to store an index group of size \p isize. */
void
gmx_ana_index_reserve(gmx_ana_index_t *g, int isize);
/** Frees any memory not necessary to hold the current contents. */
void
gmx_ana_index_squeeze(gmx_ana_index_t *g);
/** Initializes an empty index group. */
void
gmx_ana_index_clear(gmx_ana_index_t *g);
/** Constructs a \c gmx_ana_index_t from given values. */
void
gmx_ana_index_set(gmx_ana_index_t *g, int isize, atom_id *index, char *name,
                  int nalloc);
/** Creates a simple index group from the first to the \p natoms'th atom. */
void
gmx_ana_index_init_simple(gmx_ana_index_t *g, int natoms, char *name);
/** Frees memory allocated for an index group. */
void
gmx_ana_index_deinit(gmx_ana_index_t *g);
/** Copies a \c gmx_ana_index_t. */
void
gmx_ana_index_copy(gmx_ana_index_t *dest, gmx_ana_index_t *src, gmx_bool bAlloc);

/** Writes out the contents of a index group. */
void
gmx_ana_index_dump(gmx_ana_index_t *g, int i, int maxn);

/** Checks whether all indices are between 0 and \p natoms. */
void
gmx_ana_index_check(gmx_ana_index_t *g, int natoms);
/** Checks whether an index group is sorted. */
gmx_bool
gmx_ana_index_check_sorted(gmx_ana_index_t *g);
/*@}*/

/*! \name Functions for set operations on gmx_ana_index_t
 */
/*@{*/
/** Sorts the indices within an index group. */
void
gmx_ana_index_sort(gmx_ana_index_t *g);
/** Checks whether two index groups are equal. */
gmx_bool
gmx_ana_index_equals(gmx_ana_index_t *a, gmx_ana_index_t *b);
/** Checks whether a sorted index group contains another sorted index group. */
gmx_bool
gmx_ana_index_contains(gmx_ana_index_t *a, gmx_ana_index_t *b);

/** Calculates the intersection between two sorted index groups. */
void
gmx_ana_index_intersection(gmx_ana_index_t *dest,
                           gmx_ana_index_t *a, gmx_ana_index_t *b);
/** Calculates the set difference between two sorted index groups. */
void
gmx_ana_index_difference(gmx_ana_index_t *dest,
                         gmx_ana_index_t *a, gmx_ana_index_t *b);
/** Calculates the size of the difference between two sorted index groups. */
int
gmx_ana_index_difference_size(gmx_ana_index_t *a, gmx_ana_index_t *b);
/** Calculates the union of two sorted index groups. */
void
gmx_ana_index_union(gmx_ana_index_t *dest,
                    gmx_ana_index_t *a, gmx_ana_index_t *b);
/** Merges two distinct sorted index groups. */
void
gmx_ana_index_merge(gmx_ana_index_t *dest,
                    gmx_ana_index_t *a, gmx_ana_index_t *b);
/** Calculates the intersection and the difference in one call. */
void
gmx_ana_index_partition(gmx_ana_index_t *dest1, gmx_ana_index_t *dest2,
                        gmx_ana_index_t *src, gmx_ana_index_t *g);
/*@}*/

/*! \name Functions for handling gmx_ana_indexmap_t and related things
 */
/*@{*/
/** Partition a group based on topology information. */
void
gmx_ana_index_make_block(t_blocka *t, t_topology *top, gmx_ana_index_t *g,
                         e_index_t type, gmx_bool bComplete);
/** Checks whether a group consists of full blocks. */
gmx_bool
gmx_ana_index_has_full_blocks(gmx_ana_index_t *g, t_block *b);
/** Checks whether a group consists of full blocks. */
gmx_bool
gmx_ana_index_has_full_ablocks(gmx_ana_index_t *g, t_blocka *b);
/** Checks whether a group consists of full residues/molecules. */
gmx_bool
gmx_ana_index_has_complete_elems(gmx_ana_index_t *g, e_index_t type, t_topology *top);

/** Initializes an empty index group mapping. */
void
gmx_ana_indexmap_clear(gmx_ana_indexmap_t *m);
/** Reserves memory for an index group mapping. */
void
gmx_ana_indexmap_reserve(gmx_ana_indexmap_t *m, int nr, int isize);
/** Initializes an index group mapping. */
GMX_LIBGMX_EXPORT
void
gmx_ana_indexmap_init(gmx_ana_indexmap_t *m, gmx_ana_index_t *g,
                      t_topology *top, e_index_t type);
/** Sets an index group mapping to be static. */
void
gmx_ana_indexmap_set_static(gmx_ana_indexmap_t *m, t_blocka *b);
/** Frees memory allocated for index group mapping. */
void
gmx_ana_indexmap_deinit(gmx_ana_indexmap_t *m);
/** Makes a deep copy of an index group mapping. */
void
gmx_ana_indexmap_copy(gmx_ana_indexmap_t *dest, gmx_ana_indexmap_t *src, gmx_bool bFirst);
/** Updates an index group mapping. */
GMX_LIBGMX_EXPORT
void
gmx_ana_indexmap_update(gmx_ana_indexmap_t *m, gmx_ana_index_t *g, gmx_bool bMaskOnly);
/*@}*/

#ifdef __cplusplus
}
#endif

#endif
