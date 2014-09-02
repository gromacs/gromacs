/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2009,2010,2011,2012,2013,2014, by the GROMACS development team, led by
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
 * Mostly, these functions are used internally by the selection engine, but
 * it is necessary to use some of these functions in order to provide external
 * index groups to a gmx::SelectionCollection.
 * Some of the checking functions can be useful outside the selection engine to
 * check the validity of input groups.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inlibraryapi
 * \ingroup module_selection
 */
#ifndef GMX_SELECTION_INDEXUTIL_H
#define GMX_SELECTION_INDEXUTIL_H

#include <cstdio>

#include <string>

#include "gromacs/legacyheaders/types/simple.h"
#include "gromacs/topology/block.h"

struct t_topology;

/** Stores a set of index groups. */
struct gmx_ana_indexgrps_t;

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
struct gmx_ana_index_t
{
    /** Number of atoms. */
    int                 isize;
    /** List of atoms. */
    atom_id            *index;
    /** Number of items allocated for \p index. */
    int                 nalloc_index;
};

/*! \brief
 * Data structure for calculating index group mappings.
 */
struct gmx_ana_indexmap_t
{
    /** Type of the mapping. */
    e_index_t           type;
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
     * This array provides IDs for the current index group.  Instead of a
     * zero-based mapping that \p refid provides, the values from the \p orgid
     * array are used, thus allowing the mapping to be customized.
     * In other words, `mapid[i] = orgid[refid[i]]`.
     * If \p bMaskOnly is provided to gmx_ana_indexmap_update(), this array
     * equals \p orgid.
     */
    int                *mapid;
    /*! \brief
     * Mapped block structure.
     *
     * A block structure that corresponds to the current index group.
     * \c mapb.nra and \c mapb.a correspond to the last mapped index group.
     */
    t_blocka            mapb;

    /*! \brief
     * Customizable ID numbers for the original blocks.
     *
     * This array has \p b.nr elements, each defining an original ID number for
     * a block in \p b (i.e., in the original group passed to
     * gmx_ana_indexmap_init()).
     * These are initialized in gmx_ana_indexmap_init() based on the type:
     *  - \ref INDEX_ATOM : the atom indices
     *  - \ref INDEX_RES :  the residue indices
     *  - \ref INDEX_MOL :  the molecule indices
     *
     * All the above numbers are zero-based.
     * After gmx_ana_indexmap_init(), the caller is free to change these values
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
     * true if the current reference IDs are for the whole group (internal use only).
     *
     * This is used internally to optimize the evaluation such that
     * gmx_ana_indexmap_update() does not take any time if the group is
     * actually static.
     */
    bool                bStatic;
};


/*! \name Functions for handling gmx_ana_indexgrps_t
 */
/*@{*/
/** Reads index groups from a file or constructs them from topology. */
void
gmx_ana_indexgrps_init(gmx_ana_indexgrps_t **g, t_topology *top,
                       const char *fnm);
/** Frees memory allocated for index groups. */
void
gmx_ana_indexgrps_free(gmx_ana_indexgrps_t *g);
/** Returns true if the index group structure is emtpy. */
bool
gmx_ana_indexgrps_is_empty(gmx_ana_indexgrps_t *g);

/** Returns a pointer to an index group. */
gmx_ana_index_t *
gmx_ana_indexgrps_get_grp(gmx_ana_indexgrps_t *g, int n);
/** Extracts a single index group. */
bool
gmx_ana_indexgrps_extract(gmx_ana_index_t *dest, std::string *destName,
                          gmx_ana_indexgrps_t *src, int n);
/** Finds and extracts a single index group by name. */
bool
gmx_ana_indexgrps_find(gmx_ana_index_t *dest, std::string *destName,
                       gmx_ana_indexgrps_t *src, const char *name);

/** Writes out a list of index groups. */
void
gmx_ana_indexgrps_print(FILE *fp, gmx_ana_indexgrps_t *g, int maxn);
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
gmx_ana_index_set(gmx_ana_index_t *g, int isize, atom_id *index, int nalloc);
/** Creates a simple index group from the first to the \p natoms'th atom. */
void
gmx_ana_index_init_simple(gmx_ana_index_t *g, int natoms);
/** Frees memory allocated for an index group. */
void
gmx_ana_index_deinit(gmx_ana_index_t *g);
/** Copies a \c gmx_ana_index_t. */
void
gmx_ana_index_copy(gmx_ana_index_t *dest, gmx_ana_index_t *src, bool bAlloc);

/** Writes out the contents of a index group. */
void
gmx_ana_index_dump(FILE *fp, gmx_ana_index_t *g, int maxn);

/*! \brief
 * Returns maximum atom index that appears in an index group.
 *
 * \param[in]  g      Index group to query.
 * \returns    Largest atom index that appears in \p g, or zero if \p g is empty.
 */
int
gmx_ana_index_get_max_index(gmx_ana_index_t *g);
/** Checks whether an index group is sorted. */
bool
gmx_ana_index_check_sorted(gmx_ana_index_t *g);
/*! \brief
 * Checks whether an index group has atoms from a defined range.
 *
 * \param[in]  g      Index group to check.
 * \param[in]  natoms Largest atom number allowed.
 * \returns    true if all atoms in the index group are in the
 *     range 0 to \p natoms (i.e., no atoms over \p natoms are referenced).
 */
bool
gmx_ana_index_check_range(gmx_ana_index_t *g, int natoms);
/*@}*/

/*! \name Functions for set operations on gmx_ana_index_t
 */
/*@{*/
/** Sorts the indices within an index group. */
void
gmx_ana_index_sort(gmx_ana_index_t *g);
/** Checks whether two index groups are equal. */
bool
gmx_ana_index_equals(gmx_ana_index_t *a, gmx_ana_index_t *b);
/** Checks whether a sorted index group contains another sorted index group. */
bool
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
                         e_index_t type, bool bComplete);
/** Checks whether a group consists of full blocks. */
bool
gmx_ana_index_has_full_blocks(gmx_ana_index_t *g, t_block *b);
/** Checks whether a group consists of full blocks. */
bool
gmx_ana_index_has_full_ablocks(gmx_ana_index_t *g, t_blocka *b);
/** Checks whether a group consists of full residues/molecules. */
bool
gmx_ana_index_has_complete_elems(gmx_ana_index_t *g, e_index_t type, t_topology *top);

/** Initializes an empty index group mapping. */
void
gmx_ana_indexmap_clear(gmx_ana_indexmap_t *m);
/** Reserves memory for an index group mapping. */
void
gmx_ana_indexmap_reserve(gmx_ana_indexmap_t *m, int nr, int isize);
/** Initializes an index group mapping. */
void
gmx_ana_indexmap_init(gmx_ana_indexmap_t *m, gmx_ana_index_t *g,
                      t_topology *top, e_index_t type);
/*! \brief
 * Initializes `orgid` entries based on topology grouping.
 *
 * \param[in,out] m     Mapping structure to use/initialize.
 * \param[in]     top   Topology structure
 *     (can be NULL if not required for \p type).
 * \param[in]     type  Type of groups to use.
 * \returns  The number of groups of type \p type that were present in \p m.
 * \throws   InconsistentInputError if the blocks in \p m do not have a unique
 *     group (e.g., contain atoms from multiple residues with `type == INDEX_RES`).
 *
 * By default, the gmx_ana_indexmap_t::orgid fields are initialized to
 * atom/residue/molecule indices from the topology (see documentation for the
 * struct for more details).
 * This function can be used to set the field to a zero-based group index
 * instead.  The first block will always get `orgid[0] = 0`, and all following
 * blocks that belong to the same residue/molecule (\p type) will get the same
 * index.  Each time a block does not belong to the same group, it will get the
 * next available number.
 * If `type == INDEX_ATOM`, the `orgid` field is initialized as 0, 1, ...,
 * independent of whether the blocks are single atoms or not.
 *
 * Strong exception safety guarantee.
 */
int
gmx_ana_indexmap_init_orgid_group(gmx_ana_indexmap_t *m, t_topology *top,
                                  e_index_t type);
/** Sets an index group mapping to be static. */
void
gmx_ana_indexmap_set_static(gmx_ana_indexmap_t *m, t_blocka *b);
/** Frees memory allocated for index group mapping. */
void
gmx_ana_indexmap_deinit(gmx_ana_indexmap_t *m);
/** Makes a deep copy of an index group mapping. */
void
gmx_ana_indexmap_copy(gmx_ana_indexmap_t *dest, gmx_ana_indexmap_t *src, bool bFirst);
/** Updates an index group mapping. */
void
gmx_ana_indexmap_update(gmx_ana_indexmap_t *m, gmx_ana_index_t *g, bool bMaskOnly);
/*@}*/

#endif
