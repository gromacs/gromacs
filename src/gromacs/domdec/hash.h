/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2014,2015,2017,2018, by the GROMACS development team, led by
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
/*! \internal \file
 * \brief
 * This file declares functions for a simple hash map used by domain
 * decomposition.
 *
 * It is limited to integer keys and integer values. The purpose is highest
 * efficiency and lowest memory usage possible.  Thus the code is in a header,
 * so it can be inlined where it is used.
 *
 * Note that we could use std::unordered_map for this functionality,
 * but that is an order of magnitude slower and impacts the DD performance.
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_domdec
 */

#ifndef GMX_DOMDEC_HASH_H
#define GMX_DOMDEC_HASH_H

#include <cstdio>

//! Forward declation
static void gmx_hash_resize(gmx_hash_t * /*hash*/, int /*nkey_used_estimate*/);

/*! \internal \brief Hashing key-generation helper struct */
struct gmx_hash_e_t
{
    public:
        //! The (unique) key for storing/looking up a value
        int  key  = -1;
        //! The value belonging to key
        int  val  = -1;
        //! Index for the next element in the array with indentical value key%mod, -1 if there is no next element
        int  next = -1;
};

/*! \internal \brief Hashing helper struct */
struct gmx_hash_t
{
    public:
        //! Constructor
        gmx_hash_t(int numKeysUsedEstimate) :
            nkey(0)
        {
            gmx_hash_resize(this, numKeysUsedEstimate);
        }

        //! Keys are looked up by first checking array index key%mod in hash
        int                       mod;
        //! mask=log2(mod), used to replace a % by the faster & operation
        int                       mask;
        //! The actual array containing the keys, values and next indices
        std::vector<gmx_hash_e_t> hash;
        //! The number of keys stored
        int                       nkey;
        //! Index in hash where we should start searching for space to store a new key/value
        int                       start_space_search;
};

//! Clear all the entries in the hash table.
static void gmx_hash_clear(gmx_hash_t *hash)
{
    for (gmx_hash_e_t &entry : hash->hash)
    {
        /* Note: clearing entry.val is not needed */
        entry.key  = -1;
        entry.next = -1;
    }
    hash->start_space_search = hash->mod;

    hash->nkey = 0;
}

//! Reallocate hash table data structures.
static void gmx_hash_resize(gmx_hash_t *hash, int nkey_used_estimate)
{
    GMX_RELEASE_ASSERT(hash->nkey == 0, "Table needs to be empty for resize");

    /* Memory requirements:
     * nkey_used_est*(2+1-2(1-e^-1/2))*3 ints
     * where nkey_used_est is the local number of keys used.
     *
     * Make the direct list twice as long as the number of local keys.
     * The fraction of entries in the list with:
     * 0   size lists: e^-f
     * >=1 size lists: 1 - e^-f
     * where f is: the #keys / mod
     * The fraction of keys not in the direct list is: 1-1/f(1-e^-f).
     * The optimal table size is roughly double the number of keys.
     */
    /* Make the hash table a power of 2 and at least double the number of keys */
    hash->mod = 4;
    while (2*nkey_used_estimate > hash->mod)
    {
        hash->mod *= 2;
    }
    hash->mask   = hash->mod - 1;
    hash->hash.resize(hash->mod);
}

/*! \brief Clear all the entries in the hash table.
 *
 * With the current number of keys check if the table size is still
 * good, if not optimize it with the current number of keys.
 */
static inline void gmx_hash_clear_and_optimize(gmx_hash_t *hash)
{
    gmx_hash_clear(hash);

    /* Resize the hash table when the occupation is < 1/4 or > 2/3 */
    if (hash->nkey > 0 &&
        (4*hash->nkey < hash->mod || 3*hash->nkey > 2*hash->mod))
    {
        gmx_hash_resize(hash, hash->nkey);
    }
}

//! Set the hash entry for key to value.
static void gmx_hash_set(gmx_hash_t *hash, int key, int value)
{
    unsigned int ind = (key & hash->mask);

    if (hash->hash[ind].key >= 0)
    {
        /* Search the last entry in the linked list for this index */
        unsigned int ind_prev = ind;
        while (hash->hash[ind_prev].next >= 0)
        {
            ind_prev = hash->hash[ind_prev].next;
        }
        /* Search for space in the array */
        ind = hash->start_space_search;
        while (ind < hash->hash.size() && hash->hash[ind].key >= 0)
        {
            ind++;
        }
        /* If we are at the end of the list we need to increase the size */
        if (ind == hash->hash.size())
        {
            hash->hash.resize(hash->hash.size() + 1);
        }
        hash->hash[ind_prev].next = ind;

        hash->start_space_search = ind + 1;
    }
    hash->hash[ind].key = key;
    hash->hash[ind].val = value;

    hash->nkey++;
}

//! Delete the hash entry for key.
static inline void gmx_hash_del(gmx_hash_t *hash, int key)
{
    int ind, ind_prev;

    ind_prev = -1;
    ind      = key & hash->mask;
    do
    {
        if (hash->hash[ind].key == key)
        {
            if (ind_prev >= 0)
            {
                hash->hash[ind_prev].next = hash->hash[ind].next;

                /* This index is a linked entry, so we free an entry.
                 * Check if we are creating the first empty space.
                 */
                if (ind < hash->start_space_search)
                {
                    hash->start_space_search = ind;
                }
            }
            hash->hash[ind].key  = -1;
            hash->hash[ind].val  = -1;
            hash->hash[ind].next = -1;

            hash->nkey--;

            return;
        }
        ind_prev = ind;
        ind      = hash->hash[ind].next;
    }
    while (ind >= 0);
}

//! Change the value for present hash entry for key.
static inline void gmx_hash_change_value(gmx_hash_t *hash, int key, int value)
{
    int ind = (key & hash->mask);

    do
    {
        if (hash->hash[ind].key == key)
        {
            hash->hash[ind].val = value;

            return;
        }
        ind = hash->hash[ind].next;
    }
    while (ind >= 0);
}

//! Change the hash value if already set, otherwise set the hash value.
static inline void gmx_hash_change_or_set(gmx_hash_t *hash, int key, int value)
{
    int ind = (key & hash->mask);

    do
    {
        if (hash->hash[ind].key == key)
        {
            hash->hash[ind].val = value;

            return;
        }
        ind = hash->hash[ind].next;
    }
    while (ind >= 0);

    gmx_hash_set(hash, key, value);
}

//! Returns if the key is present, if the key is present *value is set.
static inline gmx_bool gmx_hash_get(const gmx_hash_t *hash, int key, int *value)
{
    int ind = (key & hash->mask);

    do
    {
        if (hash->hash[ind].key == key)
        {
            *value = hash->hash[ind].val;

            return TRUE;
        }
        ind = hash->hash[ind].next;
    }
    while (ind >= 0);

    return FALSE;
}

//! Returns the value or -1 if the key is not present.
static int gmx_hash_get_minone(const gmx_hash_t *hash, int key)
{
    int ind = (key & hash->mask);

    do
    {
        if (hash->hash[ind].key == key)
        {
            return hash->hash[ind].val;
        }
        ind = hash->hash[ind].next;
    }
    while (ind >= 0);

    return -1;
}

#endif
