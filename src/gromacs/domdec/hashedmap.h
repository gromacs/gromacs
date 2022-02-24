/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
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
/*! \libinternal \file
 * \brief
 * Defines structures and functions for mapping from keys to entries
 * indices using a hash table.
 * The functions are performance critical and should be inlined.
 *
 * \inlibraryapi
 * \ingroup module_domdec
 *
 * \author Berk Hess <hess@kth.se>
 *
 */
#ifndef GMX_DOMDEC_HASHEDMAP_H
#define GMX_DOMDEC_HASHEDMAP_H

#include <climits>

#include <algorithm>
#include <utility>
#include <vector>

#include "gromacs/compat/utility.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/exceptions.h"

namespace gmx
{

/*! \libinternal \brief Unordered key to value mapping
 *
 * Efficiently manages mapping from integer keys to values.
 * Note that this basically implements a subset of the functionality of
 * std::unordered_map, but is an order of magnitude faster.
 */
template<class T>
class HashedMap
{
private:
    /*! \libinternal \brief Structure for the key/value hash table */
    struct hashEntry
    {
        int key = -1;  /**< The key */
        T   value;     /**< The value(s) */
        int next = -1; /**< Index in the list of the next element with the same hash, -1 if none */
    };

    /*! \brief The table size is set to at least this factor time the nr of keys */
    static constexpr float c_relTableSizeSetMin = 1.5;
    /*! \brief Threshold for increasing the table size */
    static constexpr float c_relTableSizeThresholdMin = 1.3;
    /*! \brief Threshold for decreasing the table size */
    static constexpr float c_relTableSizeThresholdMax = 3.5;

    /*! \brief Resizes the table
     *
     * \param[in] numElementsEstimate  An estimate of the number of elements that will be stored
     */
    void resize(int numElementsEstimate)
    {
        GMX_RELEASE_ASSERT(numElements_ == 0, "Table needs to be empty for resize");

        /* The fraction of table entries with 0   size lists is e^-f.
         * The fraction of table entries with >=1 size lists is 1 - e^-f
         * where f is: the #elements / tableSize
         * The fraction of elements not in the direct list is: 1 - (1 - e^-f)/f.
         * Thus the optimal table size is roughly double #elements.
         */
        /* Make the hash table a power of 2 and at least 1.5 * #elements */
        int tableSize = 64;
        while (tableSize <= INT_MAX / 2
               && static_cast<float>(numElementsEstimate) * c_relTableSizeSetMin > tableSize)
        {
            tableSize *= 2;
        }
        table_.resize(tableSize);

        /* Table size is a power of 2, so a binary mask gives the hash */
        bitMask_                        = tableSize - 1;
        startIndexForSpaceForListEntry_ = tableSize;
    }

public:
    /*! \brief Constructor
     *
     * \param[in] numElementsEstimate  An estimate of the number of elements that will be stored, used for optimizing initial performance
     *
     * Note that the estimate of the number of elements is only relevant
     * for the performance up until the first call to clear(), after which
     * table size is optimized based on the actual number of elements.
     */
    HashedMap(int numElementsEstimate) { resize(numElementsEstimate); }

    /*! \brief Returns the number of elements */
    int size() const { return numElements_; }

    /*! \brief Returns the number of buckets, i.e. the number of possible hashes */
    int bucket_count() const { return bitMask_ + 1; }

private:
    /*! \brief Inserts or assigns a key and value
     *
     * \tparam    allowAssign  Sets whether assignment of a key that is present is allowed
     * \param[in] key          The key for the entry
     * \param[in] value        The value for the entry
     * \throws InvalidInputError from a debug build when attempting to insert a duplicate key with \p allowAssign=true
     */
    // cppcheck-suppress unusedPrivateFunction
    template<bool allowAssign>
    void insert_assign(int key, const T& value)
    {
        size_t ind = (key & bitMask_);

        if (table_[ind].key >= 0)
        {
            /* Loop over the entries for this hash.
             * If we find the matching key, return the value.
             */
            int ind_prev = ind;
            if (table_[ind_prev].key == key)
            {
                if (allowAssign)
                {
                    table_[ind].value = value;
                    return;
                }
                else
                {
// Note: This is performance critical, so we only throw in debug mode
#ifndef NDEBUG
                    GMX_THROW(InvalidInputError("Attempt to insert duplicate key"));
#endif
                }
            }
            while (table_[ind_prev].next >= 0)
            {
                ind_prev = table_[ind_prev].next;
                if (table_[ind_prev].key == key)
                {
                    if (allowAssign)
                    {
                        table_[ind_prev].value = value;
                        return;
                    }
                    else
                    {
#ifndef NDEBUG
                        GMX_THROW(InvalidInputError("Attempt to insert duplicate key"));
#endif
                    }
                }
            }
            /* Search for space in table_ */
            ind = startIndexForSpaceForListEntry_;
            while (ind < table_.size() && table_[ind].key >= 0)
            {
                ind++;
            }
            /* If we are at the end of the list we need to increase the size */
            if (ind == table_.size())
            {
                table_.resize(table_.size() + 1);
            }
            table_[ind_prev].next = ind;

            startIndexForSpaceForListEntry_ = ind + 1;
        }

        table_[ind].key   = key;
        table_[ind].value = value;

        numElements_ += 1;
    }

public:
    /*! \brief Inserts entry, key should not already be present
     *
     * \param[in] key    The key for the entry
     * \param[in] value  The value for the entry
     * \throws InvalidInputError from a debug build when attempting to inser         */
    void insert(int key, const T& value) { insert_assign<false>(key, value); }

    /*! \brief Inserts an entry when the key is not present, otherwise sets the value
     *
     * \param[in] key    The key for the entry
     * \param[in] value  The value for the entry
     */
    void insert_or_assign(int key, const T& value) { insert_assign<true>(key, value); }

    /*! \brief Delete the entry for key \p key, when present
     *
     * \param[in] key  The key
     */
    void erase(int key)
    {
        int ind_prev = -1;
        int ind      = (key & bitMask_);
        do
        {
            if (table_[ind].key == key)
            {
                if (ind_prev >= 0)
                {
                    table_[ind_prev].next = table_[ind].next;

                    /* This index is a linked entry, so we free an entry.
                     * Check if we are creating the first empty space.
                     */
                    if (ind < startIndexForSpaceForListEntry_)
                    {
                        startIndexForSpaceForListEntry_ = ind;
                    }
                }
                table_[ind].key  = -1;
                table_[ind].next = -1;

                numElements_ -= 1;

                return;
            }
            ind_prev = ind;
            ind      = table_[ind].next;
        } while (ind >= 0);
    }

    /*! \brief Returns a pointer to the value for the given key or nullptr when not present
     *
     * \todo Use std::as_const when CUDA 11 is a requirement.
     *
     * \param[in] key  The key
     * \return a pointer to value for the given key or nullptr when not present
     */
    T* find(int key) { return const_cast<T*>(gmx::compat::as_const(*this).find(key)); }

    /*! \brief Returns a pointer to the value for the given key or nullptr when not present
     *
     * \param[in] key  The key
     * \return a pointer to value for the given key or nullptr when not present
     */
    const T* find(int key) const
    {
        int ind = (key & bitMask_);
        do
        {
            if (table_[ind].key == key)
            {
                return &table_[ind].value;
            }
            ind = table_[ind].next;
        } while (ind >= 0);

        return nullptr;
    }

    //! Clear all the entries in the list
    void clear()
    {
        for (hashEntry& entry : table_)
        {
            entry.key  = -1;
            entry.next = -1;
        }
        startIndexForSpaceForListEntry_ = bucket_count();
        numElements_                    = 0;
    }

    /*! \brief Clear all the entries in the list and resizes the hash table
     *
     * Optimizes the size of the hash table based on the current
     * number of elements stored.
     */
    void clearAndResizeHashTable()
    {
        const int oldNumElements = numElements_;

        clear();

        /* Resize the hash table when the occupation is far from optimal.
         * Do not resize with 0 elements to avoid minimal size when clear()
         * is called twice in a row.
         */
        if (oldNumElements > 0
            && (oldNumElements * c_relTableSizeThresholdMax < bucket_count()
                || oldNumElements * c_relTableSizeThresholdMin > bucket_count()))
        {
            resize(oldNumElements);
        }
    }

private:
    /*! \brief The hash table list */
    std::vector<hashEntry> table_;
    /*! \brief The bit mask for computing the hash of a key */
    int bitMask_ = 0;
    /*! \brief Index in table_ at which to start looking for empty space for a new linked list entry */
    int startIndexForSpaceForListEntry_ = 0;
    /*! \brief The number of elements currently stored in the table */
    int numElements_ = 0;
};

} // namespace gmx

#endif
