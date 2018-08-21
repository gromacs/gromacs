/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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

#include <algorithm>
#include <vector>

#include "gromacs/compat/utility.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/exceptions.h"

namespace gmx
{

/*! \libinternal \brief Unordered key to value mapping
 *
 * Efficiently manages mapping from non-negative integer keys to values.
 * Note that this basically implements a subset of the functionality of
 * std::unordered_map, but is an order of magnitude faster.
 */
template <class T>
class HashedMap
{
    private:
        /*! \libinternal \brief Structure for the key/value hash table */
        struct hashEntry
        {
            int  key   = -1; /**< The key */
            T    value;      /**< The value(s) */
            int  next = -1;  /**< Index in the list of the next element with the same hash, -1 if none */
        };

    public:
        /*! \brief Constructor
         *
         * \param[in] baseTableSize  The size of the base table, optimal is around twice the number of expected entries
         */
        HashedMap(int baseTableSize) :
            table_(baseTableSize),
            mod_(baseTableSize),
            startSpaceSearch_(baseTableSize)
        {
        }

        /*! \brief Inserts entry, key should not already be present (throws in debug build))
         *
         * \param[in] key    The key for the entry
         * \param[in] value  The value for the entry
         */
        void insert(int      key,
                    const T &value)
        {
            // Note: This is performance critical, so we only throw in debug mode
#ifndef NDEBUG
            if (key < 0)
            {
                GMX_THROW(InvalidInputError("Invalid key value"));
            }
#endif

            size_t ind = key % mod_;

            if (table_[ind].key >= 0)
            {
                /* Loop over the entries for this hash.
                 * If we find the matching key, return the value.
                 */
                int ind_prev = ind;
#ifndef NDEBUG
                if (table_[ind_prev].key == key)
                {
                    GMX_THROW(InvalidInputError("Attempt to insert duplicate key"));
                }
#endif
                while (table_[ind_prev].next >= 0)
                {
                    ind_prev = table_[ind_prev].next;
#ifndef NDEBUG
                    if (table_[ind_prev].key == key)
                    {
                        GMX_THROW(InvalidInputError("Attempt to insert duplicate key"));
                    }
#endif
                }
                /* Search for space in the array */
                ind = startSpaceSearch_;
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

                startSpaceSearch_ = ind + 1;
            }

            table_[ind].key   = key;
            table_[ind].value = value;
        }

        /*! \brief Delete the entry for key \p key
         *
         * \param[in] key  The key
         */
        void erase(int key)
        {
            int ind_prev = -1;
            int ind      = key % mod_;
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
                        if (ind < startSpaceSearch_)
                        {
                            startSpaceSearch_ = ind;
                        }
                    }
                    table_[ind].key  = -1;
                    table_[ind].next = -1;

                    return;
                }
                ind_prev = ind;
                ind      = table_[ind].next;
            }
            while (ind >= 0);
        }

        /*! \brief Returns a pointer to the value for the given key or nullptr when not present
         *
         * \param[in] key  The key
         * \return a pointer to value for the given key or nullptr when not present
         */
        T *find(int key) { return const_cast<T*>(gmx::compat::as_const(*this).find(key)); }

        /*! \brief Returns a pointer to the value for the given key or nullptr when not present
         *
         * \param[in] key  The key
         * \return a pointer to value for the given key or nullptr when not present
         */
        const T *find(int key) const
        {
            int ind = key % mod_;
            do
            {
                if (table_[ind].key == key)
                {
                    return &table_[ind].value;
                }
                ind = table_[ind].next;
            }
            while (ind >= 0);

            return nullptr;
        }

        /*! \brief Clear all the entries in the list */
        void clear()
        {
            for (hashEntry &entry : table_)
            {
                entry.key  = -1;
                entry.next = -1;
            }
            startSpaceSearch_ = mod_;
        }

    private:
        std::vector<hashEntry> table_;            /**< The hash table list */
        int                    mod_;              /**< The hash size */
        int                    startSpaceSearch_; /**< Index in lal at which to start looking for empty space */
};

} // namespace

#endif
