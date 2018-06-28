/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2010,2014,2015,2017,2018, by the GROMACS development team, led by
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
 * Defines structures and functions for mapping from global to local atom
 * indices using a hash table.
 * The functions are performance critical and should be inlined.
 *
 * \inlibraryapi
 * \ingroup module_domdec
 *
 * \author Berk Hess <hess@kth.se>
 *
 */
#ifndef GMX_DOMDEC_GA2LA_HASHED_H
#define GMX_DOMDEC_GA2LA_HASHED_H

#include <algorithm>
#include <vector>

#include "gromacs/utility/basedefinitions.h"

/*! \libinternal \brief Global to local atom mapping
 *
 * Efficiently manages mapping from global atom indices to local atom indices
 * in the domain decomposition.
 */
class Ga2laHashed
{
    private:
        /*! \libinternal \brief Structure for the local atom info for a hash table */
        struct hashEntry
        {
            int  ga   = -1; /**< The global atom index */
            int  la   = -1; /**< The local atom index */
            int  cell = -1; /**< The DD zone index for neighboring domains, zone+zone otherwise */
            int  next = -1; /**< Index in the list of the next element with the same hash, -1 if none */
        };

    public:
        /*! \brief Constructor
         *
         * \param[in] baseTableSize  The size of the base table, optimal is around 2*numLocalAtoms
         */
        Ga2laHashed(int baseTableSize) :
            table_(baseTableSize),
            mod_(baseTableSize),
            startSpaceSearch_(baseTableSize)
        {
        };

        /*! \brief Returns whether the hashed table is available */
        bool isAvailable() const
        {
            return !table_.empty();
        };

        /*! \brief Sets the entry for global atom a_gl
         *
         * \param[in] a_gl  The global atom index
         * \param[in] a_loc The local atom index
         * \param[in] cell  The cell index
         */
        void set(int a_gl,
                 int a_loc,
                 int cell)
        {
            size_t ind = a_gl % mod_;

            if (table_[ind].ga >= 0)
            {
                /* Search the last entry in the linked list for this index */
                int ind_prev = ind;
                while (table_[ind_prev].next >= 0)
                {
                    ind_prev = table_[ind_prev].next;
                }
                /* Search for space in the array */
                ind = startSpaceSearch_;
                while (ind < table_.size() && table_[ind].ga >= 0)
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
            table_[ind].ga   = a_gl;
            table_[ind].la   = a_loc;
            table_[ind].cell = cell;
        }

        /*! \brief Delete the entry for global atom a_gl
         *
         * \param[in]     a_gl  The global atom index
         */
        void unset(int a_gl)
        {
            int ind_prev = -1;
            int ind      = a_gl % mod_;
            do
            {
                if (table_[ind].ga == a_gl)
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
                    table_[ind].ga   = -1;
                    table_[ind].cell = -1;
                    table_[ind].next = -1;

                    return;
                }
                ind_prev = ind;
                ind      = table_[ind].next;
            }
            while (ind >= 0);
        }

        /*! \brief Change the local atom for entry with global atom a_gl
         *
         * \param[in] a_gl  The global atom index
         * \param[in] a_loc The new local atom index
         */
        void changeLocalAtom(int a_gl,
                             int a_loc)
        {
            int ind = a_gl % mod_;
            do
            {
                if (table_[ind].ga == a_gl)
                {
                    table_[ind].la = a_loc;

                    return;
                }
                ind = table_[ind].next;
            }
            while (ind >= 0);
        }

        /*! \brief Returns whether the global atom a_gl available locally
         *
         * \param[in]  a_gl  The global atom index
         * \param[out] a_loc If the return value is TRUE, the local atom index
         * \param[out] cell  If the return value is TRUE, the zone or for atoms more than one cell away zone+nzone
         * \return if the global atom a_gl available locally
         */
        bool get(int  a_gl,
                 int *a_loc,
                 int *cell) const
        {
            int ind = a_gl % mod_;
            do
            {
                if (table_[ind].ga == a_gl)
                {
                    *a_loc = table_[ind].la;
                    *cell  = table_[ind].cell;

                    return true;
                }
                ind = table_[ind].next;
            }
            while (ind >= 0);

            return false;
        }

        /*! \brief Returns whether the global atom a_gl is a home atom
         *
         * \param[in]  a_gl  The global atom index
         * \param[out] a_loc If the return value is true, the local atom index
         * \return whether the global atom a_gl is a home atom
         */
        bool getHome(int  a_gl,
                     int *a_loc) const
        {
            int ind = a_gl % mod_;
            do
            {
                if (table_[ind].ga == a_gl)
                {
                    if (table_[ind].cell == 0)
                    {
                        *a_loc = table_[ind].la;

                        return true;
                    }
                    else
                    {
                        return false;
                    }
                }
                ind = table_[ind].next;
            }
            while (ind >= 0);

            return false;
        }

        /*! \brief Returns whether the global atom a_gl is a home atom
         *
         * \param[in]  a_gl  The global atom index
         * \return whether the global atom a_gl is a home atom
         */
        bool isHome(int a_gl) const
        {
            int ind = a_gl % mod_;
            do
            {
                if (table_[ind].ga == a_gl)
                {
                    return (table_[ind].cell == 0);
                }
                ind = table_[ind].next;
            }
            while (ind >= 0);

            return false;
        }

        /*! \brief Clear all the entries in the list */
        void clear()
        {
            for (hashEntry &entry : table_)
            {
                entry.ga   = -1;
                entry.next = -1;
            }
            startSpaceSearch_ = mod_;
        }

    private:
        std::vector<hashEntry> table_;            /**< The hash table list */
        const int              mod_;              /**< The hash size */
        int                    startSpaceSearch_; /**< Index in lal at which to start looking for empty space */
};

#endif
