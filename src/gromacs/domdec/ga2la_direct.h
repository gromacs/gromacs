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
 * indices using a direct list with the size of the total number of atoms.
 * The functions are performance critical and should be inlined.
 *
 * \inlibraryapi
 * \ingroup module_domdec
 *
 * \author Berk Hess <hess@kth.se>
 *
 */
#ifndef GMX_DOMDEC_GA2LA_DIRECT_H
#define GMX_DOMDEC_GA2LA_DIRECT_H

#include <algorithm>
#include <vector>

#include "gromacs/utility/basedefinitions.h"

/*! \libinternal \brief Global to local atom mapping
 *
 * Efficiently manages mapping from global atom indices to local atom indices
 * in the domain decomposition.
 */
class Ga2laDirect
{
    private:
        /*! \libinternal \brief Structure for the local atom info for a plain list */
        struct directEntry
        {
            int  la   = -1; /**< The local atom index */
            int  cell = -1; /**< The DD zone index for neighboring domains, zone+zone otherwise */
        };

    public:
        /*! \brief Constructor
         *
         * \param[in] numAtomsTotal  The total number of atoms in the system
         * \param[in] numAtomsLocal  An estimate of the number of home+communicated atoms
         */
        Ga2laDirect(int numAtomsTotal) :
            table_(numAtomsTotal)
        {
        };

        /*! \brief Returns whether the direct list is available */
        bool isAvailable() const
        {
            return !table_.empty();
        }

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
            table_[a_gl].la   = a_loc;
            table_[a_gl].cell = cell;
        }

        /*! \brief Delete the entry for global atom a_gl
         *
         * \param[in]     a_gl  The global atom index
         */
        void unset(int a_gl)
        {
            table_[a_gl].cell = -1;
        }

        /*! \brief Change the local atom for entry with global atom a_gl
         *
         * \param[in] a_gl  The global atom index
         * \param[in] a_loc The new local atom index
         */
        void changeLocalAtom(int a_gl,
                             int a_loc)
        {
            table_[a_gl].la = a_loc;
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
            *a_loc = table_[a_gl].la;
            *cell  = table_[a_gl].cell;

            return (table_[a_gl].cell >= 0);
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
            *a_loc = table_[a_gl].la;

            return (table_[a_gl].cell == 0);
        }

        /*! \brief Returns whether the global atom a_gl is a home atom
         *
         * \param[in]  a_gl  The global atom index
         * \return whether the global atom a_gl is a home atom
         */
        bool isHome(int a_gl) const
        {
            return (table_[a_gl].cell == 0);
        }

        void clear()
        {
            for (directEntry &entry : table_)
            {
                entry.cell = -1;
            }
        }

    private:
        std::vector<directEntry> table_; /**< The direct list */
};

#endif
