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
 * indices. The functions are performance critical and should be inlined.
 *
 * \inlibraryapi
 * \ingroup module_domdec
 *
 * \author Berk Hess <hess@kth.se>
 *
 */
#ifndef GMX_DOMDEC_GA2LA_H
#define GMX_DOMDEC_GA2LA_H

#include <algorithm>
#include <vector>
#include <unordered_map>

#include "gromacs/utility/basedefinitions.h"
#include "ga2la_direct.h"
#include "ga2la_hashed.h"

/*! \libinternal \brief Global to local atom mapping
 *
 * Efficiently manages mapping from global atom indices to local atom indices
 * in the domain decomposition.
 */
class gmx_ga2la_t
{
    public:
        /*! \libinternal \brief Structure for the local atom info */
        struct Entry
        {
            int  la;   /**< The local atom index */
            int  cell; /**< The DD zone index for neighboring domains, zone+zone otherwise */
        };

        /*! \brief Constructor
         *
         * \param[in] numAtomsTotal  The total number of atoms in the system
         * \param[in] numAtomsLocal  An estimate of the number of home+communicated atoms
         */
        gmx_ga2la_t(int numAtomsTotal,
                    int numAtomsLocal);
        ~gmx_ga2la_t()
        {
            bUsingDirect ? data_.direct.~vector() : data_.hashed.~unordered_map();
        }

        Entry &operator[](int a_gl)
        {
            return bUsingDirect ? data_.direct[a_gl] : data_.hashed[a_gl];
        }
        /*! \brief Delete the entry for global atom a_gl
         *
         * \param[in]     a_gl  The global atom index
         */
        void erase(int a_gl)
        {
            if (bUsingDirect)
            {
                data_.direct[a_gl].cell = -1;
            }
            else
            {
                data_.hashed.erase(a_gl);
            }
        }

        const Entry* find(int a_gl) const
        {
            if (bUsingDirect)
            {
                return (data_.direct[a_gl].cell == -1) ? nullptr : &(data_.direct[a_gl]);
            }
            else
            {
                auto it = data_.hashed.find(a_gl);
                return (it == data_.hashed.end()) ? nullptr : &(it->second);
            }
        }

        //! Returns the local atom if it is a home atom. nullptr otherwise.
        const int* findHome(int  a_gl) const
        {
            const Entry* const e = find(a_gl);
            return (e && e->cell == 0) ? &(e->la) : nullptr;
        }

        /*! \brief Clear all the entries in the list */
        void clear()
        {
            if (bUsingDirect)
            {
                for (Entry &entry : data_.direct) {entry.cell = -1; }
            }
            else
            {
                data_.hashed.clear();
            }
        }

    private:
        union Data
        {
            std::vector<Entry>             direct;
            std::unordered_map<int, Entry> hashed;
            // constructor and destructor function in parent class
            Data()  {}
            ~Data() {}
        } data_;
        bool bUsingDirect;
};

#endif
