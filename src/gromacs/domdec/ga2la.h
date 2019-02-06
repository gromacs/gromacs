/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2010,2014,2015,2017,2018,2019, by the GROMACS development team, led by
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

#include <vector>

#include "gromacs/compat/variant.h"
#include "gromacs/domdec/hashedmap.h"
#include "gromacs/utility/gmxassert.h"

/*! \libinternal \brief Global to local atom mapping
 *
 * Used for efficient mapping from global atom indices to local atom indices
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

        /*! \brief Inserts an entry, there should not already be an entry for \p a_gl
         *
         * \param[in]  a_gl   The global atom index
         * \param[in]  value  The value to set for this index
         */
        void insert(int          a_gl,
                    const Entry &value)
        {
            GMX_ASSERT(a_gl >= 0, "Only global atom indices >= 0 are supported");
            if (isDirect())
            {
                GMX_ASSERT(direct()[a_gl].cell == -1, "The key to be inserted should not be present");
                direct()[a_gl] = value;
            }
            else
            {
                hashed().insert(a_gl, value);
            }
        }

        //! Delete the entry for global atom a_gl
        void erase(int a_gl)
        {
            if (isDirect())
            {
                direct()[a_gl].cell = -1;
            }
            else
            {
                hashed().erase(a_gl);
            }
        }

        //! Returns a pointer to the entry when present, nullptr otherwise
        const Entry* find(int a_gl) const
        {
            if (isDirect())
            {
                return (direct()[a_gl].cell == -1) ? nullptr : &(direct()[a_gl]);
            }
            else
            {
                return (hashed().find(a_gl));
            }
        }

        //! Returns the local atom index if it is a home atom, nullptr otherwise
        const int* findHome(int a_gl) const
        {
            const Entry* const e = find(a_gl);
            return (e && e->cell == 0) ? &(e->la) : nullptr;
        }

        /*! \brief Returns a reference to the entry for a_gl
         *
         * A non-release assert checks that a_gl is present.
         */
        Entry &at(int a_gl)
        {
            if (isDirect())
            {
                GMX_ASSERT(direct()[a_gl].cell >= 0, "a_gl should be present");
                return direct()[a_gl];
            }
            else
            {
                Entry *search = hashed().find(a_gl);
                GMX_ASSERT(search, "a_gl should be present");
                return *search;
            }
        }

        //! Clear all the entries in the list.
        void clear()
        {
            if (isDirect())
            {
                for (Entry &entry : direct()) {entry.cell = -1; }
            }
            else
            {
                hashed().clear();
            }
        }

    private:
        bool isDirect() const { return gmx::compat::holds_alternative<std::vector<Entry> >(data_); }
        std::vector<Entry> &direct() { return gmx::compat::get<std::vector<Entry> >(data_); }
        const std::vector<Entry> &direct() const { return gmx::compat::get<std::vector<Entry> >(data_); }
        gmx::HashedMap<Entry> &hashed() { return gmx::compat::get<gmx::HashedMap<Entry> >(data_); }
        const gmx::HashedMap<Entry> &hashed() const { return gmx::compat::get<gmx::HashedMap<Entry> >(data_); }
        gmx::compat::variant<std::vector<Entry>, gmx::HashedMap<Entry> > data_;
};

#endif
