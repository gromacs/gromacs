/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
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
        int la;   /**< The local atom index */
        int cell; /**< The DD zone index for neighboring domains, zone+zone otherwise */
    };

    /*! \brief Constructor
     *
     * \param[in] numAtomsTotal  The total number of atoms in the system
     * \param[in] numAtomsLocal  An estimate of the number of home+communicated atoms
     */
    gmx_ga2la_t(int numAtomsTotal, int numAtomsLocal);
    ~gmx_ga2la_t() { usingDirect_ ? data_.direct.~vector() : data_.hashed.~HashedMap(); }

    /*! \brief Inserts an entry, there should not already be an entry for \p a_gl
     *
     * \param[in]  a_gl   The global atom index
     * \param[in]  value  The value to set for this index
     */
    void insert(int a_gl, const Entry& value)
    {
        GMX_ASSERT(a_gl >= 0, "Only global atom indices >= 0 are supported");
        if (usingDirect_)
        {
            GMX_ASSERT(data_.direct[a_gl].cell == -1,
                       "The key to be inserted should not be present");
            data_.direct[a_gl] = value;
        }
        else
        {
            data_.hashed.insert(a_gl, value);
        }
    }

    //! Delete the entry for global atom a_gl
    void erase(int a_gl)
    {
        if (usingDirect_)
        {
            data_.direct[a_gl].cell = -1;
        }
        else
        {
            data_.hashed.erase(a_gl);
        }
    }

    //! Returns a pointer to the entry when present, nullptr otherwise
    const Entry* find(int a_gl) const
    {
        if (usingDirect_)
        {
            return (data_.direct[a_gl].cell == -1) ? nullptr : &(data_.direct[a_gl]);
        }
        else
        {
            return (data_.hashed.find(a_gl));
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
    Entry& at(int a_gl)
    {
        if (usingDirect_)
        {
            GMX_ASSERT(data_.direct[a_gl].cell >= 0, "a_gl should be present");
            return data_.direct[a_gl];
        }
        else
        {
            Entry* search = data_.hashed.find(a_gl);
            GMX_ASSERT(search, "a_gl should be present");
            return *search;
        }
    }

    /*! \brief Clear all the entries in the list.
     *
     * \param[in] resizeHashTable  When true the hash table is optimized based on the current number of entries stored
     */
    void clear(const bool resizeHashTable)
    {
        if (usingDirect_)
        {
            for (Entry& entry : data_.direct)
            {
                entry.cell = -1;
            }
        }
        else if (resizeHashTable)
        {
            data_.hashed.clearAndResizeHashTable();
        }
        else
        {
            data_.hashed.clear();
        }
    }

private:
    union Data
    {
        std::vector<Entry>    direct;
        gmx::HashedMap<Entry> hashed;
        // constructor and destructor function in parent class
        Data() {}
        ~Data() {}
    } data_;
    const bool usingDirect_;
};

#endif
