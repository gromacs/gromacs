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

#include <variant>
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

    using DirectList = std::vector<Entry>;
    using HashedList = gmx::HashedMap<Entry>;

    /*! \brief Constructor
     *
     * \param[in] numAtomsTotal  The total number of atoms in the system
     * \param[in] numAtomsLocal  An estimate of the number of home+communicated atoms
     */
    gmx_ga2la_t(int numAtomsTotal, int numAtomsLocal);
    ~gmx_ga2la_t() {}

    /*! \brief Inserts an entry, there should not already be an entry for \p a_gl
     *
     * \param[in]  a_gl   The global atom index
     * \param[in]  value  The value to set for this index
     */
    void insert(int a_gl, const Entry& value)
    {
        GMX_ASSERT(a_gl >= 0, "Only global atom indices >= 0 are supported");
        if (usingDirect())
        {
            auto& directList = *std::get_if<DirectList>(&data_);
            GMX_ASSERT(directList[a_gl].cell == -1, "The key to be inserted should not be present");
            directList[a_gl] = value;
        }
        else
        {
            auto& hashedList = *std::get_if<HashedList>(&data_);
            hashedList.insert(a_gl, value);
        }
    }

    //! Delete the entry for global atom a_gl
    void erase(int a_gl)
    {
        if (usingDirect())
        {
            auto& directList      = *std::get_if<DirectList>(&data_);
            directList[a_gl].cell = -1;
        }
        else
        {
            auto& hashedList = *std::get_if<HashedList>(&data_);
            hashedList.erase(a_gl);
        }
    }

    //! Returns a pointer to the entry when present, nullptr otherwise
    const Entry* find(int a_gl) const
    {
        if (usingDirect())
        {
            const auto& directList = *std::get_if<DirectList>(&data_);
            return (directList[a_gl].cell == -1) ? nullptr : &(directList[a_gl]);
        }
        else
        {
            const auto& hashedList = *std::get_if<HashedList>(&data_);
            return hashedList.find(a_gl);
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
        if (usingDirect())
        {
            auto& directList = *std::get_if<DirectList>(&data_);
            GMX_ASSERT(directList[a_gl].cell >= 0, "a_gl should be present");
            return directList[a_gl];
        }
        else
        {
            auto&  hashedList = *std::get_if<HashedList>(&data_);
            Entry* search     = hashedList.find(a_gl);
            GMX_ASSERT(search, "a_gl should be present");
            return *search;
        }
    }

    /*! \brief Clear all the entries in the list.
     *
     * Note that this might use OpenMP threading, so it should not be called from within an OpenMP region.
     *
     * \param[in] resizeHashTable  When true the hash table is optimized based on the current number of entries stored
     */
    void clear(bool resizeHashTable);

private:
    //! Returns whether we are using the direct list
    bool usingDirect() const { return std::holds_alternative<DirectList>(data_); }

    /*! \brief Variant for storing the global atom index to local atom information map
     *
     * Direct list uses std::vector<Entry> of size number of global atoms
     * Hashed map alternative for larger lists, indexed by global atom index
     */
    std::variant<
            //! Direct list of local atom information of size number of global atoms
            std::vector<Entry>,
            //! A Hashed map of local atom information, indexed by global atom index
            gmx::HashedMap<Entry>>
            data_;
};

#endif
