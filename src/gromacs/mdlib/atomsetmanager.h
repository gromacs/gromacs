/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
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
 * \brief
 * Declares gmx::AtomSetManager
 *
 * \author Christian Blau <cblau@gwdg.de>
 */

#ifndef GMX_ATOMSETMANAGER_H
#define GMX_ATOMSETMANAGER_H

#include <memory>
#include <map>

#include "atomset.h"
#include "gromacs/utility/classhelpers.h"


struct gmx_ga2la_t;

namespace gmx
{

/*! \internal \brief
 * Manage sets of atoms.
 * A set of atoms is a collection of local, global and collective indices of atoms with a name.
 * Triggers update of indices from domain decomposition code for all registered atom sets if run in parallel.
 *
 */
class AtomSetManager
{
    public:

        /*! \brief Construct an atom set manager.
         *
         * \param[in] bParallel atom set manager needs to know if it manages sets in parallel
         */
        AtomSetManager(const bool bParallel);
        ~AtomSetManager();

        /*! \brief Add a new atom set to be managed that is identified with a name and tells which atoms should be managed.
         *
         * \param[in] atom_set_name The name of the atom set
         * \param[in] nat The number of atoms in the atom set
         * \param[in] ind The indices of the atoms in the atom set
         *
         * \throws InternalError if atom set with the requested name already exists
         */
        void add(const std::string &atom_set_name, const int number_of_atoms, const int *index);
        /*! \brief Return a constant reference to an atom set.
         *
         * Atom sets may be altered only through the atom set manager.
         * \param [in] atom_set_name Name of the atom set to be retrieved.
         *
         * \throws InternalError if no atom set with the requested name exists
         */
        const AtomSet &get(const std::string &atom_set_name) const;

        /*! \brief Erase and destroy an atom set.
         *
         * Atoms will still be part of the simulation, only their indices will no longer be managed as atom set.
         *
         * \param [in] atom_set_name Name of the atom set to be erased.
         *
         * \throws InternalError if no atom set with the requested name exists
         */
        void erase(const std::string &atom_set_name);

        void set_indices_in_domain_decomposition(const gmx_ga2la_t  *ga2la);

    private:
        class Impl;
        PrivateImplPointer<Impl> impl_;
};



} // namespace gmx

#endif
