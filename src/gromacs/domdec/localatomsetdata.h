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
 * \internal \brief
 * Declares gmx::internal::LocalAtomSetData.
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \ingroup module_domdec
 */
#ifndef GMX_LOCALATOMSETDATA_H
#define GMX_LOCALATOMSETDATA_H

#include <vector>

struct gmx_ga2la_t;

namespace gmx
{

namespace internal
{

/* \brief Internal class for storing and managing atom indices of an atom set.
 */
class LocalAtomSetData
{
    public:
        /*! \brief Store the data for an atom set with an index group.
         *
         * If called within a parallel simulation, initialization shall be finished
         * during domain decomposition with setLocalAndCollectiveIndices.
         *
         * \param[in] number_of_atoms The total number of atoms of this atom set.
         * \param[in] index An global integer index of the atoms in the index set.
         * \param[in] bParallel True if simulation is run in parallel
         */
        explicit LocalAtomSetData(const int number_of_atoms, const int *index, bool bParallel);

        /*! \brief Sets the local and collective indices from a lookup in ga2la.
         *
         * This makes only sense when atoms are indeed distributed over nodes,
         * otherwise local indices equal global indices and are set in init.
         *
         * \param[in] ga2la lookup table that reports if an atom is local.
         */
        void setLocalAndCollectiveIndices(const gmx_ga2la_t  *ga2la);


        /*! \brief Global indices of the atoms in this set. */
        std::vector<int> global_index_;
        /*! \brief Maps indices on node (0..num_atoms_local_) to global atom indicices. */
        std::vector<int> collective_index_;
        /*! \brief Local indices of the atoms.
         * Access,e.g., the i-th local atom coordinate of this set by x[local_index_[i]].
         * Constructed and updated every domain-decomposition step.
         */
        std::vector<int> local_index_;
};

} // namespace internal

} // namespace gmx




#endif
