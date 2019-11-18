/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2019, by the GROMACS development team, led by
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
/*! \internal \file
 *
 * \brief Declares the AtomDistribution struct.
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_domdec
 */
#ifndef GMX_DOMDEC_ATOMDISTRIBUTION_H
#define GMX_DOMDEC_ATOMDISTRIBUTION_H

#include <array>
#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/arrayref.h"

/*! \internal
 * \brief Distribution of atom groups over the domain (only available on the master rank)
 */
struct AtomDistribution
{
    /*! \internal
     * \brief Collection of local group and atom counts for a domain
     */
    struct DomainAtomGroups
    {
        gmx::ArrayRef<const int> atomGroups; /**< List of our atom groups */
        int                      numAtoms;   /**< Our number of local atoms */
    };

    /*! \brief Constructor */
    AtomDistribution(const ivec numCells, int numAtomGroups, int numAtoms);

    std::vector<DomainAtomGroups> domainGroups; /**< Group and atom division over ranks/domains */
    std::vector<int>              atomGroups; /**< The atom group division of the whole system, pointed into by counts[].atomGroups */

    /* Temporary buffers, stored permanently here to avoid reallocation */
    std::array<std::vector<real>, DIM> cellSizesBuffer; /**< Cell boundaries, sizes: num_cells_in_dim + 1 */
    std::vector<int>       intBuffer;  /**< Buffer for communicating cg and atom counts */
    std::vector<gmx::RVec> rvecBuffer; /**< Buffer for state scattering and gathering */
};

/*! \brief Returns state scatter/gather buffer element counts and displacements
 *
 * NOTE: Should only be called with a pointer to a valid ma struct
 *       (only available on the master rank).
 */
void get_commbuffer_counts(AtomDistribution* ma, int** counts, int** disps);

#endif
