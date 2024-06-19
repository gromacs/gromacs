/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2005- The GROMACS Authors
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
/*! \internal \file
 *
 * \brief This file declares functions for domdec to use
 * while managing inter-atomic constraints.
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_domdec
 */

#ifndef GMX_DOMDEC_DOMDEC_CONSTRAINTS_H
#define GMX_DOMDEC_DOMDEC_CONSTRAINTS_H

#include <cstdint>

#include <memory>
#include <vector>

#include "gromacs/domdec/hashedmap.h"
#include "gromacs/utility/arrayref.h"

namespace gmx
{
class Constraints;
}

struct gmx_domdec_t;
struct gmx_mtop_t;
struct InteractionList;

/*! \brief Struct used during constraint setup with domain decomposition */
struct gmx_domdec_constraints_t
{
    //! @cond Doxygen_Suppress
    std::vector<int> molb_con_offset; /**< Offset in the constraint array for each molblock */
    std::vector<int> molb_ncon_mol; /**< The number of constraints per molecule for each molblock */

    int ncon; /**< The fully local and conneced constraints */
    /* The global constraint number, only required for clearing gc_req */
    std::vector<int> con_gl;     /**< Global constraint indices for local constraints */
    std::vector<int> con_nlocat; /**< Number of local atoms (2/1/0) for each constraint */

    std::vector<bool> gc_req; /**< Boolean that tells if a global constraint index has been requested; note: size global #constraints */

    /* Hash table for keeping track of requests */
    std::unique_ptr<gmx::HashedMap<int>> ga2la; /**< Global to local communicated constraint atom only index */

    /* Multi-threading stuff */
    int                          nthread; /**< Number of threads used for DD constraint setup */
    std::vector<InteractionList> ils;     /**< Constraint ilist working arrays, size \p nthread */

    /* Buffers for requesting atoms */
    std::vector<std::vector<int>> requestedGlobalAtomIndices; /**< Buffers for requesting global atom indices, one per thread */

    //! @endcond
};

/*! \brief Clears the local indices for the constraint communication setup */
void dd_clear_local_constraint_indices(gmx_domdec_t* dd);

/*! \brief Sets up communication and atom indices for all local+connected constraints */
int dd_make_local_constraints(struct gmx_domdec_t*           dd,
                              int                            at_start,
                              const struct gmx_mtop_t&       mtop,
                              gmx::ArrayRef<const int32_t>   atomInfo,
                              gmx::Constraints*              constr,
                              int                            nrec,
                              gmx::ArrayRef<InteractionList> il_local);

/*! \brief Initializes the data structures for constraint communication */
void init_domdec_constraints(gmx_domdec_t* dd, const gmx_mtop_t& mtop);

#endif
