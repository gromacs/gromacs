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
#ifndef GMX_MDTYPES_COMMREC_H
#define GMX_MDTYPES_COMMREC_H

#include <cstddef>

#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/gmxmpi.h"

struct gmx_domdec_t;

#define DUTY_PP (1U << 0U)
#define DUTY_PME (1U << 1U)

//! Whether the current DD role is main or slave
enum class DDRole
{
    Main,
    Agent
};

//! Whether one or more ranks are used
enum class NumRanks
{
    Single,
    Multiple
};

typedef struct
{
    int      bUse;
    MPI_Comm comm_intra;
    int      rank_intra;
    MPI_Comm comm_inter;

} gmx_nodecomm_t;

struct t_commrec
{
    /* The nodeids in one sim are numbered sequentially from 0.
     * All communication within some simulation should happen
     * in mpi_comm_mysim, or its subset mpi_comm_mygroup.
     */
    //! The rank-id in mpi_comm_mysim;
    int sim_nodeid;
    //! The number of ranks in mpi_comm_mysim
    int nnodes;
    //! The number of separate PME ranks, 0 when no separate PME ranks are used
    int npmenodes;

    //! The rank-id in mpi_comm_mygroup;
    int nodeid;

    /* MPI communicators within a single simulation
     * Note: other parts of the code may further subset these communicators.
     */
    MPI_Comm mpi_comm_mysim;   /* communicator including all ranks of
                                  a single simulation */
    MPI_Comm mpi_comm_mygroup; /* subset of mpi_comm_mysim including only
                                  the ranks in the same group (PP or PME) */
    //! The number of ranks in mpi_comm_mygroup
    int sizeOfMyGroupCommunicator;

    //! The communicator used before DD was initialized
    MPI_Comm mpiDefaultCommunicator;
    int      sizeOfDefaultCommunicator;
    int      rankInDefaultCommunicator;

    gmx_nodecomm_t nc;

    /* Handle to domain decomposition manager, owned elsewhere in mdrunner. */
    gmx_domdec_t* dd;

    /* The duties of this node, see the DUTY_ defines above.
     * This should be read through thisRankHasDuty() or getThisRankDuties().
     */
    int duty;
};

/*! \brief
 * Returns the rank's duty, and asserts that it has been initialized.
 */
inline int getThisRankDuties(const t_commrec* cr)
{
    GMX_ASSERT(cr, "Invalid commrec pointer");
    GMX_ASSERT(cr->duty != 0, "Commrec duty was not initialized!");
    return cr->duty;
}

/*! \brief
 * A convenience getter for the commrec duty assignment;
 * asserts that duty is actually valid (have been initialized).
 *
 * \param[in] cr    Communication structure pointer
 * \param[in] duty  A single duty's corresponding DUTY_ flag. Combinations are not supported.
 *
 * \returns Whether this duty is assigned to this rank.
 */
inline bool thisRankHasDuty(const t_commrec* cr, int duty)
{
    GMX_ASSERT((duty == DUTY_PME) || (duty == DUTY_PP), "Invalid duty type");
    return (getThisRankDuties(cr) & duty) != 0;
}

/*! \brief True if this is a simulation with more than 1 rank
 *
 * In particular, this is true for multi-rank runs with TPI and NM, because
 * they use a decomposition that is not the domain decomposition used by
 * other simulation types. */
#define PAR(cr) ((cr)->sizeOfDefaultCommunicator > 1)

//! True of this is the main node
#define MAIN(cr) (((cr)->rankInDefaultCommunicator == 0) || !PAR(cr))

// Note that currently, main is always PP main, so this is equivalent to MAIN(cr)
//! True if this is the particle-particle main
#define SIMMAIN(cr) ((MAIN(cr) && thisRankHasDuty((cr), DUTY_PP)) || !PAR(cr))

//! The node id for this rank
#define RANK(cr, nodeid) (nodeid)

//! The node id for the main
#define MAINRANK(cr) (0)

/*! \brief Returns whether the domain decomposition machinery is active and reorders atoms
 *
 * This tells whether atoms are reordered at pair search steps. When the return value
 * is true, atoms are not in the order of the input and mtop.
 *
 * Note that when the return value is true, there are not necessarily
 * multiple domains. The domain decomposition machinery is also active and
 * reorders the atoms also with a single MPI rank, or 1 PP and 1 PME rank,
 * with most integrators. Only a few special non-integrator "integrators"
 * do not (yet) support the domain decomposition machinery and therefore
 * this function is still needed.
 */
static bool inline haveDDAtomOrdering(const t_commrec& cr)
{
    return cr.dd != nullptr;
}

/*! \brief Returns whether we have actual domain decomposition for the particle-particle interactions
 *
 * Will return false when we use 1 rank for PP and 1 for PME
 */
static bool inline havePPDomainDecomposition(const t_commrec* cr)
{
    /* NOTE: It would be better to use cr->dd->nnodes, but we do not want
     *       to pull in a dependency on domdec.h into this file.
     */
    GMX_ASSERT(cr != nullptr, "Invalid call of havePPDomainDecomposition before commrec is made");
    GMX_ASSERT(cr->npmenodes >= 0,
               "Invalid call of havePPDomainDecomposition before MPMD automated decomposition was "
               "chosen.");
    return (cr->dd != nullptr && cr->nnodes - cr->npmenodes > 1);
}

#endif
