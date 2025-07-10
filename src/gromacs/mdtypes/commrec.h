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

#include <memory>

#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/mpicomm.h"

struct gmx_domdec_t;

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

struct t_commrec
{
    //! Constructs a valid object with all communicators set to \p mpiComm
    t_commrec(const gmx::MpiComm& mpiComm);

    ~t_commrec();

    //! Returns whether this is the main rank in the simulation, i.e. the rank that does IO
    bool isSimulationMainRank() const { return (commMySim.rank() == 0); }

    //! Transfers the ownership of \p ddUniquePtr and sets \p dd
    void setDD(std::unique_ptr<gmx_domdec_t>&& ddUniquePtr);

    //! Destroys the dd object
    void destroyDD();

    //! Whether to use NVSHMEM
    bool useNvshmem = false;
    /* The nodeids in one sim are numbered sequentially from 0.
     * All communication within some simulation should happen
     * in mpi_comm_mysim, or its subset mpi_comm_mygroup.
     */
    //! The number of separate PME ranks, 0 when no separate PME ranks are used
    int npmenodes = 0;

    /* MPI communicators within a single simulation
     * Note: other parts of the code may further subset these communicators.
     */
    /* Communicator including all ranks of a single simulation */
    gmx::MpiComm commMySim;
    /* Subset of commMySim including only the ranks in the same group (PP or PME) */
    gmx::MpiComm commMyGroup;

    //! The communicator used before DD was initialized
    gmx::MpiComm mpiDefaultCommunicator;

private:
    //! Storage for the domain decomposition data
    std::unique_ptr<gmx_domdec_t> ddUniquePtr_;

public:
    //! C-pointer to ddUniquePtr (should be replaced by a getter)
    gmx_domdec_t* dd = nullptr;
};

/*! \brief True if this is a simulation with more than 1 rank
 *
 * In particular, this is true for multi-rank runs with TPI and NM, because
 * they use a decomposition that is not the domain decomposition used by
 * other simulation types. */
#define PAR(cr) ((cr)->mpiDefaultCommunicator.size() > 1)

//! True of this is the main node
#define MAIN(cr) ((cr)->mpiDefaultCommunicator.rank() == 0)

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

#endif
