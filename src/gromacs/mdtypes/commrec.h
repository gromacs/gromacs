/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2017,2018,2019, by the GROMACS development team, led by
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
#ifndef GMX_MDTYPES_COMMREC_H
#define GMX_MDTYPES_COMMREC_H

#include <stddef.h>

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/gmxmpi.h"

struct mpi_in_place_buf_t;
struct gmx_domdec_t;

#define DUTY_PP (1U << 0U)
#define DUTY_PME (1U << 1U)

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
    int sim_nodeid, nnodes, npmenodes;

    /* thread numbers: */
    /* Not used yet: int threadid, nthreads; */
    /* The nodeid in the PP/PME, PP or PME group */
    int nodeid;

    /* MPI communicators within a single simulation
     * Note: other parts of the code may further subset these communicators.
     */
    MPI_Comm mpi_comm_mysim;   /* communicator including all ranks of
                                  a single simulation */
    MPI_Comm mpi_comm_mygroup; /* subset of mpi_comm_mysim including only
                                  the ranks in the same group (PP or PME) */

    gmx_nodecomm_t nc;

    /* For domain decomposition */
    gmx_domdec_t* dd;

    /* The duties of this node, see the DUTY_ defines above.
     * This should be read through thisRankHasDuty() or getThisRankDuties().
     */
    int duty;

    /* these buffers are used as destination buffers if MPI_IN_PLACE isn't
       supported.*/
    mpi_in_place_buf_t* mpb;
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
#define PAR(cr) ((cr)->nnodes > 1)

//! True of this is the master node
#define MASTER(cr) (((cr)->nodeid == 0) || !PAR(cr))

//! True if this is the particle-particle master
#define SIMMASTER(cr) ((MASTER(cr) && thisRankHasDuty((cr), DUTY_PP)) || !PAR(cr))

//! The node id for this rank
#define RANK(cr, nodeid) (nodeid)

//! The node id for the master
#define MASTERRANK(cr) (0)

/*! \brief Do we decompose the work of this simulation?
 *
 * True if this simulation uses more than one PP rank, or if this simulation
 * uses at least one PME-only rank.
 *
 * PAR(cr) is true if this is true, but the converse does not apply (see docs
 * of PAR(cr)).
 *
 * This is true if havePPDomainDecomposition is true, but the converse does not
 * apply (see docs of havePpDomainDecomposition()).
 *
 * \todo As part of Redmine #2395, replace calls to this with
 * havePPDomainDecomposition or a call of some other/new function, as
 * appropriate to each case. Then eliminate this macro. */
#define DOMAINDECOMP(cr) (((cr)->dd != nullptr) && PAR(cr))

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
