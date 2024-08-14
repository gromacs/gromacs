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

#include "gmxpre.h"

#include "domdec.h"

#include "config.h"

#include <cassert>
#include <cinttypes>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <array>
#include <filesystem>
#include <memory>
#include <string>

#include "gromacs/domdec/builder.h"
#include "gromacs/domdec/collect.h"
#include "gromacs/domdec/computemultibodycutoffs.h"
#include "gromacs/domdec/dlb.h"
#include "gromacs/domdec/dlbtiming.h"
#include "gromacs/domdec/domdec_network.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/domdec/ga2la.h"
#include "gromacs/domdec/gpuhaloexchange.h"
#include "gromacs/domdec/hashedmap.h"
#include "gromacs/domdec/localtopologychecker.h"
#include "gromacs/domdec/options.h"
#include "gromacs/domdec/partition.h"
#include "gromacs/domdec/reversetopology.h"
#include "gromacs/domdec/utility.h"
#include "gromacs/ewald/pme.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/gpu_utils/device_stream_manager.h"
#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/hardware/hw_info.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/calc_verletbuf.h"
#include "gromacs/mdlib/constr.h"
#include "gromacs/mdlib/constraintrange.h"
#include "gromacs/mdlib/updategroups.h"
#include "gromacs/mdlib/updategroupscog.h"
#include "gromacs/mdlib/vcm.h"
#include "gromacs/mdlib/vsite.h"
#include "gromacs/mdrun/mdmodules.h"
#include "gromacs/mdrunutility/mdmodulesnotifiers.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/mdrunoptions.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pulling/pull.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/mtop_atomloops.h"
#include "gromacs/topology/mtop_lookup.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/basenetwork.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/fixedcapacityvector.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/iserializer.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/range.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/strconvert.h"
#include "gromacs/utility/stringstream.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textwriter.h"

#include "atomdistribution.h"
#include "box.h"
#include "cellsizes.h"
#include "distribute.h"
#include "domdec_constraints.h"
#include "domdec_internal.h"
#include "domdec_setup.h"
#include "domdec_specatomcomm.h"
#include "domdec_vsite.h"
#include "redistribute.h"
#include "utility.h"

namespace gmx
{
class LocalAtomSetManager;
class ObservablesReducerBuilder;
} // namespace gmx

// TODO remove this when moving domdec into gmx namespace
using gmx::ArrayRef;
using gmx::DdRankOrder;
using gmx::DlbOption;
using gmx::DomdecOptions;
using gmx::RangePartitioning;

static const char* enumValueToString(DlbState enumValue)
{
    static constexpr gmx::EnumerationArray<DlbState, const char*> dlbStateNames = {
        "off", "off", "auto", "locked", "on", "on"
    };
    return dlbStateNames[enumValue];
}

//! The minimum step interval for DD algorithms requiring global communication
static constexpr int sc_minimumGCStepInterval = 100;

static void ddindex2xyz(const ivec nc, int ind, ivec xyz)
{
    xyz[XX] = ind / (nc[YY] * nc[ZZ]);
    xyz[YY] = (ind / nc[ZZ]) % nc[YY];
    xyz[ZZ] = ind % nc[ZZ];
}

//! Returns the MPI rank for the PP domain corresponding to \p coord
static int ddRankFromDDCoord(const gmx_domdec_t& dd, const gmx::IVec& coord)
{
    int rank = -1;

    const CartesianRankSetup& cartSetup = dd.comm->cartesianRankSetup;
    const int                 ddindex   = dd_index(dd.numCells, coord);
    if (cartSetup.bCartesianPP_PME)
    {
        rank = cartSetup.ddindex2ddnodeid[ddindex];
    }
    else if (cartSetup.bCartesianPP)
    {
#if GMX_MPI
        gmx::IVec tmp = coord;
        MPI_Cart_rank(dd.mpi_comm_all, static_cast<int*>(&tmp[0]), &rank);
#endif
    }
    else
    {
        rank = ddindex;
    }

    return rank;
}

int ddglatnr(const gmx_domdec_t* dd, int i)
{
    int atnr = 0;

    if (dd == nullptr)
    {
        atnr = i + 1;
    }
    else
    {
        if (i >= dd->comm->atomRanges.numAtomsTotal())
        {
            gmx_fatal(FARGS,
                      "glatnr called with %d, which is larger than the local number of atoms (%d)",
                      i,
                      dd->comm->atomRanges.numAtomsTotal());
        }
        atnr = dd->globalAtomIndices[i] + 1;
    }

    return atnr;
}

void dd_store_state(const gmx_domdec_t& dd, t_state* state)
{
    if (state->ddp_count != dd.ddp_count)
    {
        gmx_incons("The MD state does not match the domain decomposition state");
    }

    state->cg_gl.resize(dd.numHomeAtoms);
    for (int i = 0; i < dd.numHomeAtoms; i++)
    {
        state->cg_gl[i] = dd.globalAtomIndices[i];
    }

    state->ddp_count_cg_gl = dd.ddp_count;
}

const gmx::DomdecZones& getDomdecZones(const gmx_domdec_t& dd)
{
    return dd.zones;
}

int dd_numAtomsZones(const gmx_domdec_t& dd)
{
    return dd.comm->atomRanges.end(DDAtomRanges::Type::Zones);
}

int dd_numHomeAtoms(const gmx_domdec_t& dd)
{
    return dd.comm->atomRanges.numHomeAtoms();
}

int dd_natoms_mdatoms(const gmx_domdec_t& dd)
{
    /* We currently set mdatoms entries for all atoms:
     * local + non-local + communicated for vsite + constraints
     */

    return dd.comm->atomRanges.numAtomsTotal();
}

int dd_natoms_vsite(const gmx_domdec_t& dd)
{
    return dd.comm->atomRanges.end(DDAtomRanges::Type::Vsites);
}

void dd_get_constraint_range(const gmx_domdec_t& dd, int* at_start, int* at_end)
{
    *at_start = dd.comm->atomRanges.start(DDAtomRanges::Type::Constraints);
    *at_end   = dd.comm->atomRanges.end(DDAtomRanges::Type::Constraints);
}

void dd_move_x(gmx_domdec_t* dd, const matrix box, gmx::ArrayRef<gmx::RVec> x, gmx_wallcycle* wcycle)
{
    wallcycle_start(wcycle, WallCycleCounter::MoveX);

    rvec shift = { 0, 0, 0 };

    gmx_domdec_comm_t* comm = dd->comm.get();

    int nzone   = 1;
    int nat_tot = comm->atomRanges.numHomeAtoms();
    for (int d = 0; d < dd->ndim; d++)
    {
        const bool bPBC   = (dd->ci[dd->dim[d]] == 0);
        const bool bScrew = (bPBC && dd->unitCellInfo.haveScrewPBC && dd->dim[d] == XX);
        if (bPBC)
        {
            copy_rvec(box[dd->dim[d]], shift);
        }
        gmx_domdec_comm_dim_t* cd = &comm->cd[d];
        for (const gmx_domdec_ind_t& ind : cd->ind)
        {
            DDBufferAccess<gmx::RVec> sendBufferAccess(comm->rvecBuffer, ind.nsend[nzone + 1]);
            gmx::ArrayRef<gmx::RVec>& sendBuffer = sendBufferAccess.buffer;
            int                       n          = 0;
            if (!bPBC)
            {
                for (int j : ind.index)
                {
                    sendBuffer[n] = x[j];
                    n++;
                }
            }
            else if (!bScrew)
            {
                for (int j : ind.index)
                {
                    /* We need to shift the coordinates */
                    for (int d = 0; d < DIM; d++)
                    {
                        sendBuffer[n][d] = x[j][d] + shift[d];
                    }
                    n++;
                }
            }
            else
            {
                for (int j : ind.index)
                {
                    /* Shift x */
                    sendBuffer[n][XX] = x[j][XX] + shift[XX];
                    /* Rotate y and z.
                     * This operation requires a special shift force
                     * treatment, which is performed in calc_vir.
                     */
                    sendBuffer[n][YY] = box[YY][YY] - x[j][YY];
                    sendBuffer[n][ZZ] = box[ZZ][ZZ] - x[j][ZZ];
                    n++;
                }
            }

            DDBufferAccess<gmx::RVec> receiveBufferAccess(
                    comm->rvecBuffer2, cd->receiveInPlace ? 0 : ind.nrecv[nzone + 1]);

            gmx::ArrayRef<gmx::RVec> receiveBuffer;
            if (cd->receiveInPlace)
            {
                receiveBuffer = gmx::arrayRefFromArray(x.data() + nat_tot, ind.nrecv[nzone + 1]);
            }
            else
            {
                receiveBuffer = receiveBufferAccess.buffer;
            }
            /* Send and receive the coordinates */
            ddSendrecv(dd, d, dddirBackward, sendBuffer, receiveBuffer);

            if (!cd->receiveInPlace)
            {
                int j = 0;
                for (int zone = 0; zone < nzone; zone++)
                {
                    for (int i = ind.cell2at0[zone]; i < ind.cell2at1[zone]; i++)
                    {
                        x[i] = receiveBuffer[j++];
                    }
                }
            }
            nat_tot += ind.nrecv[nzone + 1];
        }
        nzone += nzone;
    }

    wallcycle_stop(wcycle, WallCycleCounter::MoveX);
}

void dd_move_f(gmx_domdec_t* dd, gmx::ForceWithShiftForces* forceWithShiftForces, gmx_wallcycle* wcycle)
{
    wallcycle_start(wcycle, WallCycleCounter::MoveF);

    gmx::ArrayRef<gmx::RVec> f      = forceWithShiftForces->force();
    gmx::ArrayRef<gmx::RVec> fshift = forceWithShiftForces->shiftForces();

    gmx_domdec_comm_t& comm    = *dd->comm;
    int                nzone   = dd->zones.numZones() / 2;
    int                nat_tot = comm.atomRanges.end(DDAtomRanges::Type::Zones);
    for (int d = dd->ndim - 1; d >= 0; d--)
    {
        /* Only forces in domains near the PBC boundaries need to
           consider PBC in the treatment of fshift */
        const bool shiftForcesNeedPbc =
                (forceWithShiftForces->computeVirial() && dd->ci[dd->dim[d]] == 0);
        const bool applyScrewPbc = (dd->unitCellInfo.haveScrewPBC && dd->dim[d] == XX);
        /* Determine which shift vector we need */
        ivec vis        = { 0, 0, 0 };
        vis[dd->dim[d]] = 1;
        const int is    = gmx::ivecToShiftIndex(vis);

        /* Loop over the pulses */
        const gmx_domdec_comm_dim_t& cd = comm.cd[d];
        for (int p = cd.numPulses() - 1; p >= 0; p--)
        {
            const gmx_domdec_ind_t&   ind = cd.ind[p];
            DDBufferAccess<gmx::RVec> receiveBufferAccess(comm.rvecBuffer, ind.nsend[nzone + 1]);
            gmx::ArrayRef<gmx::RVec>& receiveBuffer = receiveBufferAccess.buffer;

            nat_tot -= ind.nrecv[nzone + 1];

            DDBufferAccess<gmx::RVec> sendBufferAccess(
                    comm.rvecBuffer2, cd.receiveInPlace ? 0 : ind.nrecv[nzone + 1]);

            gmx::ArrayRef<gmx::RVec> sendBuffer;
            if (cd.receiveInPlace)
            {
                sendBuffer = gmx::arrayRefFromArray(f.data() + nat_tot, ind.nrecv[nzone + 1]);
            }
            else
            {
                sendBuffer = sendBufferAccess.buffer;
                int j      = 0;
                for (int zone = 0; zone < nzone; zone++)
                {
                    for (int i = ind.cell2at0[zone]; i < ind.cell2at1[zone]; i++)
                    {
                        sendBuffer[j++] = f[i];
                    }
                }
            }
            /* Communicate the forces */
            ddSendrecv(dd, d, dddirForward, sendBuffer, receiveBuffer);
            /* Add the received forces */
            int n = 0;
            if (!applyScrewPbc && !shiftForcesNeedPbc)
            {
                for (int j : ind.index)
                {
                    for (int d = 0; d < DIM; d++)
                    {
                        f[j][d] += receiveBuffer[n][d];
                    }
                    n++;
                }
            }
            else if (!applyScrewPbc)
            {
                for (int j : ind.index)
                {
                    for (int d = 0; d < DIM; d++)
                    {
                        f[j][d] += receiveBuffer[n][d];
                    }
                    /* Add this force to the shift force */
                    for (int d = 0; d < DIM; d++)
                    {
                        fshift[is][d] += receiveBuffer[n][d];
                    }
                    n++;
                }
            }
            else
            {
                for (int j : ind.index)
                {
                    /* Rotate the force */
                    f[j][XX] += receiveBuffer[n][XX];
                    f[j][YY] -= receiveBuffer[n][YY];
                    f[j][ZZ] -= receiveBuffer[n][ZZ];
                    if (shiftForcesNeedPbc)
                    {
                        /* Add this force to the shift force */
                        for (int d = 0; d < DIM; d++)
                        {
                            fshift[is][d] += receiveBuffer[n][d];
                        }
                    }
                    n++;
                }
            }
        }
        nzone /= 2;
    }
    wallcycle_stop(wcycle, WallCycleCounter::MoveF);
}

real dd_cutoff_multibody(const gmx_domdec_t* dd)
{
    const gmx_domdec_comm_t& comm       = *dd->comm;
    const DDSystemInfo&      systemInfo = comm.systemInfo;

    real r = -1;
    if (systemInfo.haveInterDomainMultiBodyBondeds)
    {
        if (comm.cutoff_mbody > 0)
        {
            r = comm.cutoff_mbody;
        }
        else
        {
            /* cutoff_mbody=0 means we do not have DLB */
            r = comm.cellsize_min[dd->dim[0]];
            for (int di = 1; di < dd->ndim; di++)
            {
                r = std::min(r, comm.cellsize_min[dd->dim[di]]);
            }
            if (comm.systemInfo.filterBondedCommunication)
            {
                r = std::max(r, comm.cutoff_mbody);
            }
            else
            {
                r = std::min(r, systemInfo.cutoff);
            }
        }
    }

    return r;
}

real dd_cutoff_twobody(const gmx_domdec_t* dd)
{
    const real r_mb = dd_cutoff_multibody(dd);

    return std::max(dd->comm->systemInfo.cutoff, r_mb);
}


static void dd_cart_coord2pmecoord(const DDRankSetup&        ddRankSetup,
                                   const CartesianRankSetup& cartSetup,
                                   const ivec                coord,
                                   ivec                      coord_pme)
{
    const int nc   = ddRankSetup.numPPCells[cartSetup.cartpmedim];
    const int ntot = cartSetup.ntot[cartSetup.cartpmedim];
    copy_ivec(coord, coord_pme);
    coord_pme[cartSetup.cartpmedim] =
            nc + (coord[cartSetup.cartpmedim] * (ntot - nc) + (ntot - nc) / 2) / nc;
}

/* Returns the PME rank index in 0...npmenodes-1 for the PP cell with index ddCellIndex */
static int ddindex2pmeindex(const DDRankSetup& ddRankSetup, const int ddCellIndex)
{
    const int npp  = ddRankSetup.numPPRanks;
    const int npme = ddRankSetup.numRanksDoingPme;

    /* Here we assign a PME node to communicate with this DD node
     * by assuming that the major index of both is x.
     * We add npme/2 to obtain an even distribution.
     */
    return (ddCellIndex * npme + npme / 2) / npp;
}

static std::vector<int> dd_interleaved_pme_ranks(const DDRankSetup& ddRankSetup)
{
    std::vector<int> pmeRanks(ddRankSetup.numRanksDoingPme);

    int n = 0;
    for (int i = 0; i < ddRankSetup.numPPRanks; i++)
    {
        const int p0 = ddindex2pmeindex(ddRankSetup, i);
        const int p1 = ddindex2pmeindex(ddRankSetup, i + 1);
        if (i + 1 == ddRankSetup.numPPRanks || p1 > p0)
        {
            if (debug)
            {
                fprintf(debug, "pme_rank[%d] = %d\n", n, i + 1 + n);
            }
            pmeRanks[n] = i + 1 + n;
            n++;
        }
    }

    return pmeRanks;
}

static int gmx_ddcoord2pmeindex(const gmx_domdec_t& dd, int x, int y, int z)
{
    ivec coords;

    coords[XX]     = x;
    coords[YY]     = y;
    coords[ZZ]     = z;
    const int slab = ddindex2pmeindex(dd.comm->ddRankSetup, dd_index(dd.numCells, coords));

    return slab;
}

static int ddcoord2simnodeid(const t_commrec* cr, int x, int y, int z)
{
    const CartesianRankSetup& cartSetup = cr->dd->comm->cartesianRankSetup;
    ivec                      coords    = { x, y, z };
    int                       nodeid    = -1;

    if (cartSetup.bCartesianPP_PME)
    {
#if GMX_MPI
        MPI_Cart_rank(cr->mpi_comm_mysim, coords, &nodeid);
#endif
    }
    else
    {
        int ddindex = dd_index(cr->dd->numCells, coords);
        if (cartSetup.bCartesianPP)
        {
            nodeid = cartSetup.ddindex2simnodeid[ddindex];
        }
        else
        {
            const DDRankSetup& rankSetup = cr->dd->comm->ddRankSetup;
            if (rankSetup.rankOrder != DdRankOrder::pp_pme && rankSetup.usePmeOnlyRanks)
            {
                nodeid = ddindex + gmx_ddcoord2pmeindex(*cr->dd, x, y, z);
            }
            else
            {
                nodeid = ddindex;
            }
        }
    }

    return nodeid;
}

static int dd_simnode2pmenode(const DDRankSetup&        ddRankSetup,
                              const CartesianRankSetup& cartSetup,
                              gmx::ArrayRef<const int>  pmeRanks,
                              const t_commrec gmx_unused* cr,
                              const int                   sim_nodeid)
{
    int pmenode = -1;

    /* This assumes a uniform x domain decomposition grid cell size */
    if (cartSetup.bCartesianPP_PME)
    {
#if GMX_MPI
        ivec coord, coord_pme;
        MPI_Cart_coords(cr->mpi_comm_mysim, sim_nodeid, DIM, coord);
        if (coord[cartSetup.cartpmedim] < ddRankSetup.numPPCells[cartSetup.cartpmedim])
        {
            /* This is a PP rank */
            dd_cart_coord2pmecoord(ddRankSetup, cartSetup, coord, coord_pme);
            MPI_Cart_rank(cr->mpi_comm_mysim, coord_pme, &pmenode);
        }
#endif
    }
    else if (cartSetup.bCartesianPP)
    {
        if (sim_nodeid < ddRankSetup.numPPRanks)
        {
            pmenode = ddRankSetup.numPPRanks + ddindex2pmeindex(ddRankSetup, sim_nodeid);
        }
    }
    else
    {
        /* This assumes DD cells with identical x coordinates
         * are numbered sequentially.
         */
        if (pmeRanks.empty())
        {
            if (sim_nodeid < ddRankSetup.numPPRanks)
            {
                /* The DD index equals the nodeid */
                pmenode = ddRankSetup.numPPRanks + ddindex2pmeindex(ddRankSetup, sim_nodeid);
            }
        }
        else
        {
            int i = 0;
            while (sim_nodeid > pmeRanks[i])
            {
                i++;
            }
            if (sim_nodeid < pmeRanks[i])
            {
                pmenode = pmeRanks[i];
            }
        }
    }

    return pmenode;
}

NumPmeDomains getNumPmeDomains(const gmx_domdec_t* dd)
{
    if (dd != nullptr)
    {
        return { dd->comm->ddRankSetup.npmenodes_x, dd->comm->ddRankSetup.npmenodes_y };
    }
    else
    {
        return { 1, 1 };
    }
}

std::vector<int> get_pme_ddranks(const t_commrec* cr, const int pmenodeid)
{
    const DDRankSetup&        ddRankSetup = cr->dd->comm->ddRankSetup;
    const CartesianRankSetup& cartSetup   = cr->dd->comm->cartesianRankSetup;
    GMX_RELEASE_ASSERT(ddRankSetup.usePmeOnlyRanks,
                       "This function should only be called when PME-only ranks are in use");
    std::vector<int> ddranks;
    ddranks.reserve(gmx::divideRoundUp(ddRankSetup.numPPRanks, ddRankSetup.numRanksDoingPme));

    for (int x = 0; x < ddRankSetup.numPPCells[XX]; x++)
    {
        for (int y = 0; y < ddRankSetup.numPPCells[YY]; y++)
        {
            for (int z = 0; z < ddRankSetup.numPPCells[ZZ]; z++)
            {
                if (cartSetup.bCartesianPP_PME)
                {
                    ivec coord = { x, y, z };
                    ivec coord_pme;
                    dd_cart_coord2pmecoord(ddRankSetup, cartSetup, coord, coord_pme);
                    if (cr->dd->ci[XX] == coord_pme[XX] && cr->dd->ci[YY] == coord_pme[YY]
                        && cr->dd->ci[ZZ] == coord_pme[ZZ])
                    {
                        ddranks.push_back(ddcoord2simnodeid(cr, x, y, z));
                    }
                }
                else
                {
                    /* The slab corresponds to the nodeid in the PME group */
                    if (gmx_ddcoord2pmeindex(*cr->dd, x, y, z) == pmenodeid)
                    {
                        ddranks.push_back(ddcoord2simnodeid(cr, x, y, z));
                    }
                }
            }
        }
    }
    return ddranks;
}

static gmx_bool receive_vir_ener(const gmx_domdec_t* dd, gmx::ArrayRef<const int> pmeRanks, const t_commrec* cr)
{
    bool bReceive = true;

    const DDRankSetup& ddRankSetup = dd->comm->ddRankSetup;
    if (ddRankSetup.usePmeOnlyRanks)
    {
        const CartesianRankSetup& cartSetup = dd->comm->cartesianRankSetup;
        if (cartSetup.bCartesianPP_PME)
        {
#if GMX_MPI
            int  pmenode = dd_simnode2pmenode(ddRankSetup, cartSetup, pmeRanks, cr, cr->sim_nodeid);
            ivec coords;
            MPI_Cart_coords(cr->mpi_comm_mysim, cr->sim_nodeid, DIM, coords);
            coords[cartSetup.cartpmedim]++;
            if (coords[cartSetup.cartpmedim] < dd->numCells[cartSetup.cartpmedim])
            {
                int rank = 0;
                MPI_Cart_rank(cr->mpi_comm_mysim, coords, &rank);
                if (dd_simnode2pmenode(ddRankSetup, cartSetup, pmeRanks, cr, rank) == pmenode)
                {
                    /* This is not the last PP node for pmenode */
                    bReceive = FALSE;
                }
            }
#else
            GMX_RELEASE_ASSERT(
                    false,
                    "Without MPI we should not have Cartesian PP-PME with #PMEnodes < #DDnodes");
#endif
        }
        else
        {
            int pmenode = dd_simnode2pmenode(ddRankSetup, cartSetup, pmeRanks, cr, cr->sim_nodeid);
            if (cr->sim_nodeid + 1 < cr->nnodes
                && dd_simnode2pmenode(ddRankSetup, cartSetup, pmeRanks, cr, cr->sim_nodeid + 1) == pmenode)
            {
                /* This is not the last PP node for pmenode */
                bReceive = FALSE;
            }
        }
    }

    return bReceive;
}

static std::vector<real> set_slb_pme_dim_f(gmx_domdec_t* dd, int dim)
{
    gmx_domdec_comm_t* comm = dd->comm.get();

    std::vector<real> dim_f(dd->numCells[dim] + 1);
    dim_f[0] = 0;
    for (int i = 1; i < dd->numCells[dim]; i++)
    {
        if (!comm->slb_frac[dim].empty())
        {
            dim_f[i] = dim_f[i - 1] + comm->slb_frac[dim][i - 1];
        }
        else
        {
            dim_f[i] = static_cast<real>(i) / static_cast<real>(dd->numCells[dim]);
        }
    }
    dim_f[dd->numCells[dim]] = 1;
    return dim_f;
}

static void init_ddpme(gmx_domdec_t* dd, gmx_ddpme_t* ddpme, int dimind)
{
    const DDRankSetup& ddRankSetup = dd->comm->ddRankSetup;

    if (dimind == 0 && dd->dim[0] == YY && ddRankSetup.npmenodes_x == 1)
    {
        ddpme->dim = YY;
    }
    else
    {
        ddpme->dim = dimind;
    }
    ddpme->dim_match = (ddpme->dim == dd->dim[dimind]);

    ddpme->nslab = (ddpme->dim == 0 ? ddRankSetup.npmenodes_x : ddRankSetup.npmenodes_y);

    if (ddpme->nslab <= 1)
    {
        return;
    }

    const int nso = ddRankSetup.numRanksDoingPme / ddpme->nslab;
    /* Determine for each PME slab the PP location range for dimension dim */
    ddpme->pp_min.resize(ddpme->nslab);
    ddpme->pp_max.resize(ddpme->nslab);
    for (int slab = 0; slab < ddpme->nslab; slab++)
    {
        ddpme->pp_min[slab] = dd->numCells[dd->dim[dimind]] - 1;
        ddpme->pp_max[slab] = 0;
    }
    for (int i = 0; i < dd->nnodes; i++)
    {
        ivec xyz;
        ddindex2xyz(dd->numCells, i, xyz);
        /* For y only use our y/z slab.
         * This assumes that the PME x grid size matches the DD grid size.
         */
        if (dimind == 0 || xyz[XX] == dd->ci[XX])
        {
            const int pmeindex  = ddindex2pmeindex(ddRankSetup, i);
            const int slab      = (dimind == 0) ? (pmeindex / nso) : (pmeindex % ddpme->nslab);
            ddpme->pp_min[slab] = std::min(ddpme->pp_min[slab], xyz[dimind]);
            ddpme->pp_max[slab] = std::max(ddpme->pp_max[slab], xyz[dimind]);
        }
    }

    ddpme->slb_dim_f = set_slb_pme_dim_f(dd, ddpme->dim);
}

int dd_pme_maxshift_x(const gmx_domdec_t& dd)
{
    const DDRankSetup& ddRankSetup = dd.comm->ddRankSetup;

    if (ddRankSetup.ddpme[0].dim == XX)
    {
        return ddRankSetup.ddpme[0].maxshift;
    }
    else
    {
        return 0;
    }
}

int dd_pme_maxshift_y(const gmx_domdec_t& dd)
{
    const DDRankSetup& ddRankSetup = dd.comm->ddRankSetup;

    if (ddRankSetup.ddpme[0].dim == YY)
    {
        return ddRankSetup.ddpme[0].maxshift;
    }
    else if (ddRankSetup.npmedecompdim >= 2 && ddRankSetup.ddpme[1].dim == YY)
    {
        return ddRankSetup.ddpme[1].maxshift;
    }
    else
    {
        return 0;
    }
}

bool ddUsesUpdateGroups(const gmx_domdec_t& dd)
{
    return dd.comm->systemInfo.useUpdateGroups;
}

void dd_cycles_add(const gmx_domdec_t* dd, float cycles, int ddCycl)
{
    /* Note that the cycles value can be incorrect, either 0 or some
     * extremely large value, when our thread migrated to another core
     * with an unsynchronized cycle counter. If this happens less often
     * that once per nstlist steps, this will not cause issues, since
     * we later subtract the maximum value from the sum over nstlist steps.
     * A zero count will slightly lower the total, but that's a small effect.
     * Note that the main purpose of the subtraction of the maximum value
     * is to avoid throwing off the load balancing when stalls occur due
     * e.g. system activity or network congestion.
     */
    dd->comm->cycl[ddCycl] += cycles;
    dd->comm->cycl_n[ddCycl]++;
    if (cycles > dd->comm->cycl_max[ddCycl])
    {
        dd->comm->cycl_max[ddCycl] = cycles;
    }
}

#if GMX_MPI
static void make_load_communicator(gmx_domdec_t* dd, int dim_ind, ivec loc)
{
    MPI_Comm c_row = MPI_COMM_NULL;
    ivec     loc_c;
    bool     bPartOfGroup = false;

    const int dim = dd->dim[dim_ind];
    copy_ivec(loc, loc_c);
    for (int i = 0; i < dd->numCells[dim]; i++)
    {
        loc_c[dim]     = i;
        const int rank = dd_index(dd->numCells, loc_c);
        if (rank == dd->rank)
        {
            /* This process is part of the group */
            bPartOfGroup = TRUE;
        }
    }
    MPI_Comm_split(dd->mpi_comm_all, bPartOfGroup ? 0 : MPI_UNDEFINED, dd->rank, &c_row);
    if (bPartOfGroup)
    {
        dd->comm->mpi_comm_load[dim_ind] = c_row;
        if (!isDlbDisabled(dd->comm->dlbState))
        {
            DDCellsizesWithDlb& cellsizes = dd->comm->cellsizesWithDlb[dim_ind];

            if (dd->ci[dim] == dd->main_ci[dim])
            {
                /* This is the root process of this row */
                cellsizes.rowCoordinator = std::make_unique<RowCoordinator>();

                RowCoordinator& rowCoordinator = *cellsizes.rowCoordinator;
                rowCoordinator.cellFrac.resize(ddCellFractionBufferSize(dd, dim_ind));
                rowCoordinator.oldCellFrac.resize(dd->numCells[dim] + 1);
                rowCoordinator.isCellMin.resize(dd->numCells[dim]);
                if (dim_ind > 0)
                {
                    rowCoordinator.bounds.resize(dd->numCells[dim]);
                }
                rowCoordinator.buf_ncd.resize(dd->numCells[dim]);
            }
            else
            {
                /* This is not a root process, we only need to receive cell_f */
                cellsizes.fracRow.resize(ddCellFractionBufferSize(dd, dim_ind));
            }
        }
        if (dd->ci[dim] == dd->main_ci[dim])
        {
            dd->comm->load[dim_ind].load.resize(dd->numCells[dim] * DD_NLOAD_MAX);
        }
    }
}
#endif

void dd_setup_dlb_resource_sharing(const t_commrec* cr, int gpu_id)
{
#if GMX_MPI
    gmx_domdec_t* dd = cr->dd;

    if (!thisRankHasDuty(cr, DUTY_PP) || gpu_id < 0)
    {
        /* Only ranks with short-ranged tasks (currently) use GPUs.
         * If we don't have GPUs assigned, there are no resources to share.
         */
        return;
    }

    if (cr->nnodes == 1)
    {
        dd->comm->nrank_gpu_shared = 1;

        return;
    }

    const int physicalnode_id_hash = gmx_physicalnode_id_hash();

    if (debug)
    {
        fprintf(debug, "dd_setup_dd_dlb_gpu_sharing:\n");
        fprintf(debug, "DD PP rank %d physical node hash %d gpu_id %d\n", dd->rank, physicalnode_id_hash, gpu_id);
    }
    /* Split the PP communicator over the physical nodes */
    /* TODO: See if we should store this (before), as it's also used for
     * for the nodecomm summation.
     */
    // TODO PhysicalNodeCommunicator could be extended/used to handle
    // the need for per-node per-group communicators.
    MPI_Comm mpi_comm_pp_physicalnode;
    MPI_Comm_split(dd->mpi_comm_all, physicalnode_id_hash, dd->rank, &mpi_comm_pp_physicalnode);
    MPI_Comm_split(mpi_comm_pp_physicalnode, gpu_id, dd->rank, &dd->comm->mpi_comm_gpu_shared);
    MPI_Comm_free(&mpi_comm_pp_physicalnode);
    MPI_Comm_size(dd->comm->mpi_comm_gpu_shared, &dd->comm->nrank_gpu_shared);

    if (debug)
    {
        fprintf(debug, "nrank_gpu_shared %d\n", dd->comm->nrank_gpu_shared);
    }

    /* Note that some ranks could share a GPU, while others don't */

    if (dd->comm->nrank_gpu_shared == 1)
    {
        MPI_Comm_free(&dd->comm->mpi_comm_gpu_shared);
    }
#else
    GMX_UNUSED_VALUE(cr);
    GMX_UNUSED_VALUE(gpu_id);
#endif
}

static void make_load_communicators(gmx_domdec_t gmx_unused* dd)
{
#if GMX_MPI
    ivec loc;

    if (debug)
    {
        fprintf(debug, "Making load communicators\n");
    }

    dd->comm->load.resize(std::max(dd->ndim, 1));
    dd->comm->mpi_comm_load.resize(std::max(dd->ndim, 1));

    if (dd->ndim == 0)
    {
        return;
    }

    clear_ivec(loc);
    make_load_communicator(dd, 0, loc);
    if (dd->ndim > 1)
    {
        const int dim0 = dd->dim[0];
        for (int i = 0; i < dd->numCells[dim0]; i++)
        {
            loc[dim0] = i;
            make_load_communicator(dd, 1, loc);
        }
    }
    if (dd->ndim > 2)
    {
        const int dim0 = dd->dim[0];
        for (int i = 0; i < dd->numCells[dim0]; i++)
        {
            loc[dim0]      = i;
            const int dim1 = dd->dim[1];
            for (int j = 0; j < dd->numCells[dim1]; j++)
            {
                loc[dim1] = j;
                make_load_communicator(dd, 2, loc);
            }
        }
    }

    if (debug)
    {
        fprintf(debug, "Finished making load communicators\n");
    }
#endif
}

/*! \brief Sets up the relation between neighboring domains and zones */
static void setup_neighbor_relations(gmx_domdec_t* dd)
{
    GMX_ASSERT((dd->ndim >= 0) && (dd->ndim <= DIM), "Must have valid number of dimensions for DD");

    const std::array<int, 2> shifts = { 1, -1 };

    for (int d = 0; d < dd->ndim; d++)
    {
        const int dim = dd->dim[d];

        for (gmx::Index i = 0; i < gmx::ssize(shifts); i++)
        {
            gmx::IVec tmp      = dd->ci;
            tmp[dim]           = (tmp[dim] + shifts[i] + dd->numCells[dim]) % dd->numCells[dim];
            dd->neighbor[d][i] = ddRankFromDDCoord(*dd, tmp);
        }

        if (debug)
        {
            fprintf(debug,
                    "DD rank %d neighbor ranks in dir %d are + %d - %d\n",
                    dd->rank,
                    dim,
                    dd->neighbor[d][0],
                    dd->neighbor[d][1]);
        }
    }

    if (!isDlbDisabled(dd->comm->dlbState))
    {
        dd->comm->cellsizesWithDlb.resize(dd->ndim);
    }

    if (dd->comm->ddSettings.recordLoad)
    {
        make_load_communicators(dd);
    }
}

static void make_pp_communicator(const gmx::MDLogger& mdlog,
                                 gmx_domdec_t*        dd,
                                 t_commrec gmx_unused* cr,
                                 bool gmx_unused       reorder)
{
#if GMX_MPI
    gmx_domdec_comm_t*  comm      = dd->comm.get();
    CartesianRankSetup& cartSetup = comm->cartesianRankSetup;

    if (cartSetup.bCartesianPP)
    {
        /* Set up cartesian communication for the particle-particle part */
        GMX_LOG(mdlog.info)
                .appendTextFormatted("Will use a Cartesian communicator: %d x %d x %d",
                                     dd->numCells[XX],
                                     dd->numCells[YY],
                                     dd->numCells[ZZ]);

        ivec periods;
        for (int i = 0; i < DIM; i++)
        {
            periods[i] = TRUE;
        }
        MPI_Comm comm_cart = MPI_COMM_NULL;
        MPI_Cart_create(cr->mpi_comm_mygroup, DIM, dd->numCells, periods, static_cast<int>(reorder), &comm_cart);
        /* We overwrite the old communicator with the new cartesian one */
        cr->mpi_comm_mygroup = comm_cart;
    }

    dd->mpi_comm_all = cr->mpi_comm_mygroup;
    MPI_Comm_rank(dd->mpi_comm_all, &dd->rank);

    if (cartSetup.bCartesianPP_PME)
    {
        /* Since we want to use the original cartesian setup for sim,
         * and not the one after split, we need to make an index.
         */
        cartSetup.ddindex2ddnodeid.resize(dd->nnodes);
        cartSetup.ddindex2ddnodeid[dd_index(dd->numCells, dd->ci)] = dd->rank;
        gmx_sumi(dd->nnodes, cartSetup.ddindex2ddnodeid.data(), cr);
        /* Get the rank of the DD main,
         * above we made sure that the main node is a PP node.
         */
        int rank = MAIN(cr) ? dd->rank : 0;
        MPI_Allreduce(&rank, &dd->mainrank, 1, MPI_INT, MPI_SUM, dd->mpi_comm_all);
    }
    else if (cartSetup.bCartesianPP)
    {
        if (!comm->ddRankSetup.usePmeOnlyRanks)
        {
            /* The PP communicator is also
             * the communicator for this simulation
             */
            cr->mpi_comm_mysim = cr->mpi_comm_mygroup;
        }
        cr->nodeid = dd->rank;

        MPI_Cart_coords(dd->mpi_comm_all, dd->rank, DIM, dd->ci);

        /* We need to make an index to go from the coordinates
         * to the nodeid of this simulation.
         */
        cartSetup.ddindex2simnodeid.resize(dd->nnodes);
        std::vector<int> buf(dd->nnodes);
        if (thisRankHasDuty(cr, DUTY_PP))
        {
            buf[dd_index(dd->numCells, dd->ci)] = cr->sim_nodeid;
        }
        /* Communicate the ddindex to simulation nodeid index */
        MPI_Allreduce(buf.data(), cartSetup.ddindex2simnodeid.data(), dd->nnodes, MPI_INT, MPI_SUM, cr->mpi_comm_mysim);

        /* Determine the main coordinates and rank.
         * The DD main should be the same node as the main of this sim.
         */
        for (int i = 0; i < dd->nnodes; i++)
        {
            if (cartSetup.ddindex2simnodeid[i] == 0)
            {
                ddindex2xyz(dd->numCells, i, dd->main_ci);
                MPI_Cart_rank(dd->mpi_comm_all, dd->main_ci, &dd->mainrank);
            }
        }
        if (debug)
        {
            fprintf(debug, "The main rank is %d\n", dd->mainrank);
        }
    }
    else
    {
        /* No Cartesian communicators */
        /* We use the rank in dd->comm->all as DD index */
        ddindex2xyz(dd->numCells, dd->rank, dd->ci);
        /* The simulation main nodeid is 0, so the DD main rank is also 0 */
        dd->mainrank = 0;
        clear_ivec(dd->main_ci);
    }
#endif

    GMX_LOG(mdlog.info)
            .appendTextFormatted("Domain decomposition rank %d, coordinates %d %d %d\n",
                                 dd->rank,
                                 dd->ci[XX],
                                 dd->ci[YY],
                                 dd->ci[ZZ]);
    if (debug)
    {
        fprintf(debug,
                "Domain decomposition rank %d, coordinates %d %d %d\n\n",
                dd->rank,
                dd->ci[XX],
                dd->ci[YY],
                dd->ci[ZZ]);
    }
}

static void receive_ddindex2simnodeid(gmx_domdec_t* dd, t_commrec* cr)
{
#if GMX_MPI
    CartesianRankSetup& cartSetup = dd->comm->cartesianRankSetup;

    if (!cartSetup.bCartesianPP_PME && cartSetup.bCartesianPP)
    {
        cartSetup.ddindex2simnodeid.resize(dd->nnodes);
        std::vector<int> buf(dd->nnodes);
        if (thisRankHasDuty(cr, DUTY_PP))
        {
            buf[dd_index(dd->numCells, dd->ci)] = cr->sim_nodeid;
        }
        /* Communicate the ddindex to simulation nodeid index */
        MPI_Allreduce(buf.data(), cartSetup.ddindex2simnodeid.data(), dd->nnodes, MPI_INT, MPI_SUM, cr->mpi_comm_mysim);
    }
#else
    GMX_UNUSED_VALUE(dd);
    GMX_UNUSED_VALUE(cr);
#endif
}

static CartesianRankSetup split_communicator(const gmx::MDLogger& mdlog,
                                             t_commrec*           cr,
                                             const DdRankOrder    ddRankOrder,
                                             bool gmx_unused      reorder,
                                             const DDRankSetup&   ddRankSetup,
                                             ivec                 ddCellIndex,
                                             std::vector<int>*    pmeRanks)
{
    CartesianRankSetup cartSetup;

    cartSetup.bCartesianPP     = (ddRankOrder == DdRankOrder::cartesian);
    cartSetup.bCartesianPP_PME = false;

    const ivec& numDDCells = ddRankSetup.numPPCells;
    /* Initially we set ntot to the number of PP cells */
    copy_ivec(numDDCells, cartSetup.ntot);

    if (cartSetup.bCartesianPP)
    {
        const int numDDCellsTot = ddRankSetup.numPPRanks;
        bool      bDiv[DIM];
        for (int i = 1; i < DIM; i++)
        {
            bDiv[i] = ((ddRankSetup.numRanksDoingPme * numDDCells[i]) % numDDCellsTot == 0);
        }
        if (bDiv[YY] || bDiv[ZZ])
        {
            cartSetup.bCartesianPP_PME = TRUE;
            /* If we have 2D PME decomposition, which is always in x+y,
             * we stack the PME only nodes in z.
             * Otherwise we choose the direction that provides the thinnest slab
             * of PME only nodes as this will have the least effect
             * on the PP communication.
             * But for the PME communication the opposite might be better.
             */
            if (bDiv[ZZ] && (ddRankSetup.npmenodes_y > 1 || !bDiv[YY] || numDDCells[YY] > numDDCells[ZZ]))
            {
                cartSetup.cartpmedim = ZZ;
            }
            else
            {
                cartSetup.cartpmedim = YY;
            }
            cartSetup.ntot[cartSetup.cartpmedim] +=
                    (ddRankSetup.numRanksDoingPme * numDDCells[cartSetup.cartpmedim]) / numDDCellsTot;
        }
        else
        {
            GMX_LOG(mdlog.info)
                    .appendTextFormatted(
                            "Number of PME-only ranks (%d) is not a multiple of nx*ny (%d*%d) or "
                            "nx*nz (%d*%d)",
                            ddRankSetup.numRanksDoingPme,
                            numDDCells[XX],
                            numDDCells[YY],
                            numDDCells[XX],
                            numDDCells[ZZ]);
            GMX_LOG(mdlog.info)
                    .appendText("Will not use a Cartesian communicator for PP <-> PME\n");
        }
    }

    if (cartSetup.bCartesianPP_PME)
    {
#if GMX_MPI
        int  rank = 0;
        ivec periods;

        GMX_LOG(mdlog.info)
                .appendTextFormatted(
                        "Will use a Cartesian communicator for PP <-> PME: %d x %d x %d",
                        cartSetup.ntot[XX],
                        cartSetup.ntot[YY],
                        cartSetup.ntot[ZZ]);

        for (int i = 0; i < DIM; i++)
        {
            periods[i] = TRUE;
        }
        MPI_Comm comm_cart = MPI_COMM_NULL;
        MPI_Cart_create(cr->mpi_comm_mysim, DIM, cartSetup.ntot, periods, static_cast<int>(reorder), &comm_cart);
        MPI_Comm_rank(comm_cart, &rank);
        if (MAIN(cr) && rank != 0)
        {
            gmx_fatal(FARGS, "MPI rank 0 was renumbered by MPI_Cart_create, we do not allow this");
        }

        /* With this assigment we loose the link to the original communicator
         * which will usually be MPI_COMM_WORLD, unless have multisim.
         */
        cr->mpi_comm_mysim = comm_cart;
        cr->sim_nodeid     = rank;

        MPI_Cart_coords(cr->mpi_comm_mysim, cr->sim_nodeid, DIM, ddCellIndex);

        GMX_LOG(mdlog.info)
                .appendTextFormatted("Cartesian rank %d, coordinates %d %d %d\n",
                                     cr->sim_nodeid,
                                     ddCellIndex[XX],
                                     ddCellIndex[YY],
                                     ddCellIndex[ZZ]);

        if (ddCellIndex[cartSetup.cartpmedim] < numDDCells[cartSetup.cartpmedim])
        {
            cr->duty = DUTY_PP;
        }
        if (!ddRankSetup.usePmeOnlyRanks
            || ddCellIndex[cartSetup.cartpmedim] >= numDDCells[cartSetup.cartpmedim])
        {
            cr->duty = DUTY_PME;
        }

        /* Split the sim communicator into PP and PME only nodes */
        MPI_Comm_split(cr->mpi_comm_mysim,
                       getThisRankDuties(cr),
                       dd_index(cartSetup.ntot, ddCellIndex),
                       &cr->mpi_comm_mygroup);
        MPI_Comm_size(cr->mpi_comm_mygroup, &cr->sizeOfMyGroupCommunicator);
#else
        GMX_UNUSED_VALUE(ddCellIndex);
#endif
    }
    else
    {
        switch (ddRankOrder)
        {
            case DdRankOrder::pp_pme:
                GMX_LOG(mdlog.info).appendText("Order of the ranks: PP first, PME last");
                break;
            case DdRankOrder::interleave:
                /* Interleave the PP-only and PME-only ranks */
                GMX_LOG(mdlog.info).appendText("Interleaving PP and PME ranks");
                *pmeRanks = dd_interleaved_pme_ranks(ddRankSetup);
                break;
            case DdRankOrder::cartesian: break;
            default: gmx_fatal(FARGS, "Invalid ddRankOrder=%d", static_cast<int>(ddRankOrder));
        }

        if (dd_simnode2pmenode(ddRankSetup, cartSetup, *pmeRanks, cr, cr->sim_nodeid) == -1)
        {
            cr->duty = DUTY_PME;
        }
        else
        {
            cr->duty = DUTY_PP;
        }
#if GMX_MPI
        /* Split the sim communicator into PP and PME only nodes */
        MPI_Comm_split(cr->mpi_comm_mysim, getThisRankDuties(cr), cr->nodeid, &cr->mpi_comm_mygroup);
        MPI_Comm_size(cr->mpi_comm_mygroup, &cr->sizeOfMyGroupCommunicator);
        MPI_Comm_rank(cr->mpi_comm_mygroup, &cr->nodeid);
#endif
    }

    GMX_LOG(mdlog.info)
            .appendTextFormatted("This rank does only %s work.\n",
                                 thisRankHasDuty(cr, DUTY_PP) ? "particle-particle" : "PME-mesh");

    return cartSetup;
}

/*! \brief Makes the PP communicator and the PME communicator, when needed
 *
 * Returns the Cartesian rank setup.
 * Sets \p cr->mpi_comm_mygroup
 * For PP ranks, sets the DD PP cell index in \p ddCellIndex.
 * With separate PME ranks in interleaved order, set the PME ranks in \p pmeRanks.
 */
static CartesianRankSetup makeGroupCommunicators(const gmx::MDLogger& mdlog,
                                                 const DDSettings&    ddSettings,
                                                 const DdRankOrder    ddRankOrder,
                                                 const DDRankSetup&   ddRankSetup,
                                                 t_commrec*           cr,
                                                 ivec                 ddCellIndex,
                                                 std::vector<int>*    pmeRanks)
{
    CartesianRankSetup cartSetup;

    // As a default, both group and sim communicators are equal to the default communicator
    cr->mpi_comm_mygroup = cr->mpiDefaultCommunicator;
    cr->mpi_comm_mysim   = cr->mpiDefaultCommunicator;
    cr->nnodes           = cr->sizeOfDefaultCommunicator;
    cr->nodeid           = cr->rankInDefaultCommunicator;
    cr->sim_nodeid       = cr->rankInDefaultCommunicator;

    if (ddRankSetup.usePmeOnlyRanks)
    {
        /* Split the communicator into a PP and PME part */
        cartSetup = split_communicator(
                mdlog, cr, ddRankOrder, ddSettings.useCartesianReorder, ddRankSetup, ddCellIndex, pmeRanks);
    }
    else
    {
        /* All nodes do PP and PME */
        /* We do not require separate communicators */
        cartSetup.bCartesianPP     = false;
        cartSetup.bCartesianPP_PME = false;
    }

    return cartSetup;
}

/*! \brief For PP ranks, sets or makes the communicator
 *
 * For PME ranks get the rank id.
 * For PP only ranks, sets the PME-only rank.
 */
static void setupGroupCommunication(const gmx::MDLogger&     mdlog,
                                    const DDSettings&        ddSettings,
                                    gmx::ArrayRef<const int> pmeRanks,
                                    t_commrec*               cr,
                                    const int                numAtomsInSystem,
                                    gmx_domdec_t*            dd)
{
    const DDRankSetup&        ddRankSetup = dd->comm->ddRankSetup;
    const CartesianRankSetup& cartSetup   = dd->comm->cartesianRankSetup;

    if (thisRankHasDuty(cr, DUTY_PP))
    {
        if (dd->nnodes > 1)
        {
            /* Copy or make a new PP communicator */

            /* We (possibly) reordered the nodes in split_communicator,
             * so it is no longer required in make_pp_communicator.
             */
            const bool useCartesianReorder =
                    (ddSettings.useCartesianReorder && !cartSetup.bCartesianPP_PME);

            make_pp_communicator(mdlog, dd, cr, useCartesianReorder);
        }
    }
    else
    {
        receive_ddindex2simnodeid(dd, cr);
    }

    if (!thisRankHasDuty(cr, DUTY_PME))
    {
        /* Set up the commnuication to our PME node */
        dd->pme_nodeid = dd_simnode2pmenode(ddRankSetup, cartSetup, pmeRanks, cr, cr->sim_nodeid);
        dd->pme_receive_vir_ener = receive_vir_ener(dd, pmeRanks, cr);
        if (debug)
        {
            fprintf(debug,
                    "My pme_nodeid %d receive ener %s\n",
                    dd->pme_nodeid,
                    gmx::boolToString(dd->pme_receive_vir_ener));
        }
    }
    else
    {
        dd->pme_nodeid = -1;
    }

    /* We can not use DDMAIN(dd), because dd->mainrank is set later */
    if (MAIN(cr))
    {
        dd->ma = std::make_unique<AtomDistribution>(dd->numCells, numAtomsInSystem, numAtomsInSystem);
    }
}

static std::vector<real> get_slb_frac(const gmx::MDLogger& mdlog, const char* dir, int nc, const char* size_string)
{
    std::vector<real> slb_frac;
    if (nc > 1 && size_string != nullptr)
    {
        GMX_LOG(mdlog.info).appendTextFormatted("Using static load balancing for the %s direction", dir);
        slb_frac.resize(nc);
        real tot = 0;
        for (int i = 0; i < nc; i++)
        {
            double dbl = 0;
            int    n   = 0;
            sscanf(size_string, "%20lf%n", &dbl, &n);
            if (dbl == 0)
            {
                gmx_fatal(FARGS,
                          "Incorrect or not enough DD cell size entries for direction %s: '%s'",
                          dir,
                          size_string);
            }
            slb_frac[i] = dbl;
            size_string += n;
            tot += slb_frac[i];
        }
        /* Normalize */
        std::string relativeCellSizes = "Relative cell sizes:";
        for (int i = 0; i < nc; i++)
        {
            slb_frac[i] /= tot;
            relativeCellSizes += gmx::formatString(" %5.3f", slb_frac[i]);
        }
        GMX_LOG(mdlog.info).appendText(relativeCellSizes);
    }

    return slb_frac;
}

static int multi_body_bondeds_count(const gmx_mtop_t& mtop)
{
    int n = 0;
    for (const auto ilists : IListRange(mtop))
    {
        for (auto& ilist : extractILists(ilists.list(), IF_BOND))
        {
            if (NRAL(ilist.functionType) > 2)
            {
                n += ilists.nmol() * (ilist.iatoms.size() / ilistStride(ilist));
            }
        }
    }

    return n;
}

static int dd_getenv(const gmx::MDLogger& mdlog, const char* env_var, int def)
{
    int   nst = def;
    char* val = getenv(env_var);
    if (val)
    {
        if (sscanf(val, "%20d", &nst) <= 0)
        {
            nst = 1;
        }
        GMX_LOG(mdlog.info).appendTextFormatted("Found env.var. %s = %s, using value %d", env_var, val, nst);
    }

    return nst;
}

static void check_dd_restrictions(const gmx_domdec_t* dd, const t_inputrec& inputrec, const gmx::MDLogger& mdlog)
{
    if (inputrec.pbcType == PbcType::Screw
        && (dd->numCells[XX] == 1 || dd->numCells[YY] > 1 || dd->numCells[ZZ] > 1))
    {
        gmx_fatal(FARGS,
                  "With pbc=%s can only do domain decomposition in the x-direction",
                  c_pbcTypeNames[inputrec.pbcType].c_str());
    }

    if (inputrec.nstlist == 0)
    {
        gmx_fatal(FARGS, "Domain decomposition does not work with nstlist=0");
    }

    if (inputrec.comm_mode == ComRemovalAlgorithm::Angular && inputrec.pbcType != PbcType::No)
    {
        GMX_LOG(mdlog.warning)
                .appendText(
                        "comm-mode angular will give incorrect results when the comm group "
                        "partially crosses a periodic boundary");
    }
}

static real average_cellsize_min(const gmx_ddbox_t& ddbox, const ivec numDomains)
{
    real r = ddbox.box_size[XX];
    for (int d = 0; d < DIM; d++)
    {
        if (numDomains[d] > 1)
        {
            /* Check using the initial average cell size */
            r = std::min(r, ddbox.box_size[d] * ddbox.skew_fac[d] / numDomains[d]);
        }
    }

    return r;
}

/*! \brief Depending on the DLB initial value return the DLB switched off state or issue an error.
 */
static DlbState forceDlbOffOrBail(DlbState             cmdlineDlbState,
                                  const std::string&   reasonStr,
                                  const gmx::MDLogger& mdlog)
{
    std::string dlbNotSupportedErr = "Dynamic load balancing requested, but ";
    std::string dlbDisableNote     = "NOTE: disabling dynamic load balancing as ";

    if (cmdlineDlbState == DlbState::onUser)
    {
        gmx_fatal(FARGS, "%s", (dlbNotSupportedErr + reasonStr).c_str());
    }
    else if (cmdlineDlbState == DlbState::offCanTurnOn)
    {
        GMX_LOG(mdlog.info).appendText(dlbDisableNote + reasonStr);
    }
    return DlbState::offForever;
}

/*! \brief Return the dynamic load balancer's initial state based on initial conditions and user inputs.
 *
 * This function parses the parameters of "-dlb" command line option setting
 * corresponding state values. Then it checks the consistency of the determined
 * state with other run parameters and settings. As a result, the initial state
 * may be altered or an error may be thrown if incompatibility of options is detected.
 *
 * \param [in] mdlog                Logger.
 * \param [in] options              DomdecOptions.
 * \param [in] bRecordLoad          True if the load balancer is recording load information.
 * \param [in] mdrunOptions         Options for mdrun.
 * \param [in] inputrec             Pointer mdrun to input parameters.
 * \param [in] useGpuForPme         PME offloaded to GPU
 * \param [in] canUseGpuPmeDecomposition         GPU pme decomposition supported
 * \returns                         DLB initial/startup state.
 */
static DlbState determineInitialDlbState(const gmx::MDLogger&     mdlog,
                                         const DomdecOptions&     options,
                                         gmx_bool                 bRecordLoad,
                                         const gmx::MdrunOptions& mdrunOptions,
                                         const t_inputrec&        inputrec,
                                         const bool               useGpuForPme,
                                         const bool               canUseGpuPmeDecomposition)
{
    DlbState dlbState = DlbState::offCanTurnOn;

    switch (options.dlbOption)
    {
        case DlbOption::turnOnWhenUseful: dlbState = DlbState::offCanTurnOn; break;
        case DlbOption::no: dlbState = DlbState::offUser; break;
        case DlbOption::yes: dlbState = DlbState::onUser; break;
        default: gmx_incons("Invalid dlbOption enum value");
    }

    // Disable DLB when GPU PME decomposition is used
    // GPU PME decomposition uses an extended halo exchange algorithm without any X/F-redistribution
    // DLB causes changes in PP domain which doesn't work well with new PME decomposition algorithm

    // ToDo: This code has an assumption that options.numPmeRanks < 0 results in no separate PME rank (which is true currently),
    // if in future "-npme -1" option results in single separate PME rank, we shouldn't disable DLB in that case
    if (useGpuForPme && canUseGpuPmeDecomposition && (options.numPmeRanks == 0 || options.numPmeRanks > 1))
    {
        std::string reasonStr = "it is not supported with GPU PME decomposition.";
        return forceDlbOffOrBail(dlbState, reasonStr, mdlog);
    }

    /* Reruns don't support DLB: bail or override auto mode */
    if (mdrunOptions.rerun)
    {
        std::string reasonStr = "it is not supported in reruns.";
        return forceDlbOffOrBail(dlbState, reasonStr, mdlog);
    }

    /* Unsupported integrators */
    if (!EI_DYNAMICS(inputrec.eI))
    {
        auto reasonStr =
                gmx::formatString("it is only supported with dynamics, not with integrator '%s'.",
                                  enumValueToString(inputrec.eI));
        return forceDlbOffOrBail(dlbState, reasonStr, mdlog);
    }

    /* Without cycle counters we can't time work to balance on */
    if (!bRecordLoad)
    {
        std::string reasonStr =
                "cycle counters unsupported or not enabled in the operating system kernel.";
        return forceDlbOffOrBail(dlbState, reasonStr, mdlog);
    }

    if (mdrunOptions.reproducible)
    {
        std::string reasonStr = "you started a reproducible run.";
        switch (dlbState)
        {
            case DlbState::offUser: break;
            case DlbState::offForever:
                GMX_RELEASE_ASSERT(false, "DlbState::offForever is not a valid initial state");
                break;
            case DlbState::offCanTurnOn: return forceDlbOffOrBail(dlbState, reasonStr, mdlog);
            case DlbState::onCanTurnOff:
                GMX_RELEASE_ASSERT(false, "DlbState::offCanTurnOff is not a valid initial state");
                break;
            case DlbState::onUser:
                return forceDlbOffOrBail(
                        dlbState,
                        reasonStr
                                + " In load balanced runs binary reproducibility cannot be "
                                  "ensured.",
                        mdlog);
            default:
                gmx_fatal(FARGS,
                          "Death horror: undefined case (%d) for load balancing choice",
                          static_cast<int>(dlbState));
        }
    }

    return dlbState;
}

static std::unique_ptr<gmx_domdec_comm_t> init_dd_comm()
{
    auto comm = std::make_unique<gmx_domdec_comm_t>();

    comm->n_load_have    = 0;
    comm->n_load_collect = 0;

    comm->haveTurnedOffDlb = false;

    for (int i = 0; i < static_cast<int>(DDAtomRanges::Type::Number); i++)
    {
        comm->sum_nat[i] = 0;
    }
    comm->ndecomp   = 0;
    comm->nload     = 0;
    comm->load_step = 0;
    comm->load_sum  = 0;
    comm->load_max  = 0;
    clear_ivec(comm->load_lim);
    comm->load_mdf = 0;
    comm->load_pme = 0;

    return comm;
}

static void setupUpdateGroups(const gmx::MDLogger&              mdlog,
                              const gmx_mtop_t&                 mtop,
                              ArrayRef<const RangePartitioning> updateGroupingsPerMoleculeType,
                              const bool                        useUpdateGroups,
                              const real                        maxUpdateGroupRadius,
                              DDSystemInfo*                     systemInfo)
{
    systemInfo->updateGroupingsPerMoleculeType = updateGroupingsPerMoleculeType;
    systemInfo->useUpdateGroups                = useUpdateGroups;
    systemInfo->maxUpdateGroupRadius           = maxUpdateGroupRadius;

    if (systemInfo->useUpdateGroups)
    {
        int numUpdateGroups = 0;
        for (const auto& molblock : mtop.molblock)
        {
            numUpdateGroups += molblock.nmol
                               * systemInfo->updateGroupingsPerMoleculeType[molblock.type].numBlocks();
        }

        GMX_LOG(mdlog.info)
                .appendTextFormatted(
                        "Using update groups, nr %d, average size %.1f atoms, max. radius %.3f "
                        "nm\n",
                        numUpdateGroups,
                        mtop.natoms / static_cast<double>(numUpdateGroups),
                        systemInfo->maxUpdateGroupRadius);
    }
}

UnitCellInfo::UnitCellInfo(const t_inputrec& ir) :
    npbcdim(numPbcDimensions(ir.pbcType)),
    numBoundedDimensions(inputrec2nboundeddim(&ir)),
    ddBoxIsDynamic(numBoundedDimensions < DIM || inputrecDynamicBox(&ir)),
    haveScrewPBC(ir.pbcType == PbcType::Screw)
{
}

/* Returns whether molecules are always whole, i.e. not broken by PBC */
static bool moleculesAreAlwaysWhole(const gmx_mtop_t&                           mtop,
                                    const bool                                  useUpdateGroups,
                                    gmx::ArrayRef<const gmx::RangePartitioning> updateGroupingsPerMoleculeType)
{
    if (useUpdateGroups)
    {
        GMX_RELEASE_ASSERT(updateGroupingsPerMoleculeType.size() == mtop.moltype.size(),
                           "Need one grouping per moltype");
        for (size_t mol = 0; mol < mtop.moltype.size(); mol++)
        {
            if (updateGroupingsPerMoleculeType[mol].numBlocks() > 1)
            {
                return false;
            }
        }
    }
    else
    {
        for (const auto& moltype : mtop.moltype)
        {
            if (moltype.atoms.nr > 1)
            {
                return false;
            }
        }
    }

    return true;
}

/*! \brief Generate the simulation system information */
static DDSystemInfo getSystemInfo(const gmx::MDLogger&              mdlog,
                                  DDRole                            ddRole,
                                  MPI_Comm                          communicator,
                                  const DomdecOptions&              options,
                                  const gmx_mtop_t&                 mtop,
                                  const t_inputrec&                 ir,
                                  const matrix                      box,
                                  ArrayRef<const RangePartitioning> updateGroupingPerMoleculeType,
                                  const bool                        useUpdateGroups,
                                  const real                        maxUpdateGroupRadius,
                                  gmx::ArrayRef<const gmx::RVec>    xGlobal)
{
    const real tenPercentMargin = 1.1;

    DDSystemInfo systemInfo;

    setupUpdateGroups(
            mdlog, mtop, updateGroupingPerMoleculeType, useUpdateGroups, maxUpdateGroupRadius, &systemInfo);

    systemInfo.moleculesAreAlwaysWhole = moleculesAreAlwaysWhole(
            mtop, systemInfo.useUpdateGroups, systemInfo.updateGroupingsPerMoleculeType);
    systemInfo.haveInterDomainBondeds =
            (!systemInfo.moleculesAreAlwaysWhole || mtop.bIntermolecularInteractions);
    systemInfo.haveInterDomainMultiBodyBondeds =
            (systemInfo.haveInterDomainBondeds && multi_body_bondeds_count(mtop) > 0);

    if (systemInfo.useUpdateGroups)
    {
        systemInfo.mayHaveSplitConstraints = false;
        systemInfo.mayHaveSplitSettles     = false;
    }
    else
    {
        systemInfo.mayHaveSplitConstraints = (gmx_mtop_ftype_count(mtop, F_CONSTR) > 0
                                              || gmx_mtop_ftype_count(mtop, F_CONSTRNC) > 0);
        systemInfo.mayHaveSplitSettles     = (gmx_mtop_ftype_count(mtop, F_SETTLE) > 0);
    }

    if (ir.rlist == 0)
    {
        /* Set the cut-off to some very large value,
         * so we don't need if statements everywhere in the code.
         * We use sqrt, since the cut-off is squared in some places.
         */
        systemInfo.cutoff = GMX_CUTOFF_INF;
    }
    else
    {
        systemInfo.cutoff = atomToAtomIntoDomainToDomainCutoff(systemInfo, ir.rlist);
    }
    systemInfo.minCutoffForMultiBody = 0;

    /* Determine the minimum cell size limit, affected by many factors */
    systemInfo.cellsizeLimit             = 0;
    systemInfo.filterBondedCommunication = false;

    /* We do not allow home atoms to move beyond the neighboring domain
     * between domain decomposition steps, which limits the cell size.
     * Get an estimate of cell size limit due to atom displacement.
     * In most cases this is a large overestimate, because it assumes
     * non-interaction atoms.
     * We set the chance to 1 in a trillion steps.
     * Note that any atom in the system should not have a too large
     * displacement. Thus we use ChanceTarget::System. This means that
     * the minimum cell size increases (slowly) with the sytem size.
     */
    constexpr real c_chanceThatAtomMovesBeyondDomain = 1e-12;
    const real     limitForAtomDisplacement =
            minCellSizeForAtomDisplacement(mtop,
                                           ir,
                                           systemInfo.updateGroupingsPerMoleculeType,
                                           c_chanceThatAtomMovesBeyondDomain,
                                           ChanceTarget::System);
    GMX_LOG(mdlog.info).appendTextFormatted("Minimum cell size due to atom displacement: %.3f nm", limitForAtomDisplacement);

    systemInfo.cellsizeLimit = std::max(systemInfo.cellsizeLimit, limitForAtomDisplacement);

    /* TODO: PME decomposition currently requires atoms not to be more than
     *       2/3 of comm->cutoff, which is >=rlist, outside of their domain.
     *       In nearly all cases, limitForAtomDisplacement will be smaller
     *       than 2/3*rlist, so the PME requirement is satisfied.
     *       But it would be better for both correctness and performance
     *       to use limitForAtomDisplacement instead of 2/3*comm->cutoff.
     *       Note that we would need to improve the pairlist buffer case.
     */

    if (systemInfo.haveInterDomainBondeds)
    {
        if (options.minimumCommunicationRange > 0)
        {
            systemInfo.minCutoffForMultiBody =
                    atomToAtomIntoDomainToDomainCutoff(systemInfo, options.minimumCommunicationRange);
            if (options.useBondedCommunication)
            {
                systemInfo.filterBondedCommunication =
                        (systemInfo.minCutoffForMultiBody > systemInfo.cutoff);
            }
            else
            {
                systemInfo.cutoff = std::max(systemInfo.cutoff, systemInfo.minCutoffForMultiBody);
            }
        }
        else
        {
            real r_2b = 0;
            real r_mb = 0;

            if (ddRole == DDRole::Main)
            {
                dd_bonded_cg_distance(mdlog, mtop, ir, xGlobal, box, options.ddBondedChecking, &r_2b, &r_mb);
            }
            gmx_bcast(sizeof(r_2b), &r_2b, communicator);
            gmx_bcast(sizeof(r_mb), &r_mb, communicator);

            /* We use an initial margin of 10% for the minimum cell size,
             * except when we are just below the non-bonded cut-off.
             */
            if (options.useBondedCommunication)
            {
                if (std::max(r_2b, r_mb) > systemInfo.cutoff)
                {
                    const real r_bonded              = std::max(r_2b, r_mb);
                    systemInfo.minCutoffForMultiBody = tenPercentMargin * r_bonded;
                    /* This is the (only) place where we turn on the filtering */
                    systemInfo.filterBondedCommunication = true;
                }
                else
                {
                    const real r_bonded = r_mb;
                    systemInfo.minCutoffForMultiBody =
                            std::min(tenPercentMargin * r_bonded, systemInfo.cutoff);
                }
                /* We determine cutoff_mbody later */
                systemInfo.increaseMultiBodyCutoff = true;
            }
            else
            {
                /* No special bonded communication,
                 * simply increase the DD cut-off.
                 */
                systemInfo.minCutoffForMultiBody = tenPercentMargin * std::max(r_2b, r_mb);
                systemInfo.cutoff = std::max(systemInfo.cutoff, systemInfo.minCutoffForMultiBody);
            }
        }
        GMX_LOG(mdlog.info)
                .appendTextFormatted("Minimum cell size due to bonded interactions: %.3f nm",
                                     systemInfo.minCutoffForMultiBody);

        systemInfo.cellsizeLimit = std::max(systemInfo.cellsizeLimit, systemInfo.minCutoffForMultiBody);
    }

    systemInfo.constraintCommunicationRange = 0;
    if (systemInfo.mayHaveSplitConstraints && options.constraintCommunicationRange <= 0)
    {
        /* There is a cell size limit due to the constraints (P-LINCS) */
        systemInfo.constraintCommunicationRange = gmx::constr_r_max(mdlog, &mtop, &ir);
        GMX_LOG(mdlog.info)
                .appendTextFormatted("Estimated maximum distance required for P-LINCS: %.3f nm",
                                     systemInfo.constraintCommunicationRange);
        if (systemInfo.constraintCommunicationRange > systemInfo.cellsizeLimit)
        {
            GMX_LOG(mdlog.info)
                    .appendText(
                            "This distance will limit the DD cell size, you can override this with "
                            "-rcon");
        }
    }
    else if (options.constraintCommunicationRange > 0)
    {
        /* Here we do not check for dd->splitConstraints.
         * because one can also set a cell size limit for virtual sites only
         * and at this point we don't know yet if there are intercg v-sites.
         */
        GMX_LOG(mdlog.info)
                .appendTextFormatted("User supplied maximum distance required for P-LINCS: %.3f nm",
                                     options.constraintCommunicationRange);
        systemInfo.constraintCommunicationRange = options.constraintCommunicationRange;
    }
    systemInfo.cellsizeLimit = std::max(systemInfo.cellsizeLimit, systemInfo.constraintCommunicationRange);

    systemInfo.haveBoxDeformation = ir_haveBoxDeformation(ir);
    setBoxDeformationRate(ir.deform, box, systemInfo.boxDeformationRate);

    return systemInfo;
}

/*! \brief Exit with a fatal error if the DDGridSetup cannot be
 * implemented. */
static void checkDDGridSetup(const DDGridSetup&   ddGridSetup,
                             DDRole               ddRole,
                             MPI_Comm             communicator,
                             int                  numNodes,
                             const DomdecOptions& options,
                             const DDSettings&    ddSettings,
                             const DDSystemInfo&  systemInfo,
                             const real           cellsizeLimit,
                             const gmx_ddbox_t&   ddbox)
{
    if (options.numCells[XX] <= 0 && (ddGridSetup.numDomains[XX] == 0))
    {
        const bool  bC = (systemInfo.mayHaveSplitConstraints
                         && systemInfo.constraintCommunicationRange > systemInfo.minCutoffForMultiBody);
        std::string message =
                gmx::formatString("Change the number of ranks or mdrun option %s%s%s",
                                  !bC ? "-rdd" : "-rcon",
                                  ddSettings.initialDlbState != DlbState::offUser ? " or -dds" : "",
                                  bC ? " or your LINCS settings" : "");

        gmx_fatal_collective(FARGS,
                             communicator,
                             ddRole == DDRole::Main,
                             "There is no domain decomposition for %d ranks that is compatible "
                             "with the given box and a minimum cell size of %g nm\n"
                             "%s\n"
                             "Look in the log file for details on the domain decomposition",
                             numNodes - ddGridSetup.numPmeOnlyRanks,
                             cellsizeLimit,
                             message.c_str());
    }

    const real acs = average_cellsize_min(ddbox, ddGridSetup.numDomains);
    if (acs < cellsizeLimit)
    {
        if (options.numCells[XX] <= 0)
        {
            GMX_RELEASE_ASSERT(
                    false,
                    "dd_choose_grid() should return a grid that satisfies the cell size limits");
        }
        else
        {
            gmx_fatal_collective(
                    FARGS,
                    communicator,
                    ddRole == DDRole::Main,
                    "The initial cell size (%f) is smaller than the cell size limit (%f), change "
                    "options -dd, -rdd or -rcon, see the log file for details",
                    acs,
                    cellsizeLimit);
        }
    }

    const int numPPRanks =
            ddGridSetup.numDomains[XX] * ddGridSetup.numDomains[YY] * ddGridSetup.numDomains[ZZ];
    if (numNodes - numPPRanks != ddGridSetup.numPmeOnlyRanks)
    {
        gmx_fatal_collective(FARGS,
                             communicator,
                             ddRole == DDRole::Main,
                             "The size of the domain decomposition grid (%d) does not match the "
                             "number of PP ranks (%d). The total number of ranks is %d",
                             numPPRanks,
                             numNodes - ddGridSetup.numPmeOnlyRanks,
                             numNodes);
    }
    if (ddGridSetup.numPmeOnlyRanks > numPPRanks)
    {
        gmx_fatal_collective(FARGS,
                             communicator,
                             ddRole == DDRole::Main,
                             "The number of separate PME ranks (%d) is larger than the number of "
                             "PP ranks (%d), this is not supported.",
                             ddGridSetup.numPmeOnlyRanks,
                             numPPRanks);
    }
}

/*! \brief Set the cell size and interaction limits, as well as the DD grid */
static DDRankSetup getDDRankSetup(const gmx::MDLogger& mdlog,
                                  int                  numNodes,
                                  const DdRankOrder    rankOrder,
                                  const DDGridSetup&   ddGridSetup,
                                  const t_inputrec&    ir)
{
    GMX_LOG(mdlog.info)
            .appendTextFormatted("Domain decomposition grid %d x %d x %d, separate PME ranks %d",
                                 ddGridSetup.numDomains[XX],
                                 ddGridSetup.numDomains[YY],
                                 ddGridSetup.numDomains[ZZ],
                                 ddGridSetup.numPmeOnlyRanks);

    DDRankSetup ddRankSetup;

    ddRankSetup.rankOrder = rankOrder;

    ddRankSetup.numPPRanks = numNodes - ddGridSetup.numPmeOnlyRanks;
    copy_ivec(ddGridSetup.numDomains, ddRankSetup.numPPCells);

    ddRankSetup.usePmeOnlyRanks = (ddGridSetup.numPmeOnlyRanks > 0);
    if (ddRankSetup.usePmeOnlyRanks)
    {
        ddRankSetup.numRanksDoingPme = ddGridSetup.numPmeOnlyRanks;
    }
    else
    {
        ddRankSetup.numRanksDoingPme =
                ddGridSetup.numDomains[XX] * ddGridSetup.numDomains[YY] * ddGridSetup.numDomains[ZZ];
    }

    if (usingPme(ir.coulombtype) || usingLJPme(ir.vdwtype))
    {
        /* The following choices should match those
         * in comm_cost_est in domdec_setup.c.
         * Note that here the checks have to take into account
         * that the decomposition might occur in a different order than xyz
         * (for instance through the env.var. GMX_DD_ORDER_ZYX),
         * in which case they will not match those in comm_cost_est,
         * but since that is mainly for testing purposes that's fine.
         */
        if (ddGridSetup.numDDDimensions >= 2 && ddGridSetup.ddDimensions[0] == XX
            && ddGridSetup.ddDimensions[1] == YY
            && ddRankSetup.numRanksDoingPme > ddGridSetup.numDomains[XX]
            && ddRankSetup.numRanksDoingPme % ddGridSetup.numDomains[XX] == 0
            && getenv("GMX_PMEONEDD") == nullptr)
        {
            ddRankSetup.npmedecompdim = 2;
            ddRankSetup.npmenodes_x   = ddGridSetup.numDomains[XX];
            ddRankSetup.npmenodes_y   = ddRankSetup.numRanksDoingPme / ddRankSetup.npmenodes_x;
        }
        else
        {
            /* In case nc is 1 in both x and y we could still choose to
             * decompose pme in y instead of x, but we use x for simplicity.
             */
            ddRankSetup.npmedecompdim = 1;
            if (ddGridSetup.ddDimensions[0] == YY)
            {
                ddRankSetup.npmenodes_x = 1;
                ddRankSetup.npmenodes_y = ddRankSetup.numRanksDoingPme;
            }
            else
            {
                ddRankSetup.npmenodes_x = ddRankSetup.numRanksDoingPme;
                ddRankSetup.npmenodes_y = 1;
            }
        }
        GMX_LOG(mdlog.info)
                .appendTextFormatted("PME domain decomposition: %d x %d x %d",
                                     ddRankSetup.npmenodes_x,
                                     ddRankSetup.npmenodes_y,
                                     1);
    }
    else
    {
        ddRankSetup.npmedecompdim = 0;
        ddRankSetup.npmenodes_x   = 0;
        ddRankSetup.npmenodes_y   = 0;
    }

    return ddRankSetup;
}

/*! \brief Set the cell size and interaction limits */
static void set_dd_limits(const gmx::MDLogger& mdlog,
                          DDRole               ddRole,
                          gmx_domdec_t*        dd,
                          const DomdecOptions& options,
                          const DDSettings&    ddSettings,
                          const DDSystemInfo&  systemInfo,
                          const DDGridSetup&   ddGridSetup,
                          const int            numPPRanks,
                          const gmx_mtop_t&    mtop,
                          const t_inputrec&    ir,
                          const gmx_ddbox_t&   ddbox)
{
    gmx_domdec_comm_t* comm = dd->comm.get();
    comm->ddSettings        = ddSettings;

    /* Initialize to GPU share count to 0, might change later */
    comm->nrank_gpu_shared = 0;

    comm->dlbState = comm->ddSettings.initialDlbState;
    dd_dlb_set_should_check_whether_to_turn_dlb_on(dd, TRUE);
    /* To consider turning DLB on after 2*nstlist steps we need to check
     * at partitioning count 3. Thus we need to increase the first count by 2.
     */
    comm->ddPartioningCountFirstDlbOff += 2;

    comm->bPMELoadBalDLBLimits = FALSE;

    /* Allocate the charge group/atom sorting struct */
    comm->sort = std::make_unique<gmx_domdec_sort_t>();

    comm->systemInfo = systemInfo;

    if (systemInfo.useUpdateGroups)
    {
        /* Note: We would like to use dd->nnodes for the atom count estimate,
         *       but that is not yet available here. But this anyhow only
         *       affect performance up to the second dd_partition_system call.
         */
        const int homeAtomCountEstimate = mtop.natoms / numPPRanks;
        comm->updateGroupsCog           = std::make_unique<gmx::UpdateGroupsCog>(
                mtop, systemInfo.updateGroupingsPerMoleculeType, maxReferenceTemperature(ir), homeAtomCountEstimate);
    }

    /* Set the DD setup given by ddGridSetup */
    copy_ivec(ddGridSetup.numDomains, dd->numCells);
    dd->ndim = ddGridSetup.numDDDimensions;
    copy_ivec(ddGridSetup.ddDimensions, dd->dim);

    dd->nnodes = dd->numCells[XX] * dd->numCells[YY] * dd->numCells[ZZ];

    if (isDlbDisabled(comm->dlbState))
    {
        comm->slb_frac[XX] = get_slb_frac(mdlog, "x", dd->numCells[XX], options.cellSizeX);
        comm->slb_frac[YY] = get_slb_frac(mdlog, "y", dd->numCells[YY], options.cellSizeY);
        comm->slb_frac[ZZ] = get_slb_frac(mdlog, "z", dd->numCells[ZZ], options.cellSizeZ);
    }

    /* Set the multi-body cut-off and cellsize limit for DLB */
    comm->cutoff_mbody   = systemInfo.minCutoffForMultiBody;
    comm->cellsize_limit = systemInfo.cellsizeLimit;
    if (systemInfo.haveInterDomainBondeds && systemInfo.increaseMultiBodyCutoff)
    {
        if (systemInfo.filterBondedCommunication || !isDlbDisabled(comm->dlbState))
        {
            /* Set the bonded communication distance to halfway
             * the minimum and the maximum,
             * since the extra communication cost is nearly zero.
             */
            real acs           = average_cellsize_min(ddbox, dd->numCells);
            comm->cutoff_mbody = 0.5 * (systemInfo.minCutoffForMultiBody + acs);
            if (!isDlbDisabled(comm->dlbState))
            {
                /* Check if this does not limit the scaling */
                comm->cutoff_mbody = std::min(comm->cutoff_mbody, options.dlbScaling * acs);
            }
            if (!systemInfo.filterBondedCommunication)
            {
                /* Without bBondComm do not go beyond the n.b. cut-off */
                comm->cutoff_mbody = std::min(comm->cutoff_mbody, systemInfo.cutoff);
                if (comm->cellsize_limit >= systemInfo.cutoff)
                {
                    /* We don't loose a lot of efficieny
                     * when increasing it to the n.b. cut-off.
                     * It can even be slightly faster, because we need
                     * less checks for the communication setup.
                     */
                    comm->cutoff_mbody = systemInfo.cutoff;
                }
            }
            /* Check if we did not end up below our original limit */
            comm->cutoff_mbody = std::max(comm->cutoff_mbody, systemInfo.minCutoffForMultiBody);

            if (comm->cutoff_mbody > comm->cellsize_limit)
            {
                comm->cellsize_limit = comm->cutoff_mbody;
            }
        }
        /* Without DLB and cutoff_mbody<cutoff, cutoff_mbody is dynamic */
    }

    if (debug)
    {
        fprintf(debug,
                "Bonded atom communication beyond the cut-off: %s\n"
                "cellsize limit %f\n",
                gmx::boolToString(systemInfo.filterBondedCommunication),
                comm->cellsize_limit);
    }

    if (ddRole == DDRole::Main)
    {
        check_dd_restrictions(dd, ir, mdlog);
    }
}

static void writeSettings(gmx::TextWriter*   log,
                          gmx_domdec_t*      dd,
                          const gmx_mtop_t&  mtop,
                          const t_inputrec&  ir,
                          gmx_bool           bDynLoadBal,
                          real               dlb_scale,
                          const gmx_ddbox_t* ddbox)
{
    gmx_domdec_comm_t* comm = dd->comm.get();

    if (bDynLoadBal)
    {
        log->writeString("The maximum number of communication pulses is:");
        for (int d = 0; d < dd->ndim; d++)
        {
            log->writeStringFormatted(" %c %d", dim2char(dd->dim[d]), comm->cd[d].np_dlb);
        }
        log->ensureLineBreak();
        log->writeLineFormatted("The minimum size for domain decomposition cells is %.3f nm",
                                comm->cellsize_limit);
        log->writeLineFormatted("The requested allowed shrink of DD cells (option -dds) is: %.2f", dlb_scale);
        log->writeString("The allowed shrink of domain decomposition cells is:");
        for (int d = 0; d < DIM; d++)
        {
            if (dd->numCells[d] > 1)
            {
                const real shrink =
                        (d >= ddbox->npbcdim && dd->numCells[d] == 2)
                                ? 0
                                : comm->cellsize_min_dlb[d]
                                          / (ddbox->box_size[d] * ddbox->skew_fac[d] / dd->numCells[d]);
                log->writeStringFormatted(" %c %.2f", dim2char(d), shrink);
            }
        }
        log->ensureLineBreak();
    }
    else
    {
        ivec np;
        set_dd_cell_sizes_slb(dd, ddbox, setcellsizeslbPULSE_ONLY, np);
        log->writeString("The initial number of communication pulses is:");
        for (int d = 0; d < dd->ndim; d++)
        {
            log->writeStringFormatted(" %c %d", dim2char(dd->dim[d]), np[dd->dim[d]]);
        }
        log->ensureLineBreak();
        log->writeString("The initial domain decomposition cell size is:");
        for (int d = 0; d < DIM; d++)
        {
            if (dd->numCells[d] > 1)
            {
                log->writeStringFormatted(" %c %.2f nm", dim2char(d), dd->comm->cellsize_min[d]);
            }
        }
        log->ensureLineBreak();
        log->writeLine();
    }

    const bool haveInterDomainVsites =
            (countInterUpdategroupVsites(mtop, comm->systemInfo.updateGroupingsPerMoleculeType) != 0);

    if (comm->systemInfo.haveInterDomainBondeds || haveInterDomainVsites
        || comm->systemInfo.mayHaveSplitConstraints || comm->systemInfo.mayHaveSplitSettles)
    {
        std::string decompUnits;
        if (comm->systemInfo.useUpdateGroups)
        {
            decompUnits = "atom groups";
        }
        else
        {
            decompUnits = "atoms";
        }

        log->writeLineFormatted("The maximum allowed distance for %s involved in interactions is:",
                                decompUnits.c_str());
        log->writeLineFormatted(
                "%40s  %-7s %6.3f nm", "non-bonded interactions", "", comm->systemInfo.cutoff);

        real limit = 0;
        if (bDynLoadBal)
        {
            limit = dd->comm->cellsize_limit;
        }
        else
        {
            if (dd->unitCellInfo.ddBoxIsDynamic)
            {
                log->writeLine(
                        "(the following are initial values, they could change due to box "
                        "deformation)");
            }
            limit = dd->comm->cellsize_min[XX];
            for (int d = 1; d < DIM; d++)
            {
                limit = std::min(limit, dd->comm->cellsize_min[d]);
            }
        }

        if (comm->systemInfo.haveInterDomainBondeds)
        {
            log->writeLineFormatted("%40s  %-7s %6.3f nm",
                                    "two-body bonded interactions",
                                    "(-rdd)",
                                    std::max(comm->systemInfo.cutoff, comm->cutoff_mbody));
            log->writeLineFormatted(
                    "%40s  %-7s %6.3f nm",
                    "multi-body bonded interactions",
                    "(-rdd)",
                    (comm->systemInfo.filterBondedCommunication || isDlbOn(dd->comm->dlbState))
                            ? comm->cutoff_mbody
                            : std::min(comm->systemInfo.cutoff, limit));
        }
        if (haveInterDomainVsites)
        {
            log->writeLineFormatted("%40s  %-7s %6.3f nm", "virtual site constructions", "(-rcon)", limit);
        }
        if (comm->systemInfo.mayHaveSplitConstraints || comm->systemInfo.mayHaveSplitSettles)
        {
            std::string separation =
                    gmx::formatString("atoms separated by up to %d constraints", 1 + ir.nProjOrder);
            log->writeLineFormatted("%40s  %-7s %6.3f nm\n", separation.c_str(), "(-rcon)", limit);
        }
        log->ensureLineBreak();
    }
}

static void logSettings(const gmx::MDLogger& mdlog,
                        gmx_domdec_t*        dd,
                        const gmx_mtop_t&    mtop,
                        const t_inputrec&    ir,
                        real                 dlb_scale,
                        const gmx_ddbox_t*   ddbox)
{
    gmx::StringOutputStream stream;
    gmx::TextWriter         log(&stream);
    writeSettings(&log, dd, mtop, ir, isDlbOn(dd->comm->dlbState), dlb_scale, ddbox);
    if (dd->comm->dlbState == DlbState::offCanTurnOn)
    {
        {
            log.ensureEmptyLine();
            log.writeLine(
                    "When dynamic load balancing gets turned on, these settings will change to:");
        }
        writeSettings(&log, dd, mtop, ir, true, dlb_scale, ddbox);
    }
    GMX_LOG(mdlog.info).asParagraph().appendText(stream.toString());
}

static void set_cell_limits_dlb(const gmx::MDLogger& mdlog,
                                gmx_domdec_t*        dd,
                                real                 dlb_scale,
                                const t_inputrec&    inputrec,
                                const gmx_ddbox_t*   ddbox)
{
    int npulse       = 0;
    int npulse_d_max = 0;
    int npulse_d     = 0;

    gmx_domdec_comm_t* comm = dd->comm.get();

    bool bNoCutOff = (inputrec.rvdw == 0 || inputrec.rcoulomb == 0);

    /* Determine the maximum number of comm. pulses in one dimension */

    comm->cellsize_limit = std::max(comm->cellsize_limit, comm->cutoff_mbody);

    /* Determine the maximum required number of grid pulses */
    if (comm->cellsize_limit >= comm->systemInfo.cutoff)
    {
        /* Only a single pulse is required */
        npulse = 1;
    }
    else if (!bNoCutOff && comm->cellsize_limit > 0)
    {
        /* We round down slightly here to avoid overhead due to the latency
         * of extra communication calls when the cut-off
         * would be only slightly longer than the cell size.
         * Later cellsize_limit is redetermined,
         * so we can not miss interactions due to this rounding.
         */
        npulse = static_cast<int>(0.96 + comm->systemInfo.cutoff / comm->cellsize_limit);
    }
    else
    {
        /* There is no cell size limit */
        npulse = std::max(dd->numCells[XX] - 1, std::max(dd->numCells[YY] - 1, dd->numCells[ZZ] - 1));
    }

    if (!bNoCutOff && npulse > 1)
    {
        /* See if we can do with less pulses, based on dlb_scale */
        npulse_d_max = 0;
        for (int d = 0; d < dd->ndim; d++)
        {
            int dim  = dd->dim[d];
            npulse_d = static_cast<int>(
                    1
                    + dd->numCells[dim] * comm->systemInfo.cutoff
                              / (ddbox->box_size[dim] * ddbox->skew_fac[dim] * dlb_scale));
            npulse_d_max = std::max(npulse_d_max, npulse_d);
        }
        npulse = std::min(npulse, npulse_d_max);
    }

    /* This env var can override npulse */
    const int ddPulseEnv = dd_getenv(mdlog, "GMX_DD_NPULSE", 0);
    if (ddPulseEnv > 0)
    {
        npulse = ddPulseEnv;
    }

    comm->maxpulse       = 1;
    comm->bVacDLBNoLimit = (inputrec.pbcType == PbcType::No);
    for (int d = 0; d < dd->ndim; d++)
    {
        comm->cd[d].np_dlb = std::min(npulse, dd->numCells[dd->dim[d]] - 1);
        comm->maxpulse     = std::max(comm->maxpulse, comm->cd[d].np_dlb);
        if (comm->cd[d].np_dlb < dd->numCells[dd->dim[d]] - 1)
        {
            comm->bVacDLBNoLimit = FALSE;
        }
    }

    /* cellsize_limit is set for LINCS in init_domain_decomposition */
    if (!comm->bVacDLBNoLimit)
    {
        comm->cellsize_limit = std::max(comm->cellsize_limit, comm->systemInfo.cutoff / comm->maxpulse);
    }
    comm->cellsize_limit = std::max(comm->cellsize_limit, comm->cutoff_mbody);
    /* Set the minimum cell size for each DD dimension */
    for (int d = 0; d < dd->ndim; d++)
    {
        if (comm->bVacDLBNoLimit || comm->cd[d].np_dlb * comm->cellsize_limit >= comm->systemInfo.cutoff)
        {
            comm->cellsize_min_dlb[dd->dim[d]] = comm->cellsize_limit;
        }
        else
        {
            comm->cellsize_min_dlb[dd->dim[d]] = comm->systemInfo.cutoff / comm->cd[d].np_dlb;
        }
    }
    if (comm->cutoff_mbody <= 0)
    {
        comm->cutoff_mbody = std::min(comm->systemInfo.cutoff, comm->cellsize_limit);
    }
    if (isDlbOn(comm->dlbState))
    {
        set_dlb_limits(dd);
    }
}

bool dd_moleculesAreAlwaysWhole(const gmx_domdec_t& dd)
{
    return dd.comm->systemInfo.moleculesAreAlwaysWhole;
}

bool dd_bonded_molpbc(const gmx_domdec_t& dd, PbcType pbcType)
{
    /* If each molecule is a single charge group
     * or we use domain decomposition for each periodic dimension,
     * we do not need to take pbc into account for the bonded interactions.
     */
    return (pbcType != PbcType::No && dd.comm->systemInfo.haveInterDomainBondeds
            && !(dd.numCells[XX] > 1 && dd.numCells[YY] > 1
                 && (dd.numCells[ZZ] > 1 || pbcType == PbcType::XY)));
}

/*! \brief Sets grid size limits and PP-PME setup, prints settings to log */
static void set_ddgrid_parameters(const gmx::MDLogger& mdlog,
                                  gmx_domdec_t*        dd,
                                  real                 dlb_scale,
                                  const gmx_mtop_t&    mtop,
                                  const t_inputrec&    inputrec,
                                  const gmx_ddbox_t*   ddbox)
{
    gmx_domdec_comm_t* comm        = dd->comm.get();
    DDRankSetup&       ddRankSetup = comm->ddRankSetup;

    if (usingPme(inputrec.coulombtype) || usingLJPme(inputrec.vdwtype))
    {
        init_ddpme(dd, &ddRankSetup.ddpme[0], 0);
        if (ddRankSetup.npmedecompdim >= 2)
        {
            init_ddpme(dd, &ddRankSetup.ddpme[1], 1);
        }
    }
    else
    {
        ddRankSetup.numRanksDoingPme = 0;
        if (dd->pme_nodeid >= 0)
        {
            gmx_fatal_collective(FARGS,
                                 dd->mpi_comm_all,
                                 DDMAIN(dd),
                                 "Can not have separate PME ranks without PME electrostatics");
        }
    }

    if (debug)
    {
        fprintf(debug, "The DD cut-off is %f\n", comm->systemInfo.cutoff);
    }
    if (!isDlbDisabled(comm->dlbState))
    {
        set_cell_limits_dlb(mdlog, dd, dlb_scale, inputrec, ddbox);
    }

    // Make nstDDGlobalComm the first multiple of nstlist >= c_minimumGCStepInterval
    dd->comm->nstDDGlobalComm = sc_minimumGCStepInterval;
    if (dd->comm->nstDDGlobalComm % inputrec.nstlist != 0)
    {
        dd->comm->nstDDGlobalComm += inputrec.nstlist - (dd->comm->nstDDGlobalComm % inputrec.nstlist);
    }

    logSettings(mdlog, dd, mtop, inputrec, dlb_scale, ddbox);

    const real vol_frac = (inputrec.pbcType == PbcType::No)
                                  ? (1 - 1 / static_cast<double>(dd->nnodes))
                                  : ((1 + comm_box_frac(dd->numCells, comm->systemInfo.cutoff, *ddbox))
                                     / static_cast<double>(dd->nnodes));
    if (debug)
    {
        fprintf(debug, "Volume fraction for all DD zones: %f\n", vol_frac);
    }
    int natoms_tot = mtop.natoms;

    dd->ga2la = std::make_unique<gmx_ga2la_t>(natoms_tot, static_cast<int>(vol_frac * natoms_tot));
}

/*! \brief Get some important DD parameters which can be modified by env.vars */
static DDSettings getDDSettings(const gmx::MDLogger&     mdlog,
                                const DomdecOptions&     options,
                                const gmx::MdrunOptions& mdrunOptions,
                                const t_inputrec&        ir,
                                const bool               useGpuForPme,
                                const bool               canUseGpuPmeDecomposition)
{
    DDSettings ddSettings;

    ddSettings.useSendRecv2        = (dd_getenv(mdlog, "GMX_DD_USE_SENDRECV2", 0) != 0);
    ddSettings.dlb_scale_lim       = dd_getenv(mdlog, "GMX_DLB_MAX_BOX_SCALING", 10);
    ddSettings.useDDOrderZYX       = bool(dd_getenv(mdlog, "GMX_DD_ORDER_ZYX", 0));
    ddSettings.useCartesianReorder = bool(dd_getenv(mdlog, "GMX_NO_CART_REORDER", 1));
    ddSettings.eFlop               = dd_getenv(mdlog, "GMX_DLB_BASED_ON_FLOPS", 0);
    const int recload              = dd_getenv(mdlog, "GMX_DD_RECORD_LOAD", 1);
    ddSettings.nstDDDump           = dd_getenv(mdlog, "GMX_DD_NST_DUMP", 0);
    ddSettings.nstDDDumpGrid       = dd_getenv(mdlog, "GMX_DD_NST_DUMP_GRID", 0);
    ddSettings.DD_debug            = dd_getenv(mdlog, "GMX_DD_DEBUG", 0);

    if (ddSettings.useSendRecv2)
    {
        GMX_LOG(mdlog.info)
                .appendText(
                        "Will use two sequential MPI_Sendrecv calls instead of two simultaneous "
                        "non-blocking MPI_Irecv and MPI_Isend pairs for constraint and vsite "
                        "communication");
    }

    if (ddSettings.eFlop)
    {
        GMX_LOG(mdlog.info).appendText("Will load balance based on FLOP count");
        ddSettings.recordLoad = true;
    }
    else
    {
        ddSettings.recordLoad = (wallcycle_have_counter() && recload > 0);
    }

    ddSettings.initialDlbState = determineInitialDlbState(
            mdlog, options, ddSettings.recordLoad, mdrunOptions, ir, useGpuForPme, canUseGpuPmeDecomposition);
    GMX_LOG(mdlog.info)
            .appendTextFormatted("Dynamic load balancing: %s",
                                 enumValueToString(ddSettings.initialDlbState));

    return ddSettings;
}

gmx_domdec_t::gmx_domdec_t(const t_inputrec& ir, gmx::ArrayRef<const int> ddDims) :
    unitCellInfo(ir), zones(ddDims)
{
}

gmx_domdec_t::~gmx_domdec_t() = default;

namespace gmx
{

// TODO once the functionality stablizes, move this class and
// supporting functionality into builder.cpp
/*! \brief Impl class for DD builder */
class DomainDecompositionBuilder::Impl
{
public:
    //! Constructor
    Impl(const MDLogger&                   mdlog,
         t_commrec*                        cr,
         const DomdecOptions&              options,
         const MdrunOptions&               mdrunOptions,
         const gmx_mtop_t&                 mtop,
         const t_inputrec&                 ir,
         const MDModulesNotifiers&         notifiers,
         const matrix                      box,
         ArrayRef<const RangePartitioning> updateGroupingPerMoleculeType,
         bool                              useUpdateGroups,
         real                              maxUpdateGroupRadius,
         ArrayRef<const RVec>              xGlobal,
         bool                              useGpuForNonbonded,
         bool                              useGpuForPme,
         bool                              useGpuForUpdate,
         bool                              useGpuDirectHalo,
         bool                              canUseGpuPmeDecomposition);

    //! Build the resulting DD manager
    std::unique_ptr<gmx_domdec_t> build(LocalAtomSetManager*       atomSets,
                                        const gmx_localtop_t&      localTopology,
                                        const t_state&             localState,
                                        ObservablesReducerBuilder* observablesReducerBuilder);

    //! Objects used in constructing and configuring DD
    //! {
    //! Logging object
    const MDLogger& mdlog_;
    //! Communication object
    t_commrec* cr_;
    //! User-supplied options configuring DD behavior
    const DomdecOptions options_;
    //! Global system topology
    const gmx_mtop_t& mtop_;
    //! User input values from the tpr file
    const t_inputrec& ir_;
    //! MdModules object
    const MDModulesNotifiers& notifiers_;
    //! }

    //! Internal objects used in constructing DD
    //! {
    //! Settings combined from the user input
    DDSettings ddSettings_;
    //! Information derived from the simulation system
    DDSystemInfo systemInfo_;
    //! Box structure
    gmx_ddbox_t ddbox_ = { 0 };
    //! Organization of the DD grids
    DDGridSetup ddGridSetup_;
    //! Organzation of the DD ranks
    DDRankSetup ddRankSetup_;
    //! Number of DD cells in each dimension
    ivec ddCellIndex_ = { 0, 0, 0 };
    //! IDs of PME-only ranks
    std::vector<int> pmeRanks_;
    //! Contains a valid Cartesian-communicator-based setup, or defaults.
    CartesianRankSetup cartSetup_;
    //! }
};

DomainDecompositionBuilder::Impl::Impl(const MDLogger&                   mdlog,
                                       t_commrec*                        cr,
                                       const DomdecOptions&              options,
                                       const MdrunOptions&               mdrunOptions,
                                       const gmx_mtop_t&                 mtop,
                                       const t_inputrec&                 ir,
                                       const MDModulesNotifiers&         notifiers,
                                       const matrix                      box,
                                       ArrayRef<const RangePartitioning> updateGroupingPerMoleculeType,
                                       const bool                        useUpdateGroups,
                                       const real                        maxUpdateGroupRadius,
                                       ArrayRef<const RVec>              xGlobal,
                                       bool                              useGpuForNonbonded,
                                       bool                              useGpuForPme,
                                       bool                              useGpuForUpdate,
                                       bool                              useGpuDirectHalo,
                                       bool canUseGpuPmeDecomposition) :
    mdlog_(mdlog), cr_(cr), options_(options), mtop_(mtop), ir_(ir), notifiers_(notifiers)
{
    GMX_LOG(mdlog_.info).appendTextFormatted("\nInitializing Domain Decomposition on %d ranks", cr_->sizeOfDefaultCommunicator);

    ddSettings_ = getDDSettings(mdlog_, options_, mdrunOptions, ir_, useGpuForPme, canUseGpuPmeDecomposition);

    if (ddSettings_.eFlop > 1)
    {
        /* Ensure that we have different random flop counts on different ranks */
        srand(1 + cr_->rankInDefaultCommunicator);
    }

    systemInfo_ = getSystemInfo(mdlog_,
                                MAIN(cr_) ? DDRole::Main : DDRole::Agent,
                                cr->mpiDefaultCommunicator,
                                options_,
                                mtop_,
                                ir_,
                                box,
                                updateGroupingPerMoleculeType,
                                useUpdateGroups,
                                maxUpdateGroupRadius,
                                xGlobal);

    const int  numRanksRequested         = cr_->sizeOfDefaultCommunicator;
    const bool checkForLargePrimeFactors = (options_.numCells[0] <= 0);


    /* Checks for ability to use PME-only ranks */
    auto separatePmeRanksPermitted = checkForSeparatePmeRanks(
            notifiers_, options_, numRanksRequested, useGpuForNonbonded, useGpuForPme, canUseGpuPmeDecomposition);

    /* Checks for validity of requested Ranks setup */
    checkForValidRankCountRequests(numRanksRequested,
                                   usingPme(ir_.coulombtype) || usingLJPme(ir_.vdwtype),
                                   options_.numPmeRanks,
                                   separatePmeRanksPermitted,
                                   checkForLargePrimeFactors);

    // Now that we know whether GPU-direct halos actually will be used, we might have to modify DLB
    if (!isDlbDisabled(ddSettings_.initialDlbState) && useGpuForUpdate && useGpuDirectHalo)
    {
        ddSettings_.initialDlbState = DlbState::offForever;
        GMX_LOG(mdlog.info)
                .appendText(
                        "Disabling dynamic load balancing; unsupported with GPU communication + "
                        "update.");
    }

    // DD grid setup uses a more different cell size limit for
    // automated setup than the one in systemInfo_. The latter is used
    // in set_dd_limits() to configure DLB, for example.
    const real gridSetupCellsizeLimit =
            getDDGridSetupCellSizeLimit(mdlog_,
                                        !isDlbDisabled(ddSettings_.initialDlbState),
                                        options_.dlbScaling,
                                        ir_,
                                        systemInfo_.cellsizeLimit,
                                        numRanksRequested);
    ddGridSetup_ = getDDGridSetup(mdlog_,
                                  MAIN(cr_) ? DDRole::Main : DDRole::Agent,
                                  cr->mpiDefaultCommunicator,
                                  numRanksRequested,
                                  options_,
                                  ddSettings_,
                                  systemInfo_,
                                  gridSetupCellsizeLimit,
                                  mtop_,
                                  ir_,
                                  separatePmeRanksPermitted,
                                  box,
                                  xGlobal,
                                  &ddbox_);
    checkDDGridSetup(ddGridSetup_,
                     MAIN(cr_) ? DDRole::Main : DDRole::Agent,
                     cr->mpiDefaultCommunicator,
                     cr->sizeOfDefaultCommunicator,
                     options_,
                     ddSettings_,
                     systemInfo_,
                     gridSetupCellsizeLimit,
                     ddbox_);

    cr_->npmenodes = ddGridSetup_.numPmeOnlyRanks;

    ddRankSetup_ = getDDRankSetup(
            mdlog_, cr_->sizeOfDefaultCommunicator, options_.rankOrder, ddGridSetup_, ir_);

    /* Generate the group communicator, also decides the duty of each rank */
    cartSetup_ = makeGroupCommunicators(
            mdlog_, ddSettings_, options_.rankOrder, ddRankSetup_, cr_, ddCellIndex_, &pmeRanks_);
}

std::unique_ptr<gmx_domdec_t> DomainDecompositionBuilder::Impl::build(LocalAtomSetManager* atomSets,
                                                                      const gmx_localtop_t& localTopology,
                                                                      const t_state& localState,
                                                                      ObservablesReducerBuilder* observablesReducerBuilder)
{
    auto dd = std::make_unique<gmx_domdec_t>(
            ir_, arrayRefFromArray(ddGridSetup_.ddDimensions, ddGridSetup_.numDDDimensions));

    copy_ivec(ddCellIndex_, dd->ci);

    dd->comm = init_dd_comm();

    dd->comm->ddRankSetup        = ddRankSetup_;
    dd->comm->cartesianRankSetup = cartSetup_;

    set_dd_limits(mdlog_,
                  MAIN(cr_) ? DDRole::Main : DDRole::Agent,
                  dd.get(),
                  options_,
                  ddSettings_,
                  systemInfo_,
                  ddGridSetup_,
                  ddRankSetup_.numPPRanks,
                  mtop_,
                  ir_,
                  ddbox_);

    setupGroupCommunication(mdlog_, ddSettings_, pmeRanks_, cr_, mtop_.natoms, dd.get());

    if (thisRankHasDuty(cr_, DUTY_PP))
    {
        set_ddgrid_parameters(mdlog_, dd.get(), options_.dlbScaling, mtop_, ir_, &ddbox_);

        setup_neighbor_relations(dd.get());
    }

    /* Set overallocation to avoid frequent reallocation of arrays */
    set_over_alloc_dd(true);

    dd->atomSets = atomSets;

    dd->localTopologyChecker = std::make_unique<LocalTopologyChecker>(mdlog_,
                                                                      cr_,
                                                                      mtop_,
                                                                      options_.ddBondedChecking,
                                                                      localTopology,
                                                                      localState,
                                                                      dd->comm->systemInfo.useUpdateGroups,
                                                                      observablesReducerBuilder);

#if GMX_MPI
    MPI_Type_contiguous(DIM, GMX_MPI_REAL, &dd->comm->mpiRVec);
    MPI_Type_commit(&dd->comm->mpiRVec);
#endif

    return dd;
}

DomainDecompositionBuilder::DomainDecompositionBuilder(const MDLogger&           mdlog,
                                                       t_commrec*                cr,
                                                       const DomdecOptions&      options,
                                                       const MdrunOptions&       mdrunOptions,
                                                       const gmx_mtop_t&         mtop,
                                                       const t_inputrec&         ir,
                                                       const MDModulesNotifiers& notifiers,
                                                       const matrix              box,
                                                       ArrayRef<const RangePartitioning> updateGroupingPerMoleculeType,
                                                       const bool           useUpdateGroups,
                                                       const real           maxUpdateGroupRadius,
                                                       ArrayRef<const RVec> xGlobal,
                                                       const bool           useGpuForNonbonded,
                                                       const bool           useGpuForPme,
                                                       bool                 useGpuForUpdate,
                                                       bool                 useGpuDirectHalo,
                                                       const bool canUseGpuPmeDecomposition) :
    impl_(new Impl(mdlog,
                   cr,
                   options,
                   mdrunOptions,
                   mtop,
                   ir,
                   notifiers,
                   box,
                   updateGroupingPerMoleculeType,
                   useUpdateGroups,
                   maxUpdateGroupRadius,
                   xGlobal,
                   useGpuForNonbonded,
                   useGpuForPme,
                   useGpuForUpdate,
                   useGpuDirectHalo,
                   canUseGpuPmeDecomposition))
{
}

std::unique_ptr<gmx_domdec_t> DomainDecompositionBuilder::build(LocalAtomSetManager*  atomSets,
                                                                const gmx_localtop_t& localTopology,
                                                                const t_state&        localState,
                                                                ObservablesReducerBuilder* observablesReducerBuilder)
{
    return impl_->build(atomSets, localTopology, localState, observablesReducerBuilder);
}

DomainDecompositionBuilder::~DomainDecompositionBuilder() = default;

} // namespace gmx


//! Returns the number of halo communication pulses along Cartesian dimension \p dim
static int getNumCommunicationPulsesForDim(const gmx_ddbox_t& ddbox,
                                           const int          dim,
                                           const int          numDomains,
                                           const bool         ddBoxIsDynamic,
                                           const real         communicationDistance)
{
    real inverseOfDomainSize = DD_CELL_MARGIN * numDomains / ddbox.box_size[dim];

    if (ddBoxIsDynamic)
    {
        inverseOfDomainSize *= DD_PRES_SCALE_MARGIN;
    }

    // second part truncates, but since we add 1 this means we return value rounded up.
    return 1 + static_cast<int>(communicationDistance * inverseOfDomainSize * ddbox.skew_fac[dim]);
}

/* Returns whether a cutoff distance of \p cutoffRequested satisfies
 * all limitations of the domain decomposition and thus could be used
 */
static gmx_bool test_dd_cutoff(const t_commrec*               cr,
                               const matrix                   box,
                               gmx::ArrayRef<const gmx::RVec> x,
                               real                           cutoffRequested,
                               bool                           checkGpuDdLimitation)
{
    gmx_ddbox_t ddbox;
    int         LocallyLimited = 0;

    const auto* dd = cr->dd;

    set_ddbox(*dd, false, box, true, x, &ddbox);

    LocallyLimited = 0;

    for (int d = 0; d < dd->ndim; d++)
    {
        const int dim = dd->dim[d];

        const int np = getNumCommunicationPulsesForDim(
                ddbox, dim, dd->numCells[dim], dd->unitCellInfo.ddBoxIsDynamic, cutoffRequested);

        if (!isDlbDisabled(dd->comm->dlbState) && (dim < ddbox.npbcdim) && (dd->comm->cd[d].np_dlb > 0))
        {
            if (np > dd->comm->cd[d].np_dlb)
            {
                return FALSE;
            }

            /* If a current local cell size is smaller than the requested
             * cut-off, we could still fix it, but this gets very complicated.
             * Without fixing here, we might actually need more checks.
             */
            real cellSizeAlongDim =
                    (dd->comm->cell_x1[dim] - dd->comm->cell_x0[dim]) * ddbox.skew_fac[dim];
            if (cellSizeAlongDim * dd->comm->cd[d].np_dlb < cutoffRequested)
            {
                LocallyLimited = 1;
            }
        }

        /* The GPU halo communication code currently does not allow multiple
         * pulses along dimensions other than the first.
         */
        if (checkGpuDdLimitation && (!cr->dd->gpuHaloExchange[0].empty()) && d > 0 && np > 1)
        {
            return FALSE;
        }
    }

    if (!isDlbDisabled(dd->comm->dlbState))
    {
        /* If DLB is not active yet, we don't need to check the grid jumps.
         * Actually we shouldn't, because then the grid jump data is not set.
         */
        if (isDlbOn(dd->comm->dlbState) && gmx::check_grid_jump(0, dd, cutoffRequested, &ddbox, FALSE))
        {
            LocallyLimited = 1;
        }

        gmx_sumi(1, &LocallyLimited, cr);

        if (LocallyLimited > 0)
        {
            return FALSE;
        }
    }

    return TRUE;
}

bool change_dd_cutoff(t_commrec*                     cr,
                      const matrix                   box,
                      gmx::ArrayRef<const gmx::RVec> x,
                      real                           cutoffRequested,
                      bool                           checkGpuDdLimitation)
{
    bool bCutoffAllowed = test_dd_cutoff(cr, box, x, cutoffRequested, checkGpuDdLimitation);

    if (bCutoffAllowed)
    {
        cr->dd->comm->systemInfo.cutoff = cutoffRequested;
    }

    return bCutoffAllowed;
}

void constructGpuHaloExchange(const t_commrec&                cr,
                              const gmx::DeviceStreamManager& deviceStreamManager,
                              gmx_wallcycle*                  wcycle)
{
    GMX_RELEASE_ASSERT(deviceStreamManager.streamIsValid(gmx::DeviceStreamType::NonBondedLocal),
                       "Local non-bonded stream should be valid when using"
                       "GPU halo exchange.");
    GMX_RELEASE_ASSERT(deviceStreamManager.streamIsValid(gmx::DeviceStreamType::NonBondedNonLocal),
                       "Non-local non-bonded stream should be valid when using "
                       "GPU halo exchange.");

    for (int d = 0; d < cr.dd->ndim; d++)
    {
        for (int pulse = cr.dd->gpuHaloExchange[d].size(); pulse < cr.dd->comm->cd[d].numPulses(); pulse++)
        {
            cr.dd->gpuHaloExchange[d].push_back(std::make_unique<gmx::GpuHaloExchange>(
                    cr.dd, d, cr.mpi_comm_mygroup, deviceStreamManager.context(), pulse, wcycle));
        }
    }
}

void reinitGpuHaloExchange(const t_commrec&              cr,
                           const DeviceBuffer<gmx::RVec> d_coordinatesBuffer,
                           const DeviceBuffer<gmx::RVec> d_forcesBuffer)
{
    for (int d = 0; d < cr.dd->ndim; d++)
    {
        for (int pulse = 0; pulse < cr.dd->comm->cd[d].numPulses(); pulse++)
        {
            cr.dd->gpuHaloExchange[d][pulse]->reinitHalo(d_coordinatesBuffer, d_forcesBuffer);
        }
    }
}

GpuEventSynchronizer* communicateGpuHaloCoordinates(const t_commrec&      cr,
                                                    const matrix          box,
                                                    GpuEventSynchronizer* dependencyEvent)
{
    GpuEventSynchronizer* eventPtr = dependencyEvent;
    for (int d = 0; d < cr.dd->ndim; d++)
    {
        for (int pulse = 0; pulse < cr.dd->comm->cd[d].numPulses(); pulse++)
        {
            eventPtr = cr.dd->gpuHaloExchange[d][pulse]->communicateHaloCoordinates(box, eventPtr);
        }
    }
    return eventPtr;
}

void communicateGpuHaloForces(const t_commrec&                                    cr,
                              bool                                                accumulateForces,
                              gmx::FixedCapacityVector<GpuEventSynchronizer*, 2>* dependencyEvents)
{
    for (int d = cr.dd->ndim - 1; d >= 0; d--)
    {
        for (int pulse = cr.dd->comm->cd[d].numPulses() - 1; pulse >= 0; pulse--)
        {
            cr.dd->gpuHaloExchange[d][pulse]->communicateHaloForces(accumulateForces, dependencyEvents);
            dependencyEvents->push_back(cr.dd->gpuHaloExchange[d][pulse]->getForcesReadyOnDeviceEvent());
        }
    }
}

void dd_init_local_state(const gmx_domdec_t& dd, const t_state* state_global, t_state* state_local)
{
    std::array<int, 5> buf;

    if (DDMAIN(dd))
    {
        buf[0] = state_global->flags();
        buf[1] = state_global->ngtc;
        buf[2] = state_global->nnhpres;
        buf[3] = state_global->nhchainlength;
        buf[4] = state_global->dfhist ? state_global->dfhist->nlambda : 0;
    }
    dd_bcast(&dd, buf.size() * sizeof(int), buf.data());

    init_gtc_state(state_local, buf[1], buf[2], buf[3]);
    init_dfhist_state(state_local, buf[4]);
    state_local->setFlags(buf[0]);
}

void putUpdateGroupAtomsInSamePeriodicImage(const gmx_domdec_t&      dd,
                                            const gmx_mtop_t&        mtop,
                                            const matrix             box,
                                            gmx::ArrayRef<gmx::RVec> positions)
{
    int atomOffset = 0;
    for (const gmx_molblock_t& molblock : mtop.molblock)
    {
        const auto& updateGrouping = dd.comm->systemInfo.updateGroupingsPerMoleculeType[molblock.type];

        for (int mol = 0; mol < molblock.nmol; mol++)
        {
            for (int g = 0; g < updateGrouping.numBlocks(); g++)
            {
                const auto& block     = updateGrouping.block(g);
                const int   atomBegin = atomOffset + block.begin();
                const int   atomEnd   = atomOffset + block.end();
                for (int a = atomBegin + 1; a < atomEnd; a++)
                {
                    // Make sure that atoms in the same update group
                    // are in the same periodic image after restarts.
                    for (int d = DIM - 1; d >= 0; d--)
                    {
                        while (positions[a][d] - positions[atomBegin][d] > 0.5_real * box[d][d])
                        {
                            positions[a] -= box[d];
                        }
                        while (positions[a][d] - positions[atomBegin][d] < -0.5_real * box[d][d])
                        {
                            positions[a] += box[d];
                        }
                    }
                }
            }
            atomOffset += updateGrouping.fullRange().end();
        }
    }
}
