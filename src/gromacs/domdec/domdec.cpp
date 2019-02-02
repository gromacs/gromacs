/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019, by the GROMACS development team, led by
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
#include <memory>

#include "gromacs/domdec/collect.h"
#include "gromacs/domdec/dlb.h"
#include "gromacs/domdec/dlbtiming.h"
#include "gromacs/domdec/domdec_network.h"
#include "gromacs/domdec/ga2la.h"
#include "gromacs/domdec/partition.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/hardware/hw_info.h"
#include "gromacs/listed_forces/manage_threading.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/calc_verletbuf.h"
#include "gromacs/mdlib/constr.h"
#include "gromacs/mdlib/constraintrange.h"
#include "gromacs/mdlib/mdrun.h"
#include "gromacs/mdlib/updategroups.h"
#include "gromacs/mdlib/vsite.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pulling/pull.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/mtop_lookup.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/basenetwork.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/logger.h"
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
#include "domdec_vsite.h"
#include "redistribute.h"
#include "utility.h"

static const char *edlbs_names[int(DlbState::nr)] = { "off", "auto", "locked", "on", "on" };

/* The size per atom group of the cggl_flag buffer in gmx_domdec_comm_t */
#define DD_CGIBS 2

/* The flags for the cggl_flag buffer in gmx_domdec_comm_t */
#define DD_FLAG_NRCG  65535
#define DD_FLAG_FW(d) (1<<(16+(d)*2))
#define DD_FLAG_BW(d) (1<<(16+(d)*2+1))

/* The DD zone order */
static const ivec dd_zo[DD_MAXZONE] =
{{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 1, 1}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}};

/* The non-bonded zone-pair setup for domain decomposition
 * The first number is the i-zone, the second number the first j-zone seen by
 * this i-zone, the third number the last+1 j-zone seen by this i-zone.
 * As is, this is for 3D decomposition, where there are 4 i-zones.
 * With 2D decomposition use only the first 2 i-zones and a last+1 j-zone of 4.
 * With 1D decomposition use only the first i-zone and a last+1 j-zone of 2.
 */
static const int
    ddNonbondedZonePairRanges[DD_MAXIZONE][3] = {{0, 0, 8},
                                                 {1, 3, 6},
                                                 {2, 5, 6},
                                                 {3, 5, 7}};


/*
   #define dd_index(n,i) ((((i)[ZZ]*(n)[YY] + (i)[YY])*(n)[XX]) + (i)[XX])

   static void index2xyz(ivec nc,int ind,ivec xyz)
   {
   xyz[XX] = ind % nc[XX];
   xyz[YY] = (ind / nc[XX]) % nc[YY];
   xyz[ZZ] = ind / (nc[YY]*nc[XX]);
   }
 */

static void ddindex2xyz(const ivec nc, int ind, ivec xyz)
{
    xyz[XX] = ind / (nc[YY]*nc[ZZ]);
    xyz[YY] = (ind / nc[ZZ]) % nc[YY];
    xyz[ZZ] = ind % nc[ZZ];
}

static int ddcoord2ddnodeid(gmx_domdec_t *dd, ivec c)
{
    int ddindex;
    int ddnodeid = -1;

    ddindex = dd_index(dd->nc, c);
    if (dd->comm->bCartesianPP_PME)
    {
        ddnodeid = dd->comm->ddindex2ddnodeid[ddindex];
    }
    else if (dd->comm->bCartesianPP)
    {
#if GMX_MPI
        MPI_Cart_rank(dd->mpi_comm_all, c, &ddnodeid);
#endif
    }
    else
    {
        ddnodeid = ddindex;
    }

    return ddnodeid;
}

int ddglatnr(const gmx_domdec_t *dd, int i)
{
    int atnr;

    if (dd == nullptr)
    {
        atnr = i + 1;
    }
    else
    {
        if (i >= dd->comm->atomRanges.numAtomsTotal())
        {
            gmx_fatal(FARGS, "glatnr called with %d, which is larger than the local number of atoms (%d)", i, dd->comm->atomRanges.numAtomsTotal());
        }
        atnr = dd->globalAtomIndices[i] + 1;
    }

    return atnr;
}

t_block *dd_charge_groups_global(gmx_domdec_t *dd)
{
    return &dd->comm->cgs_gl;
}

void dd_store_state(gmx_domdec_t *dd, t_state *state)
{
    int i;

    if (state->ddp_count != dd->ddp_count)
    {
        gmx_incons("The MD state does not match the domain decomposition state");
    }

    state->cg_gl.resize(dd->ncg_home);
    for (i = 0; i < dd->ncg_home; i++)
    {
        state->cg_gl[i] = dd->globalAtomGroupIndices[i];
    }

    state->ddp_count_cg_gl = dd->ddp_count;
}

gmx_domdec_zones_t *domdec_zones(gmx_domdec_t *dd)
{
    return &dd->comm->zones;
}

void dd_get_ns_ranges(const gmx_domdec_t *dd, int icg,
                      int *jcg0, int *jcg1, ivec shift0, ivec shift1)
{
    gmx_domdec_zones_t *zones;
    int                 izone, d, dim;

    zones = &dd->comm->zones;

    izone = 0;
    while (icg >= zones->izone[izone].cg1)
    {
        izone++;
    }

    if (izone == 0)
    {
        *jcg0 = icg;
    }
    else if (izone < zones->nizone)
    {
        *jcg0 = zones->izone[izone].jcg0;
    }
    else
    {
        gmx_fatal(FARGS, "DD icg %d out of range: izone (%d) >= nizone (%d)",
                  icg, izone, zones->nizone);
    }

    *jcg1 = zones->izone[izone].jcg1;

    for (d = 0; d < dd->ndim; d++)
    {
        dim         = dd->dim[d];
        shift0[dim] = zones->izone[izone].shift0[dim];
        shift1[dim] = zones->izone[izone].shift1[dim];
        if (dd->comm->tric_dir[dim] || (isDlbOn(dd->comm) && d > 0))
        {
            /* A conservative approach, this can be optimized */
            shift0[dim] -= 1;
            shift1[dim] += 1;
        }
    }
}

int dd_numHomeAtoms(const gmx_domdec_t &dd)
{
    return dd.comm->atomRanges.numHomeAtoms();
}

int dd_natoms_mdatoms(const gmx_domdec_t *dd)
{
    /* We currently set mdatoms entries for all atoms:
     * local + non-local + communicated for vsite + constraints
     */

    return dd->comm->atomRanges.numAtomsTotal();
}

int dd_natoms_vsite(const gmx_domdec_t *dd)
{
    return dd->comm->atomRanges.end(DDAtomRanges::Type::Vsites);
}

void dd_get_constraint_range(const gmx_domdec_t *dd, int *at_start, int *at_end)
{
    *at_start = dd->comm->atomRanges.start(DDAtomRanges::Type::Constraints);
    *at_end   = dd->comm->atomRanges.end(DDAtomRanges::Type::Constraints);
}

void dd_move_x(gmx_domdec_t             *dd,
               matrix                    box,
               gmx::ArrayRef<gmx::RVec>  x,
               gmx_wallcycle            *wcycle)
{
    wallcycle_start(wcycle, ewcMOVEX);

    int                    nzone, nat_tot;
    gmx_domdec_comm_t     *comm;
    gmx_domdec_comm_dim_t *cd;
    rvec                   shift = {0, 0, 0};
    gmx_bool               bPBC, bScrew;

    comm = dd->comm;

    const gmx::RangePartitioning &atomGrouping = dd->atomGrouping();

    nzone   = 1;
    nat_tot = comm->atomRanges.numHomeAtoms();
    for (int d = 0; d < dd->ndim; d++)
    {
        bPBC   = (dd->ci[dd->dim[d]] == 0);
        bScrew = (bPBC && dd->bScrewPBC && dd->dim[d] == XX);
        if (bPBC)
        {
            copy_rvec(box[dd->dim[d]], shift);
        }
        cd = &comm->cd[d];
        for (const gmx_domdec_ind_t &ind : cd->ind)
        {
            DDBufferAccess<gmx::RVec>  sendBufferAccess(comm->rvecBuffer, ind.nsend[nzone + 1]);
            gmx::ArrayRef<gmx::RVec>  &sendBuffer = sendBufferAccess.buffer;
            int                        n          = 0;
            if (!bPBC)
            {
                for (int g : ind.index)
                {
                    for (int j : atomGrouping.block(g))
                    {
                        sendBuffer[n] = x[j];
                        n++;
                    }
                }
            }
            else if (!bScrew)
            {
                for (int g : ind.index)
                {
                    for (int j : atomGrouping.block(g))
                    {
                        /* We need to shift the coordinates */
                        for (int d = 0; d < DIM; d++)
                        {
                            sendBuffer[n][d] = x[j][d] + shift[d];
                        }
                        n++;
                    }
                }
            }
            else
            {
                for (int g : ind.index)
                {
                    for (int j : atomGrouping.block(g))
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
            }

            DDBufferAccess<gmx::RVec>  receiveBufferAccess(comm->rvecBuffer2, cd->receiveInPlace ? 0 : ind.nrecv[nzone + 1]);

            gmx::ArrayRef<gmx::RVec>   receiveBuffer;
            if (cd->receiveInPlace)
            {
                receiveBuffer = gmx::arrayRefFromArray(x.data() + nat_tot, ind.nrecv[nzone + 1]);
            }
            else
            {
                receiveBuffer = receiveBufferAccess.buffer;
            }
            /* Send and receive the coordinates */
            ddSendrecv(dd, d, dddirBackward,
                       sendBuffer, receiveBuffer);

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
            nat_tot += ind.nrecv[nzone+1];
        }
        nzone += nzone;
    }

    wallcycle_stop(wcycle, ewcMOVEX);
}

void dd_move_f(gmx_domdec_t             *dd,
               gmx::ArrayRef<gmx::RVec>  f,
               rvec                     *fshift,
               gmx_wallcycle            *wcycle)
{
    wallcycle_start(wcycle, ewcMOVEF);

    int                    nzone, nat_tot;
    gmx_domdec_comm_t     *comm;
    gmx_domdec_comm_dim_t *cd;
    ivec                   vis;
    int                    is;
    gmx_bool               bShiftForcesNeedPbc, bScrew;

    comm = dd->comm;

    const gmx::RangePartitioning &atomGrouping = dd->atomGrouping();

    nzone   = comm->zones.n/2;
    nat_tot = comm->atomRanges.end(DDAtomRanges::Type::Zones);
    for (int d = dd->ndim-1; d >= 0; d--)
    {
        /* Only forces in domains near the PBC boundaries need to
           consider PBC in the treatment of fshift */
        bShiftForcesNeedPbc   = (dd->ci[dd->dim[d]] == 0);
        bScrew                = (bShiftForcesNeedPbc && dd->bScrewPBC && dd->dim[d] == XX);
        if (fshift == nullptr && !bScrew)
        {
            bShiftForcesNeedPbc = FALSE;
        }
        /* Determine which shift vector we need */
        clear_ivec(vis);
        vis[dd->dim[d]] = 1;
        is              = IVEC2IS(vis);

        cd = &comm->cd[d];
        for (int p = cd->numPulses() - 1; p >= 0; p--)
        {
            const gmx_domdec_ind_t    &ind  = cd->ind[p];
            DDBufferAccess<gmx::RVec>  receiveBufferAccess(comm->rvecBuffer, ind.nsend[nzone + 1]);
            gmx::ArrayRef<gmx::RVec>  &receiveBuffer = receiveBufferAccess.buffer;

            nat_tot                        -= ind.nrecv[nzone+1];

            DDBufferAccess<gmx::RVec>  sendBufferAccess(comm->rvecBuffer2, cd->receiveInPlace ? 0 : ind.nrecv[nzone + 1]);

            gmx::ArrayRef<gmx::RVec>   sendBuffer;
            if (cd->receiveInPlace)
            {
                sendBuffer = gmx::arrayRefFromArray(f.data() + nat_tot, ind.nrecv[nzone + 1]);
            }
            else
            {
                sendBuffer = sendBufferAccess.buffer;
                int j = 0;
                for (int zone = 0; zone < nzone; zone++)
                {
                    for (int i = ind.cell2at0[zone]; i < ind.cell2at1[zone]; i++)
                    {
                        sendBuffer[j++] = f[i];
                    }
                }
            }
            /* Communicate the forces */
            ddSendrecv(dd, d, dddirForward,
                       sendBuffer, receiveBuffer);
            /* Add the received forces */
            int n = 0;
            if (!bShiftForcesNeedPbc)
            {
                for (int g : ind.index)
                {
                    for (int j : atomGrouping.block(g))
                    {
                        for (int d = 0; d < DIM; d++)
                        {
                            f[j][d] += receiveBuffer[n][d];
                        }
                        n++;
                    }
                }
            }
            else if (!bScrew)
            {
                /* fshift should always be defined if this function is
                 * called when bShiftForcesNeedPbc is true */
                assert(nullptr != fshift);
                for (int g : ind.index)
                {
                    for (int j : atomGrouping.block(g))
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
            }
            else
            {
                for (int g : ind.index)
                {
                    for (int j : atomGrouping.block(g))
                    {
                        /* Rotate the force */
                        f[j][XX] += receiveBuffer[n][XX];
                        f[j][YY] -= receiveBuffer[n][YY];
                        f[j][ZZ] -= receiveBuffer[n][ZZ];
                        if (fshift)
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
        }
        nzone /= 2;
    }
    wallcycle_stop(wcycle, ewcMOVEF);
}

/* Convenience function for extracting a real buffer from an rvec buffer
 *
 * To reduce the number of temporary communication buffers and avoid
 * cache polution, we reuse gmx::RVec buffers for storing reals.
 * This functions return a real buffer reference with the same number
 * of elements as the gmx::RVec buffer (so 1/3 of the size in bytes).
 */
static gmx::ArrayRef<real>
realArrayRefFromRvecArrayRef(gmx::ArrayRef<gmx::RVec> arrayRef)
{
    return gmx::arrayRefFromArray(as_rvec_array(arrayRef.data())[0],
                                  arrayRef.size());
}

void dd_atom_spread_real(gmx_domdec_t *dd, real v[])
{
    int                    nzone, nat_tot;
    gmx_domdec_comm_t     *comm;
    gmx_domdec_comm_dim_t *cd;

    comm = dd->comm;

    const gmx::RangePartitioning &atomGrouping = dd->atomGrouping();

    nzone   = 1;
    nat_tot = comm->atomRanges.numHomeAtoms();
    for (int d = 0; d < dd->ndim; d++)
    {
        cd = &comm->cd[d];
        for (const gmx_domdec_ind_t &ind : cd->ind)
        {
            /* Note: We provision for RVec instead of real, so a factor of 3
             * more than needed. The buffer actually already has this size
             * and we pass a plain pointer below, so this does not matter.
             */
            DDBufferAccess<gmx::RVec> sendBufferAccess(comm->rvecBuffer, ind.nsend[nzone + 1]);
            gmx::ArrayRef<real>       sendBuffer = realArrayRefFromRvecArrayRef(sendBufferAccess.buffer);
            int                       n          = 0;
            for (int g : ind.index)
            {
                for (int j : atomGrouping.block(g))
                {
                    sendBuffer[n++] = v[j];
                }
            }

            DDBufferAccess<gmx::RVec> receiveBufferAccess(comm->rvecBuffer2, cd->receiveInPlace ? 0 : ind.nrecv[nzone + 1]);

            gmx::ArrayRef<real>       receiveBuffer;
            if (cd->receiveInPlace)
            {
                receiveBuffer = gmx::arrayRefFromArray(v + nat_tot, ind.nrecv[nzone + 1]);
            }
            else
            {
                receiveBuffer = realArrayRefFromRvecArrayRef(receiveBufferAccess.buffer);
            }
            /* Send and receive the data */
            ddSendrecv(dd, d, dddirBackward,
                       sendBuffer, receiveBuffer);
            if (!cd->receiveInPlace)
            {
                int j = 0;
                for (int zone = 0; zone < nzone; zone++)
                {
                    for (int i = ind.cell2at0[zone]; i < ind.cell2at1[zone]; i++)
                    {
                        v[i] = receiveBuffer[j++];
                    }
                }
            }
            nat_tot += ind.nrecv[nzone+1];
        }
        nzone += nzone;
    }
}

void dd_atom_sum_real(gmx_domdec_t *dd, real v[])
{
    int                    nzone, nat_tot;
    gmx_domdec_comm_t     *comm;
    gmx_domdec_comm_dim_t *cd;

    comm = dd->comm;

    const gmx::RangePartitioning &atomGrouping = dd->atomGrouping();

    nzone   = comm->zones.n/2;
    nat_tot = comm->atomRanges.end(DDAtomRanges::Type::Zones);
    for (int d = dd->ndim-1; d >= 0; d--)
    {
        cd = &comm->cd[d];
        for (int p = cd->numPulses() - 1; p >= 0; p--)
        {
            const gmx_domdec_ind_t &ind = cd->ind[p];

            /* Note: We provision for RVec instead of real, so a factor of 3
             * more than needed. The buffer actually already has this size
             * and we typecast, so this works as intended.
             */
            DDBufferAccess<gmx::RVec> receiveBufferAccess(comm->rvecBuffer, ind.nsend[nzone + 1]);
            gmx::ArrayRef<real>       receiveBuffer = realArrayRefFromRvecArrayRef(receiveBufferAccess.buffer);
            nat_tot -= ind.nrecv[nzone + 1];

            DDBufferAccess<gmx::RVec> sendBufferAccess(comm->rvecBuffer2, cd->receiveInPlace ? 0 : ind.nrecv[nzone + 1]);

            gmx::ArrayRef<real>       sendBuffer;
            if (cd->receiveInPlace)
            {
                sendBuffer = gmx::arrayRefFromArray(v + nat_tot, ind.nrecv[nzone + 1]);
            }
            else
            {
                sendBuffer = realArrayRefFromRvecArrayRef(sendBufferAccess.buffer);
                int j = 0;
                for (int zone = 0; zone < nzone; zone++)
                {
                    for (int i = ind.cell2at0[zone]; i < ind.cell2at1[zone]; i++)
                    {
                        sendBuffer[j++] = v[i];
                    }
                }
            }
            /* Communicate the forces */
            ddSendrecv(dd, d, dddirForward,
                       sendBuffer, receiveBuffer);
            /* Add the received forces */
            int n = 0;
            for (int g : ind.index)
            {
                for (int j : atomGrouping.block(g))
                {
                    v[j] += receiveBuffer[n];
                    n++;
                }
            }
        }
        nzone /= 2;
    }
}

real dd_cutoff_multibody(const gmx_domdec_t *dd)
{
    gmx_domdec_comm_t *comm;
    int                di;
    real               r;

    comm = dd->comm;

    r = -1;
    if (comm->bInterCGBondeds)
    {
        if (comm->cutoff_mbody > 0)
        {
            r = comm->cutoff_mbody;
        }
        else
        {
            /* cutoff_mbody=0 means we do not have DLB */
            r = comm->cellsize_min[dd->dim[0]];
            for (di = 1; di < dd->ndim; di++)
            {
                r = std::min(r, comm->cellsize_min[dd->dim[di]]);
            }
            if (comm->bBondComm)
            {
                r = std::max(r, comm->cutoff_mbody);
            }
            else
            {
                r = std::min(r, comm->cutoff);
            }
        }
    }

    return r;
}

real dd_cutoff_twobody(const gmx_domdec_t *dd)
{
    real r_mb;

    r_mb = dd_cutoff_multibody(dd);

    return std::max(dd->comm->cutoff, r_mb);
}


static void dd_cart_coord2pmecoord(const gmx_domdec_t *dd, const ivec coord,
                                   ivec coord_pme)
{
    int nc, ntot;

    nc   = dd->nc[dd->comm->cartpmedim];
    ntot = dd->comm->ntot[dd->comm->cartpmedim];
    copy_ivec(coord, coord_pme);
    coord_pme[dd->comm->cartpmedim] =
        nc + (coord[dd->comm->cartpmedim]*(ntot - nc) + (ntot - nc)/2)/nc;
}

static int ddindex2pmeindex(const gmx_domdec_t *dd, int ddindex)
{
    int npp, npme;

    npp  = dd->nnodes;
    npme = dd->comm->npmenodes;

    /* Here we assign a PME node to communicate with this DD node
     * by assuming that the major index of both is x.
     * We add cr->npmenodes/2 to obtain an even distribution.
     */
    return (ddindex*npme + npme/2)/npp;
}

static int *dd_interleaved_pme_ranks(const gmx_domdec_t *dd)
{
    int *pme_rank;
    int  n, i, p0, p1;

    snew(pme_rank, dd->comm->npmenodes);
    n = 0;
    for (i = 0; i < dd->nnodes; i++)
    {
        p0 = ddindex2pmeindex(dd, i);
        p1 = ddindex2pmeindex(dd, i+1);
        if (i+1 == dd->nnodes || p1 > p0)
        {
            if (debug)
            {
                fprintf(debug, "pme_rank[%d] = %d\n", n, i+1+n);
            }
            pme_rank[n] = i + 1 + n;
            n++;
        }
    }

    return pme_rank;
}

static int gmx_ddcoord2pmeindex(const t_commrec *cr, int x, int y, int z)
{
    gmx_domdec_t *dd;
    ivec          coords;
    int           slab;

    dd = cr->dd;
    /*
       if (dd->comm->bCartesian) {
       gmx_ddindex2xyz(dd->nc,ddindex,coords);
       dd_coords2pmecoords(dd,coords,coords_pme);
       copy_ivec(dd->ntot,nc);
       nc[dd->cartpmedim]         -= dd->nc[dd->cartpmedim];
       coords_pme[dd->cartpmedim] -= dd->nc[dd->cartpmedim];

       slab = (coords_pme[XX]*nc[YY] + coords_pme[YY])*nc[ZZ] + coords_pme[ZZ];
       } else {
       slab = (ddindex*cr->npmenodes + cr->npmenodes/2)/dd->nnodes;
       }
     */
    coords[XX] = x;
    coords[YY] = y;
    coords[ZZ] = z;
    slab       = ddindex2pmeindex(dd, dd_index(dd->nc, coords));

    return slab;
}

static int ddcoord2simnodeid(const t_commrec *cr, int x, int y, int z)
{
    gmx_domdec_comm_t *comm;
    ivec               coords;
    int                ddindex, nodeid = -1;

    comm = cr->dd->comm;

    coords[XX] = x;
    coords[YY] = y;
    coords[ZZ] = z;
    if (comm->bCartesianPP_PME)
    {
#if GMX_MPI
        MPI_Cart_rank(cr->mpi_comm_mysim, coords, &nodeid);
#endif
    }
    else
    {
        ddindex = dd_index(cr->dd->nc, coords);
        if (comm->bCartesianPP)
        {
            nodeid = comm->ddindex2simnodeid[ddindex];
        }
        else
        {
            if (comm->pmenodes)
            {
                nodeid = ddindex + gmx_ddcoord2pmeindex(cr, x, y, z);
            }
            else
            {
                nodeid = ddindex;
            }
        }
    }

    return nodeid;
}

static int dd_simnode2pmenode(const gmx_domdec_t         *dd,
                              const t_commrec gmx_unused *cr,
                              int                         sim_nodeid)
{
    int pmenode = -1;

    const gmx_domdec_comm_t *comm = dd->comm;

    /* This assumes a uniform x domain decomposition grid cell size */
    if (comm->bCartesianPP_PME)
    {
#if GMX_MPI
        ivec coord, coord_pme;
        MPI_Cart_coords(cr->mpi_comm_mysim, sim_nodeid, DIM, coord);
        if (coord[comm->cartpmedim] < dd->nc[comm->cartpmedim])
        {
            /* This is a PP node */
            dd_cart_coord2pmecoord(dd, coord, coord_pme);
            MPI_Cart_rank(cr->mpi_comm_mysim, coord_pme, &pmenode);
        }
#endif
    }
    else if (comm->bCartesianPP)
    {
        if (sim_nodeid < dd->nnodes)
        {
            pmenode = dd->nnodes + ddindex2pmeindex(dd, sim_nodeid);
        }
    }
    else
    {
        /* This assumes DD cells with identical x coordinates
         * are numbered sequentially.
         */
        if (dd->comm->pmenodes == nullptr)
        {
            if (sim_nodeid < dd->nnodes)
            {
                /* The DD index equals the nodeid */
                pmenode = dd->nnodes + ddindex2pmeindex(dd, sim_nodeid);
            }
        }
        else
        {
            int i = 0;
            while (sim_nodeid > dd->comm->pmenodes[i])
            {
                i++;
            }
            if (sim_nodeid < dd->comm->pmenodes[i])
            {
                pmenode = dd->comm->pmenodes[i];
            }
        }
    }

    return pmenode;
}

NumPmeDomains getNumPmeDomains(const gmx_domdec_t *dd)
{
    if (dd != nullptr)
    {
        return {
                   dd->comm->npmenodes_x, dd->comm->npmenodes_y
        };
    }
    else
    {
        return {
                   1, 1
        };
    }
}

std::vector<int> get_pme_ddranks(const t_commrec *cr, int pmenodeid)
{
    gmx_domdec_t *dd;
    int           x, y, z;
    ivec          coord, coord_pme;

    dd = cr->dd;

    std::vector<int> ddranks;
    ddranks.reserve((dd->nnodes+cr->npmenodes-1)/cr->npmenodes);

    for (x = 0; x < dd->nc[XX]; x++)
    {
        for (y = 0; y < dd->nc[YY]; y++)
        {
            for (z = 0; z < dd->nc[ZZ]; z++)
            {
                if (dd->comm->bCartesianPP_PME)
                {
                    coord[XX] = x;
                    coord[YY] = y;
                    coord[ZZ] = z;
                    dd_cart_coord2pmecoord(dd, coord, coord_pme);
                    if (dd->ci[XX] == coord_pme[XX] &&
                        dd->ci[YY] == coord_pme[YY] &&
                        dd->ci[ZZ] == coord_pme[ZZ])
                    {
                        ddranks.push_back(ddcoord2simnodeid(cr, x, y, z));
                    }
                }
                else
                {
                    /* The slab corresponds to the nodeid in the PME group */
                    if (gmx_ddcoord2pmeindex(cr, x, y, z) == pmenodeid)
                    {
                        ddranks.push_back(ddcoord2simnodeid(cr, x, y, z));
                    }
                }
            }
        }
    }
    return ddranks;
}

static gmx_bool receive_vir_ener(const gmx_domdec_t *dd, const t_commrec *cr)
{
    gmx_bool bReceive = TRUE;

    if (cr->npmenodes < dd->nnodes)
    {
        gmx_domdec_comm_t *comm = dd->comm;
        if (comm->bCartesianPP_PME)
        {
#if GMX_MPI
            int  pmenode = dd_simnode2pmenode(dd, cr, cr->sim_nodeid);
            ivec coords;
            MPI_Cart_coords(cr->mpi_comm_mysim, cr->sim_nodeid, DIM, coords);
            coords[comm->cartpmedim]++;
            if (coords[comm->cartpmedim] < dd->nc[comm->cartpmedim])
            {
                int rank;
                MPI_Cart_rank(cr->mpi_comm_mysim, coords, &rank);
                if (dd_simnode2pmenode(dd, cr, rank) == pmenode)
                {
                    /* This is not the last PP node for pmenode */
                    bReceive = FALSE;
                }
            }
#else
            GMX_RELEASE_ASSERT(false, "Without MPI we should not have Cartesian PP-PME with #PMEnodes < #DDnodes");
#endif
        }
        else
        {
            int pmenode = dd_simnode2pmenode(dd, cr, cr->sim_nodeid);
            if (cr->sim_nodeid+1 < cr->nnodes &&
                dd_simnode2pmenode(dd, cr, cr->sim_nodeid+1) == pmenode)
            {
                /* This is not the last PP node for pmenode */
                bReceive = FALSE;
            }
        }
    }

    return bReceive;
}

static void set_slb_pme_dim_f(gmx_domdec_t *dd, int dim, real **dim_f)
{
    gmx_domdec_comm_t *comm;
    int                i;

    comm = dd->comm;

    snew(*dim_f, dd->nc[dim]+1);
    (*dim_f)[0] = 0;
    for (i = 1; i < dd->nc[dim]; i++)
    {
        if (comm->slb_frac[dim])
        {
            (*dim_f)[i] = (*dim_f)[i-1] + comm->slb_frac[dim][i-1];
        }
        else
        {
            (*dim_f)[i] = static_cast<real>(i)/static_cast<real>(dd->nc[dim]);
        }
    }
    (*dim_f)[dd->nc[dim]] = 1;
}

static void init_ddpme(gmx_domdec_t *dd, gmx_ddpme_t *ddpme, int dimind)
{
    int  pmeindex, slab, nso, i;
    ivec xyz;

    if (dimind == 0 && dd->dim[0] == YY && dd->comm->npmenodes_x == 1)
    {
        ddpme->dim = YY;
    }
    else
    {
        ddpme->dim = dimind;
    }
    ddpme->dim_match = (ddpme->dim == dd->dim[dimind]);

    ddpme->nslab = (ddpme->dim == 0 ?
                    dd->comm->npmenodes_x :
                    dd->comm->npmenodes_y);

    if (ddpme->nslab <= 1)
    {
        return;
    }

    nso = dd->comm->npmenodes/ddpme->nslab;
    /* Determine for each PME slab the PP location range for dimension dim */
    snew(ddpme->pp_min, ddpme->nslab);
    snew(ddpme->pp_max, ddpme->nslab);
    for (slab = 0; slab < ddpme->nslab; slab++)
    {
        ddpme->pp_min[slab] = dd->nc[dd->dim[dimind]] - 1;
        ddpme->pp_max[slab] = 0;
    }
    for (i = 0; i < dd->nnodes; i++)
    {
        ddindex2xyz(dd->nc, i, xyz);
        /* For y only use our y/z slab.
         * This assumes that the PME x grid size matches the DD grid size.
         */
        if (dimind == 0 || xyz[XX] == dd->ci[XX])
        {
            pmeindex = ddindex2pmeindex(dd, i);
            if (dimind == 0)
            {
                slab = pmeindex/nso;
            }
            else
            {
                slab = pmeindex % ddpme->nslab;
            }
            ddpme->pp_min[slab] = std::min(ddpme->pp_min[slab], xyz[dimind]);
            ddpme->pp_max[slab] = std::max(ddpme->pp_max[slab], xyz[dimind]);
        }
    }

    set_slb_pme_dim_f(dd, ddpme->dim, &ddpme->slb_dim_f);
}

int dd_pme_maxshift_x(const gmx_domdec_t *dd)
{
    if (dd->comm->ddpme[0].dim == XX)
    {
        return dd->comm->ddpme[0].maxshift;
    }
    else
    {
        return 0;
    }
}

int dd_pme_maxshift_y(const gmx_domdec_t *dd)
{
    if (dd->comm->ddpme[0].dim == YY)
    {
        return dd->comm->ddpme[0].maxshift;
    }
    else if (dd->comm->npmedecompdim >= 2 && dd->comm->ddpme[1].dim == YY)
    {
        return dd->comm->ddpme[1].maxshift;
    }
    else
    {
        return 0;
    }
}

void dd_cycles_add(const gmx_domdec_t *dd, float cycles, int ddCycl)
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
static void make_load_communicator(gmx_domdec_t *dd, int dim_ind, ivec loc)
{
    MPI_Comm           c_row;
    int                dim, i, rank;
    ivec               loc_c;
    gmx_bool           bPartOfGroup = FALSE;

    dim = dd->dim[dim_ind];
    copy_ivec(loc, loc_c);
    for (i = 0; i < dd->nc[dim]; i++)
    {
        loc_c[dim] = i;
        rank       = dd_index(dd->nc, loc_c);
        if (rank == dd->rank)
        {
            /* This process is part of the group */
            bPartOfGroup = TRUE;
        }
    }
    MPI_Comm_split(dd->mpi_comm_all, bPartOfGroup ? 0 : MPI_UNDEFINED, dd->rank,
                   &c_row);
    if (bPartOfGroup)
    {
        dd->comm->mpi_comm_load[dim_ind] = c_row;
        if (!isDlbDisabled(dd->comm))
        {
            DDCellsizesWithDlb &cellsizes = dd->comm->cellsizesWithDlb[dim_ind];

            if (dd->ci[dim] == dd->master_ci[dim])
            {
                /* This is the root process of this row */
                cellsizes.rowMaster  = std::make_unique<RowMaster>();

                RowMaster &rowMaster = *cellsizes.rowMaster;
                rowMaster.cellFrac.resize(ddCellFractionBufferSize(dd, dim_ind));
                rowMaster.oldCellFrac.resize(dd->nc[dim] + 1);
                rowMaster.isCellMin.resize(dd->nc[dim]);
                if (dim_ind > 0)
                {
                    rowMaster.bounds.resize(dd->nc[dim]);
                }
                rowMaster.buf_ncd.resize(dd->nc[dim]);
            }
            else
            {
                /* This is not a root process, we only need to receive cell_f */
                cellsizes.fracRow.resize(ddCellFractionBufferSize(dd, dim_ind));
            }
        }
        if (dd->ci[dim] == dd->master_ci[dim])
        {
            snew(dd->comm->load[dim_ind].load, dd->nc[dim]*DD_NLOAD_MAX);
        }
    }
}
#endif

void dd_setup_dlb_resource_sharing(t_commrec            *cr,
                                   int                   gpu_id)
{
#if GMX_MPI
    int           physicalnode_id_hash;
    gmx_domdec_t *dd;
    MPI_Comm      mpi_comm_pp_physicalnode;

    if (!thisRankHasDuty(cr, DUTY_PP) || gpu_id < 0)
    {
        /* Only ranks with short-ranged tasks (currently) use GPUs.
         * If we don't have GPUs assigned, there are no resources to share.
         */
        return;
    }

    physicalnode_id_hash = gmx_physicalnode_id_hash();

    dd = cr->dd;

    if (debug)
    {
        fprintf(debug, "dd_setup_dd_dlb_gpu_sharing:\n");
        fprintf(debug, "DD PP rank %d physical node hash %d gpu_id %d\n",
                dd->rank, physicalnode_id_hash, gpu_id);
    }
    /* Split the PP communicator over the physical nodes */
    /* TODO: See if we should store this (before), as it's also used for
     * for the nodecomm summation.
     */
    // TODO PhysicalNodeCommunicator could be extended/used to handle
    // the need for per-node per-group communicators.
    MPI_Comm_split(dd->mpi_comm_all, physicalnode_id_hash, dd->rank,
                   &mpi_comm_pp_physicalnode);
    MPI_Comm_split(mpi_comm_pp_physicalnode, gpu_id, dd->rank,
                   &dd->comm->mpi_comm_gpu_shared);
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

static void make_load_communicators(gmx_domdec_t gmx_unused *dd)
{
#if GMX_MPI
    int  dim0, dim1, i, j;
    ivec loc;

    if (debug)
    {
        fprintf(debug, "Making load communicators\n");
    }

    snew(dd->comm->load,          std::max(dd->ndim, 1));
    snew(dd->comm->mpi_comm_load, std::max(dd->ndim, 1));

    if (dd->ndim == 0)
    {
        return;
    }

    clear_ivec(loc);
    make_load_communicator(dd, 0, loc);
    if (dd->ndim > 1)
    {
        dim0 = dd->dim[0];
        for (i = 0; i < dd->nc[dim0]; i++)
        {
            loc[dim0] = i;
            make_load_communicator(dd, 1, loc);
        }
    }
    if (dd->ndim > 2)
    {
        dim0 = dd->dim[0];
        for (i = 0; i < dd->nc[dim0]; i++)
        {
            loc[dim0] = i;
            dim1      = dd->dim[1];
            for (j = 0; j < dd->nc[dim1]; j++)
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
static void setup_neighbor_relations(gmx_domdec_t *dd)
{
    int                     d, dim, i, j, m;
    ivec                    tmp, s;
    gmx_domdec_zones_t     *zones;
    gmx_domdec_ns_ranges_t *izone;

    for (d = 0; d < dd->ndim; d++)
    {
        dim = dd->dim[d];
        copy_ivec(dd->ci, tmp);
        tmp[dim]           = (tmp[dim] + 1) % dd->nc[dim];
        dd->neighbor[d][0] = ddcoord2ddnodeid(dd, tmp);
        copy_ivec(dd->ci, tmp);
        tmp[dim]           = (tmp[dim] - 1 + dd->nc[dim]) % dd->nc[dim];
        dd->neighbor[d][1] = ddcoord2ddnodeid(dd, tmp);
        if (debug)
        {
            fprintf(debug, "DD rank %d neighbor ranks in dir %d are + %d - %d\n",
                    dd->rank, dim,
                    dd->neighbor[d][0],
                    dd->neighbor[d][1]);
        }
    }

    int nzone  = (1 << dd->ndim);
    int nizone = (1 << std::max(dd->ndim - 1, 0));
    assert(nizone >= 1 && nizone <= DD_MAXIZONE);

    zones = &dd->comm->zones;

    for (i = 0; i < nzone; i++)
    {
        m = 0;
        clear_ivec(zones->shift[i]);
        for (d = 0; d < dd->ndim; d++)
        {
            zones->shift[i][dd->dim[d]] = dd_zo[i][m++];
        }
    }

    zones->n = nzone;
    for (i = 0; i < nzone; i++)
    {
        for (d = 0; d < DIM; d++)
        {
            s[d] = dd->ci[d] - zones->shift[i][d];
            if (s[d] < 0)
            {
                s[d] += dd->nc[d];
            }
            else if (s[d] >= dd->nc[d])
            {
                s[d] -= dd->nc[d];
            }
        }
    }
    zones->nizone = nizone;
    for (i = 0; i < zones->nizone; i++)
    {
        assert(ddNonbondedZonePairRanges[i][0] == i);

        izone     = &zones->izone[i];
        /* dd_zp3 is for 3D decomposition, for fewer dimensions use only
         * j-zones up to nzone.
         */
        izone->j0 = std::min(ddNonbondedZonePairRanges[i][1], nzone);
        izone->j1 = std::min(ddNonbondedZonePairRanges[i][2], nzone);
        for (dim = 0; dim < DIM; dim++)
        {
            if (dd->nc[dim] == 1)
            {
                /* All shifts should be allowed */
                izone->shift0[dim] = -1;
                izone->shift1[dim] = 1;
            }
            else
            {
                /* Determine the min/max j-zone shift wrt the i-zone */
                izone->shift0[dim] = 1;
                izone->shift1[dim] = -1;
                for (j = izone->j0; j < izone->j1; j++)
                {
                    int shift_diff = zones->shift[j][dim] - zones->shift[i][dim];
                    if (shift_diff < izone->shift0[dim])
                    {
                        izone->shift0[dim] = shift_diff;
                    }
                    if (shift_diff > izone->shift1[dim])
                    {
                        izone->shift1[dim] = shift_diff;
                    }
                }
            }
        }
    }

    if (!isDlbDisabled(dd->comm))
    {
        dd->comm->cellsizesWithDlb.resize(dd->ndim);
    }

    if (dd->comm->bRecordLoad)
    {
        make_load_communicators(dd);
    }
}

static void make_pp_communicator(const gmx::MDLogger  &mdlog,
                                 gmx_domdec_t         *dd,
                                 t_commrec gmx_unused *cr,
                                 bool gmx_unused       reorder)
{
#if GMX_MPI
    gmx_domdec_comm_t *comm;
    int                rank, *buf;
    ivec               periods;
    MPI_Comm           comm_cart;

    comm = dd->comm;

    if (comm->bCartesianPP)
    {
        /* Set up cartesian communication for the particle-particle part */
        GMX_LOG(mdlog.info).appendTextFormatted(
                "Will use a Cartesian communicator: %d x %d x %d",
                dd->nc[XX], dd->nc[YY], dd->nc[ZZ]);

        for (int i = 0; i < DIM; i++)
        {
            periods[i] = TRUE;
        }
        MPI_Cart_create(cr->mpi_comm_mygroup, DIM, dd->nc, periods, static_cast<int>(reorder),
                        &comm_cart);
        /* We overwrite the old communicator with the new cartesian one */
        cr->mpi_comm_mygroup = comm_cart;
    }

    dd->mpi_comm_all = cr->mpi_comm_mygroup;
    MPI_Comm_rank(dd->mpi_comm_all, &dd->rank);

    if (comm->bCartesianPP_PME)
    {
        /* Since we want to use the original cartesian setup for sim,
         * and not the one after split, we need to make an index.
         */
        snew(comm->ddindex2ddnodeid, dd->nnodes);
        comm->ddindex2ddnodeid[dd_index(dd->nc, dd->ci)] = dd->rank;
        gmx_sumi(dd->nnodes, comm->ddindex2ddnodeid, cr);
        /* Get the rank of the DD master,
         * above we made sure that the master node is a PP node.
         */
        if (MASTER(cr))
        {
            rank = dd->rank;
        }
        else
        {
            rank = 0;
        }
        MPI_Allreduce(&rank, &dd->masterrank, 1, MPI_INT, MPI_SUM, dd->mpi_comm_all);
    }
    else if (comm->bCartesianPP)
    {
        if (cr->npmenodes == 0)
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
        snew(comm->ddindex2simnodeid, dd->nnodes);
        snew(buf, dd->nnodes);
        if (thisRankHasDuty(cr, DUTY_PP))
        {
            buf[dd_index(dd->nc, dd->ci)] = cr->sim_nodeid;
        }
        /* Communicate the ddindex to simulation nodeid index */
        MPI_Allreduce(buf, comm->ddindex2simnodeid, dd->nnodes, MPI_INT, MPI_SUM,
                      cr->mpi_comm_mysim);
        sfree(buf);

        /* Determine the master coordinates and rank.
         * The DD master should be the same node as the master of this sim.
         */
        for (int i = 0; i < dd->nnodes; i++)
        {
            if (comm->ddindex2simnodeid[i] == 0)
            {
                ddindex2xyz(dd->nc, i, dd->master_ci);
                MPI_Cart_rank(dd->mpi_comm_all, dd->master_ci, &dd->masterrank);
            }
        }
        if (debug)
        {
            fprintf(debug, "The master rank is %d\n", dd->masterrank);
        }
    }
    else
    {
        /* No Cartesian communicators */
        /* We use the rank in dd->comm->all as DD index */
        ddindex2xyz(dd->nc, dd->rank, dd->ci);
        /* The simulation master nodeid is 0, so the DD master rank is also 0 */
        dd->masterrank = 0;
        clear_ivec(dd->master_ci);
    }
#endif

    GMX_LOG(mdlog.info).appendTextFormatted(
            "Domain decomposition rank %d, coordinates %d %d %d\n",
            dd->rank, dd->ci[XX], dd->ci[YY], dd->ci[ZZ]);
    if (debug)
    {
        fprintf(debug,
                "Domain decomposition rank %d, coordinates %d %d %d\n\n",
                dd->rank, dd->ci[XX], dd->ci[YY], dd->ci[ZZ]);
    }
}

static void receive_ddindex2simnodeid(gmx_domdec_t         *dd,
                                      t_commrec            *cr)
{
#if GMX_MPI
    gmx_domdec_comm_t *comm = dd->comm;

    if (!comm->bCartesianPP_PME && comm->bCartesianPP)
    {
        int *buf;
        snew(comm->ddindex2simnodeid, dd->nnodes);
        snew(buf, dd->nnodes);
        if (thisRankHasDuty(cr, DUTY_PP))
        {
            buf[dd_index(dd->nc, dd->ci)] = cr->sim_nodeid;
        }
        /* Communicate the ddindex to simulation nodeid index */
        MPI_Allreduce(buf, comm->ddindex2simnodeid, dd->nnodes, MPI_INT, MPI_SUM,
                      cr->mpi_comm_mysim);
        sfree(buf);
    }
#else
    GMX_UNUSED_VALUE(dd);
    GMX_UNUSED_VALUE(cr);
#endif
}

static void split_communicator(const gmx::MDLogger &mdlog,
                               t_commrec *cr, gmx_domdec_t *dd,
                               DdRankOrder gmx_unused rankOrder,
                               bool gmx_unused reorder)
{
    gmx_domdec_comm_t *comm;
    int                i;
    gmx_bool           bDiv[DIM];
#if GMX_MPI
    MPI_Comm           comm_cart;
#endif

    comm = dd->comm;

    if (comm->bCartesianPP)
    {
        for (i = 1; i < DIM; i++)
        {
            bDiv[i] = ((cr->npmenodes*dd->nc[i]) % (dd->nnodes) == 0);
        }
        if (bDiv[YY] || bDiv[ZZ])
        {
            comm->bCartesianPP_PME = TRUE;
            /* If we have 2D PME decomposition, which is always in x+y,
             * we stack the PME only nodes in z.
             * Otherwise we choose the direction that provides the thinnest slab
             * of PME only nodes as this will have the least effect
             * on the PP communication.
             * But for the PME communication the opposite might be better.
             */
            if (bDiv[ZZ] && (comm->npmenodes_y > 1 ||
                             !bDiv[YY] ||
                             dd->nc[YY] > dd->nc[ZZ]))
            {
                comm->cartpmedim = ZZ;
            }
            else
            {
                comm->cartpmedim = YY;
            }
            comm->ntot[comm->cartpmedim]
                += (cr->npmenodes*dd->nc[comm->cartpmedim])/dd->nnodes;
        }
        else
        {
            GMX_LOG(mdlog.info).appendTextFormatted(
                    "Number of PME-only ranks (%d) is not a multiple of nx*ny (%d*%d) or nx*nz (%d*%d)",
                    cr->npmenodes, dd->nc[XX], dd->nc[YY], dd->nc[XX], dd->nc[ZZ]);
            GMX_LOG(mdlog.info).appendText("Will not use a Cartesian communicator for PP <-> PME\n");
        }
    }

    if (comm->bCartesianPP_PME)
    {
#if GMX_MPI
        int  rank;
        ivec periods;

        GMX_LOG(mdlog.info).appendTextFormatted(
                "Will use a Cartesian communicator for PP <-> PME: %d x %d x %d",
                comm->ntot[XX], comm->ntot[YY], comm->ntot[ZZ]);

        for (i = 0; i < DIM; i++)
        {
            periods[i] = TRUE;
        }
        MPI_Cart_create(cr->mpi_comm_mysim, DIM, comm->ntot, periods, static_cast<int>(reorder),
                        &comm_cart);
        MPI_Comm_rank(comm_cart, &rank);
        if (MASTER(cr) && rank != 0)
        {
            gmx_fatal(FARGS, "MPI rank 0 was renumbered by MPI_Cart_create, we do not allow this");
        }

        /* With this assigment we loose the link to the original communicator
         * which will usually be MPI_COMM_WORLD, unless have multisim.
         */
        cr->mpi_comm_mysim = comm_cart;
        cr->sim_nodeid     = rank;

        MPI_Cart_coords(cr->mpi_comm_mysim, cr->sim_nodeid, DIM, dd->ci);

        GMX_LOG(mdlog.info).appendTextFormatted(
                "Cartesian rank %d, coordinates %d %d %d\n",
                cr->sim_nodeid, dd->ci[XX], dd->ci[YY], dd->ci[ZZ]);

        if (dd->ci[comm->cartpmedim] < dd->nc[comm->cartpmedim])
        {
            cr->duty = DUTY_PP;
        }
        if (cr->npmenodes == 0 ||
            dd->ci[comm->cartpmedim] >= dd->nc[comm->cartpmedim])
        {
            cr->duty = DUTY_PME;
        }

        /* Split the sim communicator into PP and PME only nodes */
        MPI_Comm_split(cr->mpi_comm_mysim,
                       getThisRankDuties(cr),
                       dd_index(comm->ntot, dd->ci),
                       &cr->mpi_comm_mygroup);
#endif
    }
    else
    {
        switch (rankOrder)
        {
            case DdRankOrder::pp_pme:
                GMX_LOG(mdlog.info).appendText("Order of the ranks: PP first, PME last");
                break;
            case DdRankOrder::interleave:
                /* Interleave the PP-only and PME-only ranks */
                GMX_LOG(mdlog.info).appendText("Interleaving PP and PME ranks");
                comm->pmenodes = dd_interleaved_pme_ranks(dd);
                break;
            case DdRankOrder::cartesian:
                break;
            default:
                gmx_fatal(FARGS, "Invalid ddRankOrder=%d", static_cast<int>(rankOrder));
        }

        if (dd_simnode2pmenode(dd, cr, cr->sim_nodeid) == -1)
        {
            cr->duty = DUTY_PME;
        }
        else
        {
            cr->duty = DUTY_PP;
        }
#if GMX_MPI
        /* Split the sim communicator into PP and PME only nodes */
        MPI_Comm_split(cr->mpi_comm_mysim,
                       getThisRankDuties(cr),
                       cr->nodeid,
                       &cr->mpi_comm_mygroup);
        MPI_Comm_rank(cr->mpi_comm_mygroup, &cr->nodeid);
#endif
    }

    GMX_LOG(mdlog.info).appendTextFormatted(
            "This rank does only %s work.\n",
            thisRankHasDuty(cr, DUTY_PP) ? "particle-particle" : "PME-mesh");
}

/*! \brief Generates the MPI communicators for domain decomposition */
static void make_dd_communicators(const gmx::MDLogger &mdlog,
                                  t_commrec *cr,
                                  gmx_domdec_t *dd, DdRankOrder ddRankOrder)
{
    gmx_domdec_comm_t *comm;
    bool               CartReorder;

    comm = dd->comm;

    copy_ivec(dd->nc, comm->ntot);

    comm->bCartesianPP     = (ddRankOrder == DdRankOrder::cartesian);
    comm->bCartesianPP_PME = FALSE;

    /* Reorder the nodes by default. This might change the MPI ranks.
     * Real reordering is only supported on very few architectures,
     * Blue Gene is one of them.
     */
    CartReorder = getenv("GMX_NO_CART_REORDER") == nullptr;

    if (cr->npmenodes > 0)
    {
        /* Split the communicator into a PP and PME part */
        split_communicator(mdlog, cr, dd, ddRankOrder, CartReorder);
        if (comm->bCartesianPP_PME)
        {
            /* We (possibly) reordered the nodes in split_communicator,
             * so it is no longer required in make_pp_communicator.
             */
            CartReorder = FALSE;
        }
    }
    else
    {
        /* All nodes do PP and PME */
#if GMX_MPI
        /* We do not require separate communicators */
        cr->mpi_comm_mygroup = cr->mpi_comm_mysim;
#endif
    }

    if (thisRankHasDuty(cr, DUTY_PP))
    {
        /* Copy or make a new PP communicator */
        make_pp_communicator(mdlog, dd, cr, CartReorder);
    }
    else
    {
        receive_ddindex2simnodeid(dd, cr);
    }

    if (!thisRankHasDuty(cr, DUTY_PME))
    {
        /* Set up the commnuication to our PME node */
        dd->pme_nodeid           = dd_simnode2pmenode(dd, cr, cr->sim_nodeid);
        dd->pme_receive_vir_ener = receive_vir_ener(dd, cr);
        if (debug)
        {
            fprintf(debug, "My pme_nodeid %d receive ener %s\n",
                    dd->pme_nodeid, gmx::boolToString(dd->pme_receive_vir_ener));
        }
    }
    else
    {
        dd->pme_nodeid = -1;
    }

    /* We can not use DDMASTER(dd), because dd->masterrank is set later */
    if (MASTER(cr))
    {
        dd->ma = std::make_unique<AtomDistribution>(dd->nc,
                                                    comm->cgs_gl.nr,
                                                    comm->cgs_gl.index[comm->cgs_gl.nr]);
    }
}

static real *get_slb_frac(const gmx::MDLogger &mdlog,
                          const char *dir, int nc, const char *size_string)
{
    real  *slb_frac, tot;
    int    i, n;
    double dbl;

    slb_frac = nullptr;
    if (nc > 1 && size_string != nullptr)
    {
        GMX_LOG(mdlog.info).appendTextFormatted(
                "Using static load balancing for the %s direction", dir);
        snew(slb_frac, nc);
        tot = 0;
        for (i = 0; i < nc; i++)
        {
            dbl = 0;
            sscanf(size_string, "%20lf%n", &dbl, &n);
            if (dbl == 0)
            {
                gmx_fatal(FARGS, "Incorrect or not enough DD cell size entries for direction %s: '%s'", dir, size_string);
            }
            slb_frac[i]  = dbl;
            size_string += n;
            tot         += slb_frac[i];
        }
        /* Normalize */
        std::string relativeCellSizes = "Relative cell sizes:";
        for (i = 0; i < nc; i++)
        {
            slb_frac[i]       /= tot;
            relativeCellSizes += gmx::formatString(" %5.3f", slb_frac[i]);
        }
        GMX_LOG(mdlog.info).appendText(relativeCellSizes);
    }

    return slb_frac;
}

static int multi_body_bondeds_count(const gmx_mtop_t *mtop)
{
    int                  n     = 0;
    gmx_mtop_ilistloop_t iloop = gmx_mtop_ilistloop_init(mtop);
    int                  nmol;
    while (const InteractionLists *ilists = gmx_mtop_ilistloop_next(iloop, &nmol))
    {
        for (auto &ilist : extractILists(*ilists, IF_BOND))
        {
            if (NRAL(ilist.functionType) >  2)
            {
                n += nmol*(ilist.iatoms.size()/ilistStride(ilist));
            }
        }
    }

    return n;
}

static int dd_getenv(const gmx::MDLogger &mdlog,
                     const char *env_var, int def)
{
    char *val;
    int   nst;

    nst = def;
    val = getenv(env_var);
    if (val)
    {
        if (sscanf(val, "%20d", &nst) <= 0)
        {
            nst = 1;
        }
        GMX_LOG(mdlog.info).appendTextFormatted(
                "Found env.var. %s = %s, using value %d",
                env_var, val, nst);
    }

    return nst;
}

static void check_dd_restrictions(const gmx_domdec_t  *dd,
                                  const t_inputrec    *ir,
                                  const gmx::MDLogger &mdlog)
{
    if (ir->ePBC == epbcSCREW &&
        (dd->nc[XX] == 1 || dd->nc[YY] > 1 || dd->nc[ZZ] > 1))
    {
        gmx_fatal(FARGS, "With pbc=%s can only do domain decomposition in the x-direction", epbc_names[ir->ePBC]);
    }

    if (ir->ns_type == ensSIMPLE)
    {
        gmx_fatal(FARGS, "Domain decomposition does not support simple neighbor searching, use grid searching or run with one MPI rank");
    }

    if (ir->nstlist == 0)
    {
        gmx_fatal(FARGS, "Domain decomposition does not work with nstlist=0");
    }

    if (ir->comm_mode == ecmANGULAR && ir->ePBC != epbcNONE)
    {
        GMX_LOG(mdlog.warning).appendText("comm-mode angular will give incorrect results when the comm group partially crosses a periodic boundary");
    }
}

static real average_cellsize_min(gmx_domdec_t *dd, gmx_ddbox_t *ddbox)
{
    int  di, d;
    real r;

    r = ddbox->box_size[XX];
    for (di = 0; di < dd->ndim; di++)
    {
        d = dd->dim[di];
        /* Check using the initial average cell size */
        r = std::min(r, ddbox->box_size[d]*ddbox->skew_fac[d]/dd->nc[d]);
    }

    return r;
}

/*! \brief Depending on the DLB initial value return the DLB switched off state or issue an error.
 */
static DlbState forceDlbOffOrBail(DlbState             cmdlineDlbState,
                                  const std::string   &reasonStr,
                                  const gmx::MDLogger &mdlog)
{
    std::string dlbNotSupportedErr  = "Dynamic load balancing requested, but ";
    std::string dlbDisableNote      = "NOTE: disabling dynamic load balancing as ";

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
 * \param [in] mdlog       Logger.
 * \param [in] dlbOption   Enum value for the DLB option.
 * \param [in] bRecordLoad True if the load balancer is recording load information.
 * \param [in] mdrunOptions  Options for mdrun.
 * \param [in] ir          Pointer mdrun to input parameters.
 * \returns                DLB initial/startup state.
 */
static DlbState determineInitialDlbState(const gmx::MDLogger &mdlog,
                                         DlbOption dlbOption, gmx_bool bRecordLoad,
                                         const MdrunOptions &mdrunOptions,
                                         const t_inputrec *ir)
{
    DlbState dlbState = DlbState::offCanTurnOn;

    switch (dlbOption)
    {
        case DlbOption::turnOnWhenUseful: dlbState = DlbState::offCanTurnOn; break;
        case DlbOption::no:               dlbState = DlbState::offUser;      break;
        case DlbOption::yes:              dlbState = DlbState::onUser;       break;
        default: gmx_incons("Invalid dlbOption enum value");
    }

    /* Reruns don't support DLB: bail or override auto mode */
    if (mdrunOptions.rerun)
    {
        std::string reasonStr = "it is not supported in reruns.";
        return forceDlbOffOrBail(dlbState, reasonStr, mdlog);
    }

    /* Unsupported integrators */
    if (!EI_DYNAMICS(ir->eI))
    {
        auto reasonStr = gmx::formatString("it is only supported with dynamics, not with integrator '%s'.", EI(ir->eI));
        return forceDlbOffOrBail(dlbState, reasonStr, mdlog);
    }

    /* Without cycle counters we can't time work to balance on */
    if (!bRecordLoad)
    {
        std::string reasonStr = "cycle counters unsupported or not enabled in the operating system kernel.";
        return forceDlbOffOrBail(dlbState, reasonStr, mdlog);
    }

    if (mdrunOptions.reproducible)
    {
        std::string reasonStr = "you started a reproducible run.";
        switch (dlbState)
        {
            case DlbState::offUser:
                break;
            case DlbState::offForever:
                GMX_RELEASE_ASSERT(false, "DlbState::offForever is not a valid initial state");
                break;
            case DlbState::offCanTurnOn:
                return forceDlbOffOrBail(dlbState, reasonStr, mdlog);
            case DlbState::onCanTurnOff:
                GMX_RELEASE_ASSERT(false, "DlbState::offCanTurnOff is not a valid initial state");
                break;
            case DlbState::onUser:
                return forceDlbOffOrBail(dlbState, reasonStr + " In load balanced runs binary reproducibility cannot be ensured.", mdlog);
            default:
                gmx_fatal(FARGS, "Death horror: undefined case (%d) for load balancing choice", static_cast<int>(dlbState));
        }
    }

    return dlbState;
}

static void set_dd_dim(const gmx::MDLogger &mdlog, gmx_domdec_t *dd)
{
    dd->ndim = 0;
    if (getenv("GMX_DD_ORDER_ZYX") != nullptr)
    {
        /* Decomposition order z,y,x */
        GMX_LOG(mdlog.info).appendText("Using domain decomposition order z, y, x");
        for (int dim = DIM-1; dim >= 0; dim--)
        {
            if (dd->nc[dim] > 1)
            {
                dd->dim[dd->ndim++] = dim;
            }
        }
    }
    else
    {
        /* Decomposition order x,y,z */
        for (int dim = 0; dim < DIM; dim++)
        {
            if (dd->nc[dim] > 1)
            {
                dd->dim[dd->ndim++] = dim;
            }
        }
    }

    if (dd->ndim == 0)
    {
        /* Set dim[0] to avoid extra checks on ndim in several places */
        dd->dim[0] = XX;
    }
}

static gmx_domdec_comm_t *init_dd_comm()
{
    gmx_domdec_comm_t *comm = new gmx_domdec_comm_t;

    comm->n_load_have      = 0;
    comm->n_load_collect   = 0;

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
    comm->load_mdf  = 0;
    comm->load_pme  = 0;

    /* This should be replaced by a unique pointer */
    comm->balanceRegion = ddBalanceRegionAllocate();

    return comm;
}

/* Returns whether mtop contains constraints and/or vsites */
static bool systemHasConstraintsOrVsites(const gmx_mtop_t &mtop)
{
    auto ilistLoop = gmx_mtop_ilistloop_init(mtop);
    int  nmol;
    while (const InteractionLists *ilists = gmx_mtop_ilistloop_next(ilistLoop, &nmol))
    {
        if (!extractILists(*ilists, IF_CONSTRAINT | IF_VSITE).empty())
        {
            return true;
        }
    }

    return false;
}

static void setupUpdateGroups(const gmx::MDLogger &mdlog,
                              const gmx_mtop_t    &mtop,
                              const t_inputrec    &inputrec,
                              real                 cutoffMargin,
                              int                  numMpiRanksTotal,
                              gmx_domdec_comm_t   *comm)
{
    /* When we have constraints and/or vsites, it is beneficial to use
     * update groups (when possible) to allow independent update of groups.
     */
    if (!systemHasConstraintsOrVsites(mtop))
    {
        /* No constraints or vsites, atoms can be updated independently */
        return;
    }

    comm->updateGroupingPerMoleculetype = gmx::makeUpdateGroups(mtop);
    comm->useUpdateGroups               =
        (!comm->updateGroupingPerMoleculetype.empty() &&
         getenv("GMX_NO_UPDATEGROUPS") == nullptr);

    if (comm->useUpdateGroups)
    {
        int numUpdateGroups = 0;
        for (const auto &molblock : mtop.molblock)
        {
            numUpdateGroups += molblock.nmol*comm->updateGroupingPerMoleculetype[molblock.type].numBlocks();
        }

        /* Note: We would like to use dd->nnodes for the atom count estimate,
         *       but that is not yet available here. But this anyhow only
         *       affect performance up to the second dd_partition_system call.
         */
        int homeAtomCountEstimate =  mtop.natoms/numMpiRanksTotal;
        comm->updateGroupsCog =
            std::make_unique<gmx::UpdateGroupsCog>(mtop,
                                                   comm->updateGroupingPerMoleculetype,
                                                   maxReferenceTemperature(inputrec),
                                                   homeAtomCountEstimate);

        /* To use update groups, the large domain-to-domain cutoff distance
         * should be compatible with the box size.
         */
        comm->useUpdateGroups = (atomToAtomIntoDomainToDomainCutoff(*comm, 0) < cutoffMargin);

        if (comm->useUpdateGroups)
        {
            GMX_LOG(mdlog.info).appendTextFormatted(
                    "Using update groups, nr %d, average size %.1f atoms, max. radius %.3f nm\n",
                    numUpdateGroups,
                    mtop.natoms/static_cast<double>(numUpdateGroups),
                    comm->updateGroupsCog->maxUpdateGroupRadius());
        }
        else
        {
            GMX_LOG(mdlog.info).appendTextFormatted("The combination of rlist and box size prohibits the use of update groups\n");
            comm->updateGroupingPerMoleculetype.clear();
            comm->updateGroupsCog.reset(nullptr);
        }
    }
}

/*! \brief Set the cell size and interaction limits, as well as the DD grid */
static void set_dd_limits_and_grid(const gmx::MDLogger &mdlog,
                                   t_commrec *cr, gmx_domdec_t *dd,
                                   const DomdecOptions &options,
                                   const MdrunOptions &mdrunOptions,
                                   const gmx_mtop_t *mtop,
                                   const t_inputrec *ir,
                                   const matrix box,
                                   gmx::ArrayRef<const gmx::RVec> xGlobal,
                                   gmx_ddbox_t *ddbox)
{
    real               r_bonded         = -1;
    real               r_bonded_limit   = -1;
    const real         tenPercentMargin = 1.1;
    gmx_domdec_comm_t *comm             = dd->comm;

    dd->npbcdim              = ePBC2npbcdim(ir->ePBC);
    dd->numBoundedDimensions = inputrec2nboundeddim(ir);
    dd->haveDynamicBox       = inputrecDynamicBox(ir);
    dd->bScrewPBC            = (ir->ePBC == epbcSCREW);

    dd->pme_recv_f_alloc = 0;
    dd->pme_recv_f_buf   = nullptr;

    /* Initialize to GPU share count to 0, might change later */
    comm->nrank_gpu_shared = 0;

    comm->dlbState         = determineInitialDlbState(mdlog, options.dlbOption, comm->bRecordLoad, mdrunOptions, ir);
    dd_dlb_set_should_check_whether_to_turn_dlb_on(dd, TRUE);
    /* To consider turning DLB on after 2*nstlist steps we need to check
     * at partitioning count 3. Thus we need to increase the first count by 2.
     */
    comm->ddPartioningCountFirstDlbOff += 2;

    GMX_LOG(mdlog.info).appendTextFormatted(
            "Dynamic load balancing: %s", edlbs_names[int(comm->dlbState)]);

    comm->bPMELoadBalDLBLimits = FALSE;

    /* Allocate the charge group/atom sorting struct */
    comm->sort = std::make_unique<gmx_domdec_sort_t>();

    comm->bCGs = (ncg_mtop(mtop) < mtop->natoms);

    /* We need to decide on update groups early, as this affects communication distances */
    comm->useUpdateGroups = false;
    if (ir->cutoff_scheme == ecutsVERLET)
    {
        real cutoffMargin = std::sqrt(max_cutoff2(ir->ePBC, box)) - ir->rlist;
        setupUpdateGroups(mdlog, *mtop, *ir, cutoffMargin, cr->nnodes, comm);
    }

    comm->bInterCGBondeds = ((ncg_mtop(mtop) > gmx_mtop_num_molecules(*mtop)) ||
                             mtop->bIntermolecularInteractions);
    if (comm->bInterCGBondeds)
    {
        comm->bInterCGMultiBody = (multi_body_bondeds_count(mtop) > 0);
    }
    else
    {
        comm->bInterCGMultiBody = FALSE;
    }

    if (comm->useUpdateGroups)
    {
        dd->splitConstraints = false;
        dd->splitSettles     = false;
    }
    else
    {
        dd->splitConstraints = gmx::inter_charge_group_constraints(*mtop);
        dd->splitSettles     = gmx::inter_charge_group_settles(*mtop);
    }

    if (ir->rlist == 0)
    {
        /* Set the cut-off to some very large value,
         * so we don't need if statements everywhere in the code.
         * We use sqrt, since the cut-off is squared in some places.
         */
        comm->cutoff   = GMX_CUTOFF_INF;
    }
    else
    {
        comm->cutoff   = atomToAtomIntoDomainToDomainCutoff(*comm, ir->rlist);
    }
    comm->cutoff_mbody = 0;

    /* Determine the minimum cell size limit, affected by many factors */
    comm->cellsize_limit = 0;
    comm->bBondComm      = FALSE;

    /* We do not allow home atoms to move beyond the neighboring domain
     * between domain decomposition steps, which limits the cell size.
     * Get an estimate of cell size limit due to atom displacement.
     * In most cases this is a large overestimate, because it assumes
     * non-interaction atoms.
     * We set the chance to 1 in a trillion steps.
     */
    constexpr real c_chanceThatAtomMovesBeyondDomain = 1e-12;
    const real     limitForAtomDisplacement          =
        minCellSizeForAtomDisplacement(*mtop, *ir,
                                       comm->updateGroupingPerMoleculetype,
                                       c_chanceThatAtomMovesBeyondDomain);
    GMX_LOG(mdlog.info).appendTextFormatted(
            "Minimum cell size due to atom displacement: %.3f nm",
            limitForAtomDisplacement);

    comm->cellsize_limit = std::max(comm->cellsize_limit,
                                    limitForAtomDisplacement);

    /* TODO: PME decomposition currently requires atoms not to be more than
     *       2/3 of comm->cutoff, which is >=rlist, outside of their domain.
     *       In nearly all cases, limitForAtomDisplacement will be smaller
     *       than 2/3*rlist, so the PME requirement is satisfied.
     *       But it would be better for both correctness and performance
     *       to use limitForAtomDisplacement instead of 2/3*comm->cutoff.
     *       Note that we would need to improve the pairlist buffer case.
     */

    if (comm->bInterCGBondeds)
    {
        if (options.minimumCommunicationRange > 0)
        {
            comm->cutoff_mbody = atomToAtomIntoDomainToDomainCutoff(*comm, options.minimumCommunicationRange);
            if (options.useBondedCommunication)
            {
                comm->bBondComm = (comm->cutoff_mbody > comm->cutoff);
            }
            else
            {
                comm->cutoff = std::max(comm->cutoff, comm->cutoff_mbody);
            }
            r_bonded_limit = comm->cutoff_mbody;
        }
        else if (ir->bPeriodicMols)
        {
            /* Can not easily determine the required cut-off */
            GMX_LOG(mdlog.warning).appendText("NOTE: Periodic molecules are present in this system. Because of this, the domain decomposition algorithm cannot easily determine the minimum cell size that it requires for treating bonded interactions. Instead, domain decomposition will assume that half the non-bonded cut-off will be a suitable lower bound.");
            comm->cutoff_mbody = comm->cutoff/2;
            r_bonded_limit     = comm->cutoff_mbody;
        }
        else
        {
            real r_2b, r_mb;

            if (MASTER(cr))
            {
                dd_bonded_cg_distance(mdlog, mtop, ir, as_rvec_array(xGlobal.data()), box,
                                      options.checkBondedInteractions,
                                      &r_2b, &r_mb);
            }
            gmx_bcast(sizeof(r_2b), &r_2b, cr);
            gmx_bcast(sizeof(r_mb), &r_mb, cr);

            /* We use an initial margin of 10% for the minimum cell size,
             * except when we are just below the non-bonded cut-off.
             */
            if (options.useBondedCommunication)
            {
                if (std::max(r_2b, r_mb) > comm->cutoff)
                {
                    r_bonded        = std::max(r_2b, r_mb);
                    r_bonded_limit  = tenPercentMargin*r_bonded;
                    comm->bBondComm = TRUE;
                }
                else
                {
                    r_bonded       = r_mb;
                    r_bonded_limit = std::min(tenPercentMargin*r_bonded, comm->cutoff);
                }
                /* We determine cutoff_mbody later */
            }
            else
            {
                /* No special bonded communication,
                 * simply increase the DD cut-off.
                 */
                r_bonded_limit     = tenPercentMargin*std::max(r_2b, r_mb);
                comm->cutoff_mbody = r_bonded_limit;
                comm->cutoff       = std::max(comm->cutoff, comm->cutoff_mbody);
            }
        }
        GMX_LOG(mdlog.info).appendTextFormatted(
                "Minimum cell size due to bonded interactions: %.3f nm",
                r_bonded_limit);

        comm->cellsize_limit = std::max(comm->cellsize_limit, r_bonded_limit);
    }

    real rconstr = 0;
    if (dd->splitConstraints && options.constraintCommunicationRange <= 0)
    {
        /* There is a cell size limit due to the constraints (P-LINCS) */
        rconstr = gmx::constr_r_max(mdlog, mtop, ir);
        GMX_LOG(mdlog.info).appendTextFormatted(
                "Estimated maximum distance required for P-LINCS: %.3f nm",
                rconstr);
        if (rconstr > comm->cellsize_limit)
        {
            GMX_LOG(mdlog.info).appendText("This distance will limit the DD cell size, you can override this with -rcon");
        }
    }
    else if (options.constraintCommunicationRange > 0)
    {
        /* Here we do not check for dd->splitConstraints.
         * because one can also set a cell size limit for virtual sites only
         * and at this point we don't know yet if there are intercg v-sites.
         */
        GMX_LOG(mdlog.info).appendTextFormatted(
                "User supplied maximum distance required for P-LINCS: %.3f nm",
                options.constraintCommunicationRange);
        rconstr = options.constraintCommunicationRange;
    }
    comm->cellsize_limit = std::max(comm->cellsize_limit, rconstr);

    comm->cgs_gl = gmx_mtop_global_cgs(mtop);

    if (options.numCells[XX] > 0)
    {
        copy_ivec(options.numCells, dd->nc);
        set_dd_dim(mdlog, dd);
        set_ddbox_cr(*cr, &dd->nc, *ir, box, xGlobal, ddbox);

        if (options.numPmeRanks >= 0)
        {
            cr->npmenodes = options.numPmeRanks;
        }
        else
        {
            /* When the DD grid is set explicitly and -npme is set to auto,
             * don't use PME ranks. We check later if the DD grid is
             * compatible with the total number of ranks.
             */
            cr->npmenodes = 0;
        }

        real acs = average_cellsize_min(dd, ddbox);
        if (acs < comm->cellsize_limit)
        {
            gmx_fatal_collective(FARGS, cr->mpi_comm_mysim, MASTER(cr),
                                 "The initial cell size (%f) is smaller than the cell size limit (%f), change options -dd, -rdd or -rcon, see the log file for details",
                                 acs, comm->cellsize_limit);
        }
    }
    else
    {
        set_ddbox_cr(*cr, nullptr, *ir, box, xGlobal, ddbox);

        /* We need to choose the optimal DD grid and possibly PME nodes */
        real limit =
            dd_choose_grid(mdlog, cr, dd, ir, mtop, box, ddbox,
                           options.numPmeRanks,
                           !isDlbDisabled(comm),
                           options.dlbScaling,
                           comm->cellsize_limit, comm->cutoff,
                           comm->bInterCGBondeds);

        if (dd->nc[XX] == 0)
        {
            char     buf[STRLEN];
            gmx_bool bC = (dd->splitConstraints && rconstr > r_bonded_limit);
            sprintf(buf, "Change the number of ranks or mdrun option %s%s%s",
                    !bC ? "-rdd" : "-rcon",
                    comm->dlbState != DlbState::offUser ? " or -dds" : "",
                    bC ? " or your LINCS settings" : "");

            gmx_fatal_collective(FARGS, cr->mpi_comm_mysim, MASTER(cr),
                                 "There is no domain decomposition for %d ranks that is compatible with the given box and a minimum cell size of %g nm\n"
                                 "%s\n"
                                 "Look in the log file for details on the domain decomposition",
                                 cr->nnodes-cr->npmenodes, limit, buf);
        }
        set_dd_dim(mdlog, dd);
    }

    GMX_LOG(mdlog.info).appendTextFormatted(
            "Domain decomposition grid %d x %d x %d, separate PME ranks %d",
            dd->nc[XX], dd->nc[YY], dd->nc[ZZ], cr->npmenodes);

    dd->nnodes = dd->nc[XX]*dd->nc[YY]*dd->nc[ZZ];
    if (cr->nnodes - dd->nnodes != cr->npmenodes)
    {
        gmx_fatal_collective(FARGS, cr->mpi_comm_mysim, MASTER(cr),
                             "The size of the domain decomposition grid (%d) does not match the number of ranks (%d). The total number of ranks is %d",
                             dd->nnodes, cr->nnodes - cr->npmenodes, cr->nnodes);
    }
    if (cr->npmenodes > dd->nnodes)
    {
        gmx_fatal_collective(FARGS, cr->mpi_comm_mysim, MASTER(cr),
                             "The number of separate PME ranks (%d) is larger than the number of PP ranks (%d), this is not supported.", cr->npmenodes, dd->nnodes);
    }
    if (cr->npmenodes > 0)
    {
        comm->npmenodes = cr->npmenodes;
    }
    else
    {
        comm->npmenodes = dd->nnodes;
    }

    if (EEL_PME(ir->coulombtype) || EVDW_PME(ir->vdwtype))
    {
        /* The following choices should match those
         * in comm_cost_est in domdec_setup.c.
         * Note that here the checks have to take into account
         * that the decomposition might occur in a different order than xyz
         * (for instance through the env.var. GMX_DD_ORDER_ZYX),
         * in which case they will not match those in comm_cost_est,
         * but since that is mainly for testing purposes that's fine.
         */
        if (dd->ndim >= 2 && dd->dim[0] == XX && dd->dim[1] == YY &&
            comm->npmenodes > dd->nc[XX] && comm->npmenodes % dd->nc[XX] == 0 &&
            getenv("GMX_PMEONEDD") == nullptr)
        {
            comm->npmedecompdim = 2;
            comm->npmenodes_x   = dd->nc[XX];
            comm->npmenodes_y   = comm->npmenodes/comm->npmenodes_x;
        }
        else
        {
            /* In case nc is 1 in both x and y we could still choose to
             * decompose pme in y instead of x, but we use x for simplicity.
             */
            comm->npmedecompdim = 1;
            if (dd->dim[0] == YY)
            {
                comm->npmenodes_x = 1;
                comm->npmenodes_y = comm->npmenodes;
            }
            else
            {
                comm->npmenodes_x = comm->npmenodes;
                comm->npmenodes_y = 1;
            }
        }
        GMX_LOG(mdlog.info).appendTextFormatted(
                "PME domain decomposition: %d x %d x %d",
                comm->npmenodes_x, comm->npmenodes_y, 1);
    }
    else
    {
        comm->npmedecompdim = 0;
        comm->npmenodes_x   = 0;
        comm->npmenodes_y   = 0;
    }

    snew(comm->slb_frac, DIM);
    if (isDlbDisabled(comm))
    {
        comm->slb_frac[XX] = get_slb_frac(mdlog, "x", dd->nc[XX], options.cellSizeX);
        comm->slb_frac[YY] = get_slb_frac(mdlog, "y", dd->nc[YY], options.cellSizeY);
        comm->slb_frac[ZZ] = get_slb_frac(mdlog, "z", dd->nc[ZZ], options.cellSizeZ);
    }

    if (comm->bInterCGBondeds && comm->cutoff_mbody == 0)
    {
        if (comm->bBondComm || !isDlbDisabled(comm))
        {
            /* Set the bonded communication distance to halfway
             * the minimum and the maximum,
             * since the extra communication cost is nearly zero.
             */
            real acs           = average_cellsize_min(dd, ddbox);
            comm->cutoff_mbody = 0.5*(r_bonded + acs);
            if (!isDlbDisabled(comm))
            {
                /* Check if this does not limit the scaling */
                comm->cutoff_mbody = std::min(comm->cutoff_mbody,
                                              options.dlbScaling*acs);
            }
            if (!comm->bBondComm)
            {
                /* Without bBondComm do not go beyond the n.b. cut-off */
                comm->cutoff_mbody = std::min(comm->cutoff_mbody, comm->cutoff);
                if (comm->cellsize_limit >= comm->cutoff)
                {
                    /* We don't loose a lot of efficieny
                     * when increasing it to the n.b. cut-off.
                     * It can even be slightly faster, because we need
                     * less checks for the communication setup.
                     */
                    comm->cutoff_mbody = comm->cutoff;
                }
            }
            /* Check if we did not end up below our original limit */
            comm->cutoff_mbody = std::max(comm->cutoff_mbody, r_bonded_limit);

            if (comm->cutoff_mbody > comm->cellsize_limit)
            {
                comm->cellsize_limit = comm->cutoff_mbody;
            }
        }
        /* Without DLB and cutoff_mbody<cutoff, cutoff_mbody is dynamic */
    }

    if (debug)
    {
        fprintf(debug, "Bonded atom communication beyond the cut-off: %s\n"
                "cellsize limit %f\n",
                gmx::boolToString(comm->bBondComm), comm->cellsize_limit);
    }

    if (MASTER(cr))
    {
        check_dd_restrictions(dd, ir, mdlog);
    }
}

static char *init_bLocalCG(const gmx_mtop_t *mtop)
{
    int   ncg, cg;
    char *bLocalCG;

    ncg = ncg_mtop(mtop);
    snew(bLocalCG, ncg);
    for (cg = 0; cg < ncg; cg++)
    {
        bLocalCG[cg] = FALSE;
    }

    return bLocalCG;
}

void dd_init_bondeds(FILE *fplog,
                     gmx_domdec_t *dd,
                     const gmx_mtop_t *mtop,
                     const gmx_vsite_t *vsite,
                     const t_inputrec *ir,
                     gmx_bool bBCheck, cginfo_mb_t *cginfo_mb)
{
    gmx_domdec_comm_t *comm;

    dd_make_reverse_top(fplog, dd, mtop, vsite, ir, bBCheck);

    comm = dd->comm;

    if (comm->bBondComm)
    {
        /* Communicate atoms beyond the cut-off for bonded interactions */
        comm = dd->comm;

        comm->cglink = make_charge_group_links(mtop, dd, cginfo_mb);

        comm->bLocalCG = init_bLocalCG(mtop);
    }
    else
    {
        /* Only communicate atoms based on cut-off */
        comm->cglink   = nullptr;
        comm->bLocalCG = nullptr;
    }
}

static void writeSettings(gmx::TextWriter       *log,
                          gmx_domdec_t          *dd,
                          const gmx_mtop_t      *mtop,
                          const t_inputrec      *ir,
                          gmx_bool               bDynLoadBal,
                          real                   dlb_scale,
                          const gmx_ddbox_t     *ddbox)
{
    gmx_domdec_comm_t *comm;
    int                d;
    ivec               np;
    real               limit, shrink;

    comm = dd->comm;

    if (bDynLoadBal)
    {
        log->writeString("The maximum number of communication pulses is:");
        for (d = 0; d < dd->ndim; d++)
        {
            log->writeStringFormatted(" %c %d", dim2char(dd->dim[d]), comm->cd[d].np_dlb);
        }
        log->ensureLineBreak();
        log->writeLineFormatted("The minimum size for domain decomposition cells is %.3f nm", comm->cellsize_limit);
        log->writeLineFormatted("The requested allowed shrink of DD cells (option -dds) is: %.2f", dlb_scale);
        log->writeString("The allowed shrink of domain decomposition cells is:");
        for (d = 0; d < DIM; d++)
        {
            if (dd->nc[d] > 1)
            {
                if (d >= ddbox->npbcdim && dd->nc[d] == 2)
                {
                    shrink = 0;
                }
                else
                {
                    shrink =
                        comm->cellsize_min_dlb[d]/
                        (ddbox->box_size[d]*ddbox->skew_fac[d]/dd->nc[d]);
                }
                log->writeStringFormatted(" %c %.2f", dim2char(d), shrink);
            }
        }
        log->ensureLineBreak();
    }
    else
    {
        set_dd_cell_sizes_slb(dd, ddbox, setcellsizeslbPULSE_ONLY, np);
        log->writeString("The initial number of communication pulses is:");
        for (d = 0; d < dd->ndim; d++)
        {
            log->writeStringFormatted(" %c %d", dim2char(dd->dim[d]), np[dd->dim[d]]);
        }
        log->ensureLineBreak();
        log->writeString("The initial domain decomposition cell size is:");
        for (d = 0; d < DIM; d++)
        {
            if (dd->nc[d] > 1)
            {
                log->writeStringFormatted(" %c %.2f nm",
                                          dim2char(d), dd->comm->cellsize_min[d]);
            }
        }
        log->ensureLineBreak();
        log->writeLine();
    }

    gmx_bool bInterCGVsites = count_intercg_vsites(mtop) != 0;

    if (comm->bInterCGBondeds ||
        bInterCGVsites ||
        dd->splitConstraints || dd->splitSettles)
    {
        std::string decompUnits;
        if (comm->bCGs)
        {
            decompUnits = "charge groups";
        }
        else if (comm->useUpdateGroups)
        {
            decompUnits = "atom groups";
        }
        else
        {
            decompUnits = "atoms";
        }

        log->writeLineFormatted("The maximum allowed distance for %s involved in interactions is:", decompUnits.c_str());
        log->writeLineFormatted("%40s  %-7s %6.3f nm", "non-bonded interactions", "", comm->cutoff);

        if (bDynLoadBal)
        {
            limit = dd->comm->cellsize_limit;
        }
        else
        {
            if (dynamic_dd_box(*dd))
            {
                log->writeLine("(the following are initial values, they could change due to box deformation)");
            }
            limit = dd->comm->cellsize_min[XX];
            for (d = 1; d < DIM; d++)
            {
                limit = std::min(limit, dd->comm->cellsize_min[d]);
            }
        }

        if (comm->bInterCGBondeds)
        {
            log->writeLineFormatted("%40s  %-7s %6.3f nm",
                                    "two-body bonded interactions", "(-rdd)",
                                    std::max(comm->cutoff, comm->cutoff_mbody));
            log->writeLineFormatted("%40s  %-7s %6.3f nm",
                                    "multi-body bonded interactions", "(-rdd)",
                                    (comm->bBondComm || isDlbOn(dd->comm)) ? comm->cutoff_mbody : std::min(comm->cutoff, limit));
        }
        if (bInterCGVsites)
        {
            log->writeLineFormatted("%40s  %-7s %6.3f nm",
                                    "virtual site constructions", "(-rcon)", limit);
        }
        if (dd->splitConstraints || dd->splitSettles)
        {
            std::string separation = gmx::formatString("atoms separated by up to %d constraints",
                                                       1+ir->nProjOrder);
            log->writeLineFormatted("%40s  %-7s %6.3f nm\n",
                                    separation.c_str(), "(-rcon)", limit);
        }
        log->ensureLineBreak();
    }
}

static void logSettings(const gmx::MDLogger &mdlog,
                        gmx_domdec_t        *dd,
                        const gmx_mtop_t    *mtop,
                        const t_inputrec    *ir,
                        real                 dlb_scale,
                        const gmx_ddbox_t   *ddbox)
{
    gmx::StringOutputStream stream;
    gmx::TextWriter         log(&stream);
    writeSettings(&log, dd, mtop, ir, isDlbOn(dd->comm), dlb_scale, ddbox);
    if (dd->comm->dlbState == DlbState::offCanTurnOn)
    {
        {
            log.ensureEmptyLine();
            log.writeLine("When dynamic load balancing gets turned on, these settings will change to:");
        }
        writeSettings(&log, dd, mtop, ir, true, dlb_scale, ddbox);
    }
    GMX_LOG(mdlog.info).asParagraph().appendText(stream.toString());
}

static void set_cell_limits_dlb(const gmx::MDLogger &mdlog,
                                gmx_domdec_t        *dd,
                                real                 dlb_scale,
                                const t_inputrec    *ir,
                                const gmx_ddbox_t   *ddbox)
{
    gmx_domdec_comm_t *comm;
    int                d, dim, npulse, npulse_d_max, npulse_d;
    gmx_bool           bNoCutOff;

    comm = dd->comm;

    bNoCutOff = (ir->rvdw == 0 || ir->rcoulomb == 0);

    /* Determine the maximum number of comm. pulses in one dimension */

    comm->cellsize_limit = std::max(comm->cellsize_limit, comm->cutoff_mbody);

    /* Determine the maximum required number of grid pulses */
    if (comm->cellsize_limit >= comm->cutoff)
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
        npulse = static_cast<int>(0.96 + comm->cutoff/comm->cellsize_limit);
    }
    else
    {
        /* There is no cell size limit */
        npulse = std::max(dd->nc[XX]-1, std::max(dd->nc[YY]-1, dd->nc[ZZ]-1));
    }

    if (!bNoCutOff && npulse > 1)
    {
        /* See if we can do with less pulses, based on dlb_scale */
        npulse_d_max = 0;
        for (d = 0; d < dd->ndim; d++)
        {
            dim      = dd->dim[d];
            npulse_d = static_cast<int>(1 + dd->nc[dim]*comm->cutoff
                                        /(ddbox->box_size[dim]*ddbox->skew_fac[dim]*dlb_scale));
            npulse_d_max = std::max(npulse_d_max, npulse_d);
        }
        npulse = std::min(npulse, npulse_d_max);
    }

    /* This env var can override npulse */
    d = dd_getenv(mdlog, "GMX_DD_NPULSE", 0);
    if (d > 0)
    {
        npulse = d;
    }

    comm->maxpulse       = 1;
    comm->bVacDLBNoLimit = (ir->ePBC == epbcNONE);
    for (d = 0; d < dd->ndim; d++)
    {
        comm->cd[d].np_dlb    = std::min(npulse, dd->nc[dd->dim[d]]-1);
        comm->maxpulse        = std::max(comm->maxpulse, comm->cd[d].np_dlb);
        if (comm->cd[d].np_dlb < dd->nc[dd->dim[d]]-1)
        {
            comm->bVacDLBNoLimit = FALSE;
        }
    }

    /* cellsize_limit is set for LINCS in init_domain_decomposition */
    if (!comm->bVacDLBNoLimit)
    {
        comm->cellsize_limit = std::max(comm->cellsize_limit,
                                        comm->cutoff/comm->maxpulse);
    }
    comm->cellsize_limit = std::max(comm->cellsize_limit, comm->cutoff_mbody);
    /* Set the minimum cell size for each DD dimension */
    for (d = 0; d < dd->ndim; d++)
    {
        if (comm->bVacDLBNoLimit ||
            comm->cd[d].np_dlb*comm->cellsize_limit >= comm->cutoff)
        {
            comm->cellsize_min_dlb[dd->dim[d]] = comm->cellsize_limit;
        }
        else
        {
            comm->cellsize_min_dlb[dd->dim[d]] =
                comm->cutoff/comm->cd[d].np_dlb;
        }
    }
    if (comm->cutoff_mbody <= 0)
    {
        comm->cutoff_mbody = std::min(comm->cutoff, comm->cellsize_limit);
    }
    if (isDlbOn(comm))
    {
        set_dlb_limits(dd);
    }
}

gmx_bool dd_bonded_molpbc(const gmx_domdec_t *dd, int ePBC)
{
    /* If each molecule is a single charge group
     * or we use domain decomposition for each periodic dimension,
     * we do not need to take pbc into account for the bonded interactions.
     */
    return (ePBC != epbcNONE && dd->comm->bInterCGBondeds &&
            !(dd->nc[XX] > 1 &&
              dd->nc[YY] > 1 &&
              (dd->nc[ZZ] > 1 || ePBC == epbcXY)));
}

/*! \brief Sets grid size limits and PP-PME setup, prints settings to log */
static void set_ddgrid_parameters(const gmx::MDLogger &mdlog,
                                  gmx_domdec_t *dd, real dlb_scale,
                                  const gmx_mtop_t *mtop, const t_inputrec *ir,
                                  const gmx_ddbox_t *ddbox)
{
    gmx_domdec_comm_t *comm;
    int                natoms_tot;
    real               vol_frac;

    comm = dd->comm;

    if (EEL_PME(ir->coulombtype) || EVDW_PME(ir->vdwtype))
    {
        init_ddpme(dd, &comm->ddpme[0], 0);
        if (comm->npmedecompdim >= 2)
        {
            init_ddpme(dd, &comm->ddpme[1], 1);
        }
    }
    else
    {
        comm->npmenodes = 0;
        if (dd->pme_nodeid >= 0)
        {
            gmx_fatal_collective(FARGS, dd->mpi_comm_all, DDMASTER(dd),
                                 "Can not have separate PME ranks without PME electrostatics");
        }
    }

    if (debug)
    {
        fprintf(debug, "The DD cut-off is %f\n", comm->cutoff);
    }
    if (!isDlbDisabled(comm))
    {
        set_cell_limits_dlb(mdlog, dd, dlb_scale, ir, ddbox);
    }

    logSettings(mdlog, dd, mtop, ir, dlb_scale, ddbox);

    if (ir->ePBC == epbcNONE)
    {
        vol_frac = 1 - 1/static_cast<double>(dd->nnodes);
    }
    else
    {
        vol_frac =
            (1 + comm_box_frac(dd->nc, comm->cutoff, ddbox))/static_cast<double>(dd->nnodes);
    }
    if (debug)
    {
        fprintf(debug, "Volume fraction for all DD zones: %f\n", vol_frac);
    }
    natoms_tot = comm->cgs_gl.index[comm->cgs_gl.nr];

    dd->ga2la  = new gmx_ga2la_t(natoms_tot,
                                 static_cast<int>(vol_frac*natoms_tot));
}

/*! \brief Set some important DD parameters that can be modified by env.vars */
static void set_dd_envvar_options(const gmx::MDLogger &mdlog,
                                  gmx_domdec_t *dd, int rank_mysim)
{
    gmx_domdec_comm_t *comm = dd->comm;

    dd->bSendRecv2      = (dd_getenv(mdlog, "GMX_DD_USE_SENDRECV2", 0) != 0);
    comm->dlb_scale_lim = dd_getenv(mdlog, "GMX_DLB_MAX_BOX_SCALING", 10);
    comm->eFlop         = dd_getenv(mdlog, "GMX_DLB_BASED_ON_FLOPS", 0);
    int recload         = dd_getenv(mdlog, "GMX_DD_RECORD_LOAD", 1);
    comm->nstDDDump     = dd_getenv(mdlog, "GMX_DD_NST_DUMP", 0);
    comm->nstDDDumpGrid = dd_getenv(mdlog, "GMX_DD_NST_DUMP_GRID", 0);
    comm->DD_debug      = dd_getenv(mdlog, "GMX_DD_DEBUG", 0);

    if (dd->bSendRecv2)
    {
        GMX_LOG(mdlog.info).appendText("Will use two sequential MPI_Sendrecv calls instead of two simultaneous non-blocking MPI_Irecv and MPI_Isend pairs for constraint and vsite communication");
    }

    if (comm->eFlop)
    {
        GMX_LOG(mdlog.info).appendText("Will load balance based on FLOP count");
        if (comm->eFlop > 1)
        {
            srand(1 + rank_mysim);
        }
        comm->bRecordLoad = TRUE;
    }
    else
    {
        comm->bRecordLoad = (wallcycle_have_counter() && recload > 0);
    }
}

gmx_domdec_t *init_domain_decomposition(const gmx::MDLogger           &mdlog,
                                        t_commrec                     *cr,
                                        const DomdecOptions           &options,
                                        const MdrunOptions            &mdrunOptions,
                                        const gmx_mtop_t              *mtop,
                                        const t_inputrec              *ir,
                                        const matrix                   box,
                                        gmx::ArrayRef<const gmx::RVec> xGlobal,
                                        gmx::LocalAtomSetManager      *atomSets)
{
    gmx_domdec_t      *dd;

    GMX_LOG(mdlog.info).appendTextFormatted(
            "\nInitializing Domain Decomposition on %d ranks", cr->nnodes);

    dd = new gmx_domdec_t;

    dd->comm = init_dd_comm();

    /* Initialize DD paritioning counters */
    dd->comm->partition_step = INT_MIN;
    dd->ddp_count            = 0;

    set_dd_envvar_options(mdlog, dd, cr->nodeid);

    gmx_ddbox_t ddbox = {0};
    set_dd_limits_and_grid(mdlog, cr, dd, options, mdrunOptions,
                           mtop, ir,
                           box, xGlobal,
                           &ddbox);

    make_dd_communicators(mdlog, cr, dd, options.rankOrder);

    if (thisRankHasDuty(cr, DUTY_PP))
    {
        set_ddgrid_parameters(mdlog, dd, options.dlbScaling, mtop, ir, &ddbox);

        setup_neighbor_relations(dd);
    }

    /* Set overallocation to avoid frequent reallocation of arrays */
    set_over_alloc_dd(TRUE);

    clear_dd_cycle_counts(dd);

    dd->atomSets = atomSets;

    return dd;
}

static gmx_bool test_dd_cutoff(t_commrec     *cr,
                               const t_state &state,
                               real           cutoffRequested)
{
    gmx_domdec_t *dd;
    gmx_ddbox_t   ddbox;
    int           d, dim, np;
    real          inv_cell_size;
    int           LocallyLimited;

    dd = cr->dd;

    set_ddbox(*dd, false, state.box, true, state.x, &ddbox);

    LocallyLimited = 0;

    for (d = 0; d < dd->ndim; d++)
    {
        dim = dd->dim[d];

        inv_cell_size = DD_CELL_MARGIN*dd->nc[dim]/ddbox.box_size[dim];
        if (dynamic_dd_box(*dd))
        {
            inv_cell_size *= DD_PRES_SCALE_MARGIN;
        }

        np = 1 + static_cast<int>(cutoffRequested*inv_cell_size*ddbox.skew_fac[dim]);

        if (!isDlbDisabled(dd->comm) && (dim < ddbox.npbcdim) && (dd->comm->cd[d].np_dlb > 0))
        {
            if (np > dd->comm->cd[d].np_dlb)
            {
                return FALSE;
            }

            /* If a current local cell size is smaller than the requested
             * cut-off, we could still fix it, but this gets very complicated.
             * Without fixing here, we might actually need more checks.
             */
            real cellSizeAlongDim = (dd->comm->cell_x1[dim] - dd->comm->cell_x0[dim])*ddbox.skew_fac[dim];
            if (cellSizeAlongDim*dd->comm->cd[d].np_dlb < cutoffRequested)
            {
                LocallyLimited = 1;
            }
        }
    }

    if (!isDlbDisabled(dd->comm))
    {
        /* If DLB is not active yet, we don't need to check the grid jumps.
         * Actually we shouldn't, because then the grid jump data is not set.
         */
        if (isDlbOn(dd->comm) &&
            check_grid_jump(0, dd, cutoffRequested, &ddbox, FALSE))
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

gmx_bool change_dd_cutoff(t_commrec     *cr,
                          const t_state &state,
                          real           cutoffRequested)
{
    gmx_bool bCutoffAllowed;

    bCutoffAllowed = test_dd_cutoff(cr, state, cutoffRequested);

    if (bCutoffAllowed)
    {
        cr->dd->comm->cutoff = cutoffRequested;
    }

    return bCutoffAllowed;
}
